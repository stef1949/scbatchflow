#!/usr/bin/env nextflow

/*
 * scbatchflow
 * Raw 10x Chromium reads -> counts -> QC -> integration benchmark -> annotation -> PostgreSQL
 */

nextflow.enable.dsl = 2

include { SIMULATE; KB_REF } from './modules/reference'
include { FASTQC; MULTIQC; KB_COUNT } from './modules/quant'
include { CELL_QC; MERGE; INTEGRATE; METRICS; FINALISE; LOAD_DB } from './modules/analysis'

def checkParams() {
    if (!params.simulate) {
        if (!params.input) {
            error "Provide a samplesheet with --input (columns: sample,batch,fastq_1,fastq_2 and optionally technology)"
        }
        def prebuilt = params.kb_index && params.t2g
        def buildable = params.genome_fasta && params.gtf
        if (!prebuilt && !buildable) {
            error "Provide a reference: --kb_index and --t2g, or --genome_fasta and --gtf"
        }
    }
    if (params.methods.tokenize(',').size() < 1) {
        error "--methods needs at least one integration method"
    }
}

workflow {
    checkParams()

    // ---- inputs ----
    if (params.simulate) {
        SIMULATE()
        ch_sheet     = SIMULATE.out.samplesheet
        ch_fasta     = SIMULATE.out.genome
        ch_gtf       = SIMULATE.out.gtf
        ch_whitelist = SIMULATE.out.whitelist
    }
    else {
        ch_sheet     = Channel.fromPath(params.input, checkIfExists: true)
        ch_fasta     = params.genome_fasta ? Channel.value(file(params.genome_fasta, checkIfExists: true)) : Channel.empty()
        ch_gtf       = params.gtf ? Channel.value(file(params.gtf, checkIfExists: true)) : Channel.empty()
        ch_whitelist = params.whitelist ? Channel.value(file(params.whitelist, checkIfExists: true)) : Channel.value([])
    }

    // one entry per sample; several rows with the same sample are treated as lanes
    ch_reads = ch_sheet
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.batch ?: row.sample, row.technology ?: params.technology, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)) }
        .groupTuple(by: [0, 1, 2])

    // ---- reference ----
    if (!params.simulate && params.kb_index && params.t2g) {
        ch_index = Channel.value(tuple(file(params.kb_index, checkIfExists: true), file(params.t2g, checkIfExists: true)))
    }
    else {
        KB_REF(ch_fasta, ch_gtf)
        ch_index = KB_REF.out.index
    }

    // ---- read QC ----
    if (!params.skip_fastqc) {
        FASTQC(ch_reads)
        MULTIQC(FASTQC.out.zip.collect())
    }

    // ---- quantification and cell QC ----
    KB_COUNT(ch_reads, ch_index, ch_whitelist)
    CELL_QC(KB_COUNT.out.counts)

    // ---- merge, then benchmark every integration method in parallel ----
    ch_markers = Channel.value(file(params.markers, checkIfExists: true))
    MERGE(CELL_QC.out.h5ad.collect(sort: { a, b -> a.name <=> b.name }), ch_markers)

    ch_methods = Channel.fromList(params.methods.tokenize(',').collect { m -> m.trim() })
    INTEGRATE(ch_methods.combine(MERGE.out.h5ad))
    METRICS(INTEGRATE.out.embedding)

    def run_id = params.run_id ?: workflow.runName
    FINALISE(
        MERGE.out.h5ad,
        INTEGRATE.out.embedding.map { method, emb, merged -> emb }.collect(),
        METRICS.out.json.collect(),
        CELL_QC.out.summary.collect(),
        CELL_QC.out.plot.collect(),
        params.run_id ?: '',   // report title; a fixed value keeps -resume able to reuse this step
    )

    // ---- results database ----
    if (params.db_url) {
        def recorded = params.findAll { k, v -> !(k in ['db_url']) && !(v instanceof Map) }
        LOAD_DB(FINALISE.out.results, run_id, groovy.json.JsonOutput.toJson(recorded))
    }
}

workflow.onComplete {
    log.info(workflow.success ? "Done. Open ${params.outdir}/results/report.html" : "Pipeline failed.")
}
