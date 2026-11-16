process FASTQC {
    tag "${sample}"
    label 'process_low'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample), val(batch), val(technology), path(r1), path(r2)

    output:
    path '*_fastqc.zip',  emit: zip
    path '*_fastqc.html', emit: html

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${r1} ${r2}
    """
}

process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path reports

    output:
    path 'multiqc_report.html'
    path 'multiqc_data'

    script:
    """
    multiqc --force .
    """
}

process KB_COUNT {
    tag "${sample}"
    label 'process_medium'
    publishDir "${params.outdir}/quant/${sample}", mode: 'copy', pattern: '*.json'

    input:
    tuple val(sample), val(batch), val(technology), path(r1), path(r2)
    tuple path(index), path(t2g)
    path whitelist

    output:
    tuple val(sample), path("${sample}.raw.h5ad"), path("${sample}.quant.json"), emit: counts
    path "${sample}.run_info.json", emit: log

    script:
    // lanes of one sample are passed as interleaved R1/R2 pairs
    def lanes_r1 = r1 instanceof List ? r1 : [r1]
    def lanes_r2 = r2 instanceof List ? r2 : [r2]
    def pairs = [lanes_r1, lanes_r2].transpose().collect { pair -> "${pair[0]} ${pair[1]}" }.join(' ')
    def onlist = whitelist ? "-w ${whitelist}" : ''
    // bustools sort allocates its whole -m buffer up front; asking for all of the task memory gets
    // the container killed at its memory limit, so give it half
    def mem = task.memory ? "-m ${Math.max(1, task.memory.toGiga().intdiv(2))}G" : ''
    """
    kb count -i ${index} -g ${t2g} -x ${technology} ${onlist} \\
        -t ${task.cpus} ${mem} -o kb_out ${pairs}
    test -s kb_out/counts_unfiltered/cells_x_genes.mtx || { echo "kb count produced no matrix" >&2; exit 1; }
    cp kb_out/run_info.json ${sample}.run_info.json
    scbf kb2h5ad --kb-dir kb_out --sample ${sample} --batch ${batch} \\
        --out ${sample}.raw.h5ad --stats ${sample}.quant.json
    """
}
