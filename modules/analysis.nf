process CELL_QC {
    tag "${sample}"
    label 'process_low'
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.json'

    input:
    tuple val(sample), path(raw), path(stats)

    output:
    path "${sample}.qc.h5ad",      emit: h5ad
    path "${sample}.summary.json", emit: summary
    path "${sample}.qc.png",       emit: plot

    script:
    def doublets = params.skip_doublets ? '--no-doublets' : ''
    """
    scbf qc --input ${raw} --stats ${stats} \\
        --out ${sample}.qc.h5ad --summary ${sample}.summary.json \\
        --expected-cells ${params.expected_cells} --min-genes ${params.min_genes} \\
        --max-pct-mito ${params.max_pct_mito} --mito-mads ${params.mito_mads} \\
        ${doublets} --plot ${sample}.qc.png --seed ${params.seed}
    """
}

process MERGE {
    label 'process_medium'

    input:
    path h5ads
    path markers

    output:
    path 'merged.h5ad', emit: h5ad

    script:
    """
    scbf merge --inputs ${h5ads} --markers ${markers} --out merged.h5ad \\
        --n-hvg ${params.n_hvg} --n-pcs ${params.n_pcs} --seed ${params.seed}
    """
}

process INTEGRATE {
    tag "${method}"
    label 'process_medium'

    input:
    tuple val(method), path(merged)

    output:
    tuple val(method), path("${method}.embedding.npy"), path(merged), emit: embedding

    script:
    """
    scbf integrate --input ${merged} --method ${method} --batch-key ${params.batch_key} \\
        --out ${method}.embedding.npy --seed ${params.seed}
    """
}

process METRICS {
    tag "${method}"
    label 'process_medium'

    input:
    tuple val(method), path(embedding), path(merged)

    output:
    path "${method}.metrics.json", emit: json

    script:
    """
    scbf metrics --input ${merged} --embedding ${embedding} --method ${method} \\
        --batch-key ${params.batch_key} --out ${method}.metrics.json --seed ${params.seed}
    """
}

process FINALISE {
    label 'process_medium'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path merged
    path embeddings
    path metrics
    path summaries
    path qc_plots
    val run_label

    output:
    path 'results', emit: results

    script:
    """
    scbf finalise --input ${merged} --embeddings ${embeddings} --metrics ${metrics} \\
        --sample-summaries ${summaries} --qc-plots ${qc_plots} \\
        --batch-key ${params.batch_key} --run-label '${run_label}' \\
        --final-method ${params.final_method} \\
        --resolution ${params.resolution} --outdir results --seed ${params.seed}
    """
}

process LOAD_DB {
    label 'process_low'
    // the database is shared state: never run two loads at once, never cache one
    maxForks 1
    cache false

    input:
    path results
    val run_id
    val params_json

    output:
    path 'db_load.json'

    script:
    """
    cat > run_params.json <<'PARAMS'
    ${params_json}
    PARAMS
    export SCBF_DB_URL='${params.db_url}'
    scbf load-db --results ${results} --run-id '${run_id}' \\
        --params-file run_params.json > db_load.json
    """
}
