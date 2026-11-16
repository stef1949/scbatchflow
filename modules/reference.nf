process SIMULATE {
    label 'process_low'

    output:
    path 'sim/genome.fa',       emit: genome
    path 'sim/genes.gtf',       emit: gtf
    path 'sim/whitelist.txt',   emit: whitelist
    path 'sim/samplesheet.csv', emit: samplesheet
    path 'sim/truth.csv',       emit: truth
    path 'sim/*.fastq.gz',      emit: fastq

    script:
    """
    scbf simulate --outdir sim --seed ${params.seed + 1}
    """
}

process KB_REF {
    label 'process_medium'
    publishDir "${params.outdir}/reference", mode: 'copy', pattern: 't2g.txt'

    input:
    path fasta
    path gtf

    output:
    tuple path('index.idx'), path('t2g.txt'), emit: index

    script:
    """
    kb ref -i index.idx -g t2g.txt -f1 cdna.fa ${fasta} ${gtf}
    """
}
