process BWA_INDEX {
    storeDir "${projectDir}/data/reference/bwa_index"

    input:
    path reference

    output:
    path "${reference}.*"

    script:
    """
    bwa index ${reference}
    """
}

process BWA_MEM {
    tag "${sample}"
    publishDir { "${params.outdir}/alignment/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(reads)
    path reference
    path index

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")

    script:
    """
    bwa mem -t ${task.cpus} -R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA' ${reference} ${reads} \\
        | samtools sort -@ 2 -o ${sample}.bam -
    samtools index ${sample}.bam
    """
}

process ALIGNMENT_STATS {
    tag "${sample}"
    publishDir { "${params.outdir}/alignment/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path "${sample}.flagstat.txt", emit: flagstat
    path "${sample}.stats.txt", emit: stats

    script:
    """
    samtools flagstat ${bam} > ${sample}.flagstat.txt
    samtools stats ${bam} > ${sample}.stats.txt
    """
}

process MOSDEPTH {
    tag "${sample}"
    publishDir { "${params.outdir}/alignment/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path "${sample}.mosdepth.summary.txt", emit: summary
    path "${sample}.mosdepth.global.dist.txt", emit: dist
    path "${sample}.per-base.bed.gz", emit: per_base

    script:
    """
    mosdepth --threads ${task.cpus} ${sample} ${bam}
    """
}
