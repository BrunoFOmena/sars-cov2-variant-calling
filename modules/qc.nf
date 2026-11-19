process FASTQC {
    tag "${sample}"
    publishDir { "${params.outdir}/qc/fastqc_${stage}" }, mode: 'copy'

    input:
    tuple val(sample), path(reads)
    val stage

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    fastqc -t ${task.cpus} --outdir . ${reads}
    """
}

process FASTP {
    tag "${sample}"
    publishDir { "${params.outdir}/trimmed/${sample}" }, mode: 'copy', pattern: '*.{json,html}'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_trim_{1,2}.fastq.gz"), emit: reads
    path "${sample}.fastp.json", emit: json
    path "${sample}.fastp.html"

    script:
    """
    fastp \\
        -i ${reads[0]} -I ${reads[1]} \\
        -o ${sample}_trim_1.fastq.gz -O ${sample}_trim_2.fastq.gz \\
        --detect_adapter_for_pe \\
        --qualified_quality_phred 20 --cut_tail --cut_tail_mean_quality 20 \\
        --length_required 50 \\
        --thread ${task.cpus} \\
        --json ${sample}.fastp.json --html ${sample}.fastp.html
    """
}

process MULTIQC {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path files

    output:
    path "multiqc_report.html"
    path "multiqc_report_data"

    script:
    """
    multiqc --force --filename multiqc_report .
    """
}
