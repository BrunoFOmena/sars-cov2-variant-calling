process VALIDATE_AND_PLOT {
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path files
    path gff

    output:
    path "variants_validated.tsv"
    path "sample_comparison.tsv"
    path "*.png"

    script:
    """
    validate_and_plot.py
    """
}
