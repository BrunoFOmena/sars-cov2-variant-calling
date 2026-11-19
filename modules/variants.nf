process LOFREQ {
    tag "${sample}"
    publishDir { "${params.outdir}/variants/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path reference

    output:
    tuple val(sample), path("${sample}.lofreq.vcf.gz"), path("${sample}.lofreq.vcf.gz.tbi")

    script:
    """
    samtools faidx ${reference}
    lofreq indelqual --dindel -f ${reference} -o indelqual.bam ${bam}
    samtools index indelqual.bam
    lofreq call --call-indels -f ${reference} -o raw.vcf indelqual.bam
    bgzip -c raw.vcf > ${sample}.lofreq.vcf.gz
    tabix -p vcf ${sample}.lofreq.vcf.gz
    """
}

process FILTER_VARIANTS {
    tag "${sample}"
    publishDir { "${params.outdir}/variants/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(tbi)

    output:
    tuple val(sample), path("${sample}.filtered.vcf.gz"), path("${sample}.filtered.vcf.gz.tbi")

    script:
    """
    bcftools view -i 'QUAL>=${params.min_qual} && DP>=${params.min_depth} && INFO/AF>=${params.min_af}' \\
        ${vcf} -Oz -o ${sample}.filtered.vcf.gz
    tabix -p vcf ${sample}.filtered.vcf.gz
    """
}

process SNPEFF {
    tag "${sample}"
    publishDir { "${params.outdir}/variants/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(tbi)

    output:
    tuple val(sample), path("${sample}.annotated.vcf.gz"), emit: vcf
    path "${sample}.snpeff.csv", emit: stats
    path "${sample}.variants.tsv", emit: table

    script:
    """
    snpEff download -v ${params.snpeff_db} -dataDir \$PWD/snpeff_data
    snpEff ann -dataDir \$PWD/snpeff_data -csvStats ${sample}.snpeff.csv -stats ${sample}.snpeff.html ${params.snpeff_db} ${vcf} \\
        | bgzip -c > ${sample}.annotated.vcf.gz

    {
      printf 'sample\\tposition\\tref\\talt\\tdepth\\tallele_freq\\tqual\\teffect\\timpact\\tgene\\tprotein_change\\n'
      bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/DP\\t%INFO/AF\\t%QUAL\\t%INFO/ANN\\n' ${sample}.annotated.vcf.gz \\
        | awk -F'\\t' -v s=${sample} '{split(\$8,a,","); split(a[1],f,"|");
            printf "%s\\t%s\\t%s\\t%s\\t%s\\t%.3f\\t%s\\t%s\\t%s\\t%s\\t%s\\n", s,\$2,\$3,\$4,\$5,\$6,\$7,f[2],f[3],f[4],f[11]}'
    } > ${sample}.variants.tsv
    """
}

process CONSENSUS {
    tag "${sample}"
    publishDir { "${params.outdir}/consensus/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(tbi), path(per_base)
    path reference

    output:
    path "${sample}.consensus.fasta"
    path "${sample}.low_coverage.bed"

    script:
    """
    # Positions below min_depth are masked as N
    zcat ${per_base} | awk -v d=${params.min_depth} '\$4 < d' > ${sample}.low_coverage.bed

    # Only majority variants are applied to the consensus
    bcftools view -i 'INFO/AF>=${params.consensus_af}' ${vcf} -Oz -o consensus_vars.vcf.gz
    tabix -p vcf consensus_vars.vcf.gz

    bcftools consensus -f ${reference} -m ${sample}.low_coverage.bed -c ${sample}.chain consensus_vars.vcf.gz \\
        | sed "1s/.*/>${sample}/" > ${sample}.consensus.fasta
    """
}

// Independent allele counts from the raw pileup at every called SNV position.
// Uses different filters than LoFreq (no BAQ, MAPQ>=20, baseQ>=20) and reports
// per-strand counts so strand bias can be checked.
process PILEUP_CHECK {
    tag "${sample}"
    publishDir { "${params.outdir}/validation/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(tbi), path(bam), path(bai)
    path reference

    output:
    tuple val(sample), path("${sample}.pileup_counts.tsv")

    script:
    """
    bcftools query -f '%CHROM\\t%POS\\n' ${vcf} > positions.tsv
    if [ -s positions.tsv ]; then
      bcftools mpileup -f ${reference} -R positions.tsv -I -B -q 20 -Q 20 -d 100000 -a AD,ADF,ADR -Ou ${bam} \\
        | bcftools query -f '%POS\\t%REF\\t%ALT\\t[%AD]\\t[%ADF]\\t[%ADR]\\n' > ${sample}.pileup_counts.tsv
    else
      : > ${sample}.pileup_counts.tsv
    fi
    """
}
