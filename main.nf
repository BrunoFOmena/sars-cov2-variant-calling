#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW     } from './modules/qc'
include { FASTQC as FASTQC_TRIMMED } from './modules/qc'
include { FASTP; MULTIQC           } from './modules/qc'
include { BWA_INDEX; BWA_MEM; ALIGNMENT_STATS; MOSDEPTH } from './modules/align'
include { LOFREQ; FILTER_VARIANTS; SNPEFF; CONSENSUS; PILEUP_CHECK } from './modules/variants'
include { VALIDATE_AND_PLOT } from './modules/report'

workflow {
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
    reference = file(params.reference, checkIfExists: true)

    FASTQC_RAW(reads_ch, 'raw')
    FASTP(reads_ch)
    FASTQC_TRIMMED(FASTP.out.reads, 'trimmed')

    BWA_INDEX(reference)
    BWA_MEM(FASTP.out.reads, reference, BWA_INDEX.out)
    ALIGNMENT_STATS(BWA_MEM.out)
    MOSDEPTH(BWA_MEM.out)

    LOFREQ(BWA_MEM.out, reference)
    FILTER_VARIANTS(LOFREQ.out)
    SNPEFF(FILTER_VARIANTS.out)

    consensus_in = FILTER_VARIANTS.out
        .map { sample, vcf, tbi -> [sample, vcf, tbi] }
        .join(MOSDEPTH.out.per_base.map { f -> [f.name.replace('.per-base.bed.gz', ''), f] })
    CONSENSUS(consensus_in, reference)

    PILEUP_CHECK(FILTER_VARIANTS.out.join(BWA_MEM.out), reference)

    VALIDATE_AND_PLOT(
        SNPEFF.out.table.mix(
            PILEUP_CHECK.out.map { sample, f -> f },
            MOSDEPTH.out.per_base
        ).collect(),
        file(params.gff, checkIfExists: true)
    )

    MULTIQC(
        FASTQC_RAW.out.mix(FASTQC_TRIMMED.out)
            .mix(FASTP.out.json)
            .mix(ALIGNMENT_STATS.out.flagstat, ALIGNMENT_STATS.out.stats)
            .mix(MOSDEPTH.out.summary, MOSDEPTH.out.dist)
            .mix(SNPEFF.out.stats)
            .collect()
    )
}
