# SARS-CoV-2 Variant Calling Pipeline

A reproducible, laptop-friendly pipeline that goes from raw paired-end Illumina reads to an annotated list of
SARS-CoV-2 variants, a consensus genome and a lineage assignment.

> **Status:** planning phase. This README documents the strategy; results and figures will be added as each
> stage is implemented.

## Goals

1. Build an end-to-end variant calling workflow with standard, widely used bioinformatics tools.
2. Keep it fully reproducible (pinned environment, workflow manager, scripted downloads).
3. Run on a consumer laptop in minutes, not hours.
4. Produce portfolio-ready outputs: QC reports, coverage plots, an annotated variant table and a consensus genome.

## Data

| Item | Accession | Details |
|---|---|---|
| Sample 1 | [SRR10903401](https://www.ncbi.nlm.nih.gov/sra/SRR10903401) | Illumina MiSeq, paired-end, 476,632 read pairs, ~144 Mbp, ~31 MB FASTQ.gz |
| Sample 2 | [SRR10903402](https://www.ncbi.nlm.nih.gov/sra/SRR10903402) | Illumina MiSeq, paired-end, 676,694 read pairs, ~204 Mbp, ~47 MB FASTQ.gz |
| Reference | [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) | Wuhan-Hu-1, complete genome, 29,903 bp (accession verified against NCBI) |

The reads come from the study published in *National Science Review* ([doi:10.1093/nsr/nwaa036](https://doi.org/10.1093/nsr/nwaa036))
and are also the example dataset of the [V-pipe](https://cbg-ethz.github.io/V-pipe/) tutorial.

**Caveat:** the library strategy is RNA-Seq (metatranscriptomics) from a human clinical sample, so a large share
of reads may be human/host. The mapping rate and the depth across the genome are therefore a QC result to be
reported, not an assumption. Any region with low coverage will be masked in the consensus.

## Compute budget

Measured on the development machine (WSL2 on Windows):

| Resource | Available | Needed by this project |
|---|---|---|
| CPU | Intel i5-12450HX, 12 threads | 4-8 threads is plenty |
| RAM | 7.6 GB (+2 GB swap) | < 2 GB (a 30 kb genome index is tiny) |
| Disk | ~940 GB free in WSL | < 1 GB including raw data and intermediates |
| GPU | RTX 3050 6 GB | Not used (no step benefits from it) |
| Network | OK | ~80 MB of downloads |

**Conclusion:** the whole dataset fits comfortably on this laptop. No cloud or HPC resources are required, and
the full pipeline should finish in a few minutes. Since compute is not a constraint, the effort goes into
analysis quality (QC, filtering, annotation, interpretation) rather than into scaling.

## Strategy

```
SRA (ENA mirror) ──► raw FASTQ (R1/R2)
                         │
                         ▼
              FastQC / MultiQC (raw)
                         │
                         ▼
          fastp: paired-end trimming + filtering
                         │
                         ▼
              FastQC / MultiQC (trimmed)
                         │
                         ▼
        BWA-MEM2 (or bwa mem) ──► samtools sort/index
                         │
                         ▼
     samtools stats / flagstat / depth, mosdepth
                         │
                         ▼
        LoFreq variant calling (low-frequency aware)
                         │
                         ▼
    bcftools filter ──► SnpEff annotation (NC_045512.2)
                         │
                    ┌────┴─────┐
                    ▼          ▼
        bcftools consensus   Pangolin / Nextclade
        (mask depth < 10x)   (lineage / clade)
                         │
                         ▼
        Summary tables + plots (Python / R)
```

### Stages

1. **Data retrieval.** Download FASTQ from the ENA mirror of the SRA accessions (plain HTTPS, no `sra-tools`
   needed) and the reference FASTA/GFF from NCBI. Checksums are recorded.
2. **Raw QC.** FastQC per file, aggregated with MultiQC.
3. **Trimming and filtering.** `fastp` in paired-end mode: adapter auto-detection, quality trimming, length
   filter, pair-aware. Its JSON/HTML report is kept.
4. **Alignment.** Map to NC_045512.2, sort and index. Duplicates are *not* blindly removed, since amplicon-like or
   very high coverage in a tiny genome makes position-based duplicate removal misleading; this decision will be
   justified in the report.
5. **Alignment QC.** Mapping rate, mean depth, % genome at >=10x/>=30x, depth profile along the genome.
6. **Variant calling.** `LoFreq` for SNVs and indels with allele frequencies (suited to viral populations where
   minor variants matter). Optionally cross-checked with `bcftools mpileup | call` as an independent caller to
   report concordance.
7. **Filtering and annotation.** Filter on depth, quality and allele frequency; annotate effects on genes with
   SnpEff using the Wuhan-Hu-1 annotation (e.g. S, ORF1ab, N).
8. **Consensus and lineage.** Build a consensus genome with low-coverage positions masked as `N`, then assign
   lineage/clade with Pangolin and Nextclade.
9. **Reporting.** Per-sample variant tables, allele frequency distribution, genome-wide coverage plot and a final
   summary table.

### Tools

| Purpose | Tool |
|---|---|
| Workflow manager | Nextflow (DSL2) |
| Environment | Conda (via micromamba, already installed) with bioconda/conda-forge, pinned in `environment.yml` |
| QC | FastQC, MultiQC, fastp |
| Alignment | BWA-MEM2 / BWA, samtools, mosdepth |
| Variants | LoFreq, bcftools |
| Annotation | SnpEff |
| Lineage | Pangolin, Nextclade |
| Plots | Python (matplotlib/pandas) |

### Planned repository layout

```
.
├── README.md
├── environment.yml        # pinned conda environment
├── main.nf                # Nextflow workflow definition
├── nextflow.config        # accessions, paths, thresholds, profiles
├── modules/               # one process module per tool
├── bin/                   # plotting and summary helper scripts
├── data/                  # raw reads and reference (git-ignored)
├── results/               # pipeline outputs (large files git-ignored)
└── docs/                  # figures and report used in the README
```

## Expected deliverables (for the portfolio)

- MultiQC report before and after trimming.
- Genome coverage/depth plot for each sample.
- Annotated VCF and a readable variant table (position, ref/alt, gene, amino-acid change, allele frequency).
- Consensus FASTA and lineage/clade call.
- Comparison of the two samples (shared vs. private variants).
- A short interpretation section, including limitations (host contamination, low-frequency variant calls,
  no biological replicates).

## Roadmap

- [x] Assess hardware and dataset size
- [x] Verify accessions (SRR10903401, SRR10903402, NC_045512.2)
- [ ] Set up the conda environment
- [ ] Data download script
- [ ] QC, trimming, alignment
- [ ] Variant calling and annotation
- [ ] Consensus and lineage
- [ ] Figures and results section in this README

## Reproducibility

Once implemented, the whole analysis will be run with:

```bash
micromamba env create -f environment.yml
micromamba activate sars2-vc
nextflow run main.nf -profile conda
```

## References

- Dataset publication: [doi:10.1093/nsr/nwaa036](https://doi.org/10.1093/nsr/nwaa036)
- Reference genome: [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2)
- V-pipe tutorial dataset: <https://cbg-ethz.github.io/V-pipe/tutorial/sars-cov2/>
