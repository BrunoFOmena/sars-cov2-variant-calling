#!/usr/bin/env bash
# Exploratory check of the raw reads: FastQC + MultiQC, and a quick BWA mapping
# against the reference to measure the viral read fraction and genome coverage.
# This is a sanity check to decide on the dataset; the full pipeline will be in Nextflow.
#
# Requires the conda env: micromamba activate sars2-vc
# Usage: bash scripts/qc_raw_and_mapping_check.sh [threads]
set -euo pipefail
THREADS="${1:-8}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF="$ROOT/data/reference/NC_045512.2.fasta"
OUT="$ROOT/results/exploratory"
mkdir -p "$OUT/fastqc" "$OUT/mapping"

fastqc -t "$THREADS" --outdir "$OUT/fastqc" "$ROOT"/data/raw/*/*.fastq.gz
multiqc -f -o "$OUT/multiqc" -n multiqc_raw "$OUT/fastqc"

[[ -f "$REF.bwt" ]] || bwa index "$REF"

summary="$OUT/mapping/mapping_summary.tsv"
printf 'sample\ttotal_reads\tmapped_reads\tmapped_pct\tmean_depth\tpct_genome_ge10x\tpct_genome_ge30x\n' > "$summary"
for dir in "$ROOT"/data/raw/*/; do
  s=$(basename "$dir")
  bam="$OUT/mapping/$s.sorted.bam"
  bwa mem -t "$THREADS" -R "@RG\tID:$s\tSM:$s\tPL:ILLUMINA" "$REF" "$dir/${s}_1.fastq.gz" "$dir/${s}_2.fastq.gz" 2> "$OUT/mapping/$s.bwa.log" \
    | samtools sort -@ 4 -o "$bam" -
  samtools index "$bam"
  samtools flagstat "$bam" > "$OUT/mapping/$s.flagstat.txt"
  total=$(samtools view -c -F 0x900 "$bam")
  mapped=$(samtools view -c -F 0x904 "$bam")
  samtools depth -a "$bam" | awk -v s="$s" -v t="$total" -v m="$mapped" \
    '{n++; sum+=$3; if($3>=10)a++; if($3>=30)b++} END{printf "%s\t%d\t%d\t%.2f\t%.1f\t%.1f\t%.1f\n", s,t,m,100*m/t,sum/n,100*a/n,100*b/n}' >> "$summary"
done
cat "$summary"
