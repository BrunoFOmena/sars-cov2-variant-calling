#!/usr/bin/env bash
# Download raw reads (ENA) and the reference genome (NCBI) and record provenance.
#
# Outputs:
#   data/raw/<SRR>/<SRR>_{1,2}.fastq.gz   raw paired-end reads (git-ignored)
#   data/reference/NC_045512.2.{fasta,gff3}  reference genome and annotation (tracked)
#   data/MANIFEST.tsv                      source URL, size, md5/sha256, date (tracked)
#
# Usage: bash scripts/download_data.sh
set -euo pipefail

ACCESSIONS=(SRR10903401 SRR10903402)
REF_ACC="NC_045512.2"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW="$ROOT/data/raw"
REF="$ROOT/data/reference"
MANIFEST="$ROOT/data/MANIFEST.tsv"
mkdir -p "$RAW" "$REF"

printf 'file\tsource_url\tbytes\tchecksum_type\tchecksum\tretrieved_utc\n' > "$MANIFEST"
# Download to a temp file and rename only on success, so partial files never look complete.
fetch() {
  [[ -s "$2" ]] && return 0
  curl -fsS --retry 5 --retry-all-errors -o "$2.part" "$1" && mv "$2.part" "$2"
}
now() { date -u +%Y-%m-%dT%H:%M:%SZ; }

# ---- Reads: ENA provides the md5 for each file, so we verify against it ----
for acc in "${ACCESSIONS[@]}"; do
  mkdir -p "$RAW/$acc"
  report=$(curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" | tail -n1)
  IFS=$'\t' read -r _run ftp_list md5_list <<< "$report"
  IFS=';' read -ra urls <<< "$ftp_list"
  IFS=';' read -ra md5s <<< "$md5_list"
  for i in "${!urls[@]}"; do
    url="https://${urls[$i]}"
    out="$RAW/$acc/$(basename "${urls[$i]}")"
    fetch "$url" "$out"
    got=$(md5sum "$out" | cut -d' ' -f1)
    if [[ "$got" != "${md5s[$i]}" ]]; then
      echo "MD5 mismatch for $out (expected ${md5s[$i]}, got $got)" >&2; exit 1
    fi
    printf '%s\t%s\t%s\tmd5\t%s\t%s\n' "${out#$ROOT/}" "$url" "$(stat -c %s "$out")" "$got" "$(now)" >> "$MANIFEST"
  done
done

# ---- Reference: NCBI E-utilities (no server-side checksum, so we record sha256) ----
EFETCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${REF_ACC}"
for spec in "fasta:fasta" "gff3:gff3"; do
  rettype="${spec%%:*}"; ext="${spec##*:}"
  out="$REF/${REF_ACC}.${ext}"
  url="${EFETCH}&rettype=${rettype}&retmode=text"
  fetch "$url" "$out"
  printf '%s\t%s\t%s\tsha256\t%s\t%s\n' "${out#$ROOT/}" "$url" "$(stat -c %s "$out")" "$(sha256sum "$out" | cut -d' ' -f1)" "$(now)" >> "$MANIFEST"
done

echo "Done. Provenance written to ${MANIFEST#$ROOT/}"
