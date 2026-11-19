#!/usr/bin/env python3
"""Merge LoFreq calls with independent pileup counts, flag suspicious variants, and draw figures.

Inputs are discovered by file name in the working directory (Nextflow stages them):
  <sample>.variants.tsv          annotated LoFreq variants
  <sample>.pileup_counts.tsv     bcftools mpileup counts (pos, ref, alt, AD, ADF, ADR)
  <sample>.per-base.bed.gz       mosdepth per-base depth
  NC_045512.2.gff3               gene annotation
"""
import glob, gzip, os, sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

samples = sorted(f.replace(".variants.tsv", "") for f in glob.glob("*.variants.tsv"))
MIN_ALT_PILEUP = 3          # alt reads required in the independent count
STRAND_P = 0.001            # Fisher p-value below which strand imbalance is flagged

def parse_pileup(path, pos_set):
    rows = {}
    if os.path.getsize(path) == 0:
        return rows
    for line in open(path):
        pos, ref, alt, ad, adf, adr = line.rstrip("\n").split("\t")
        rows[int(pos)] = (alt.split(","), [int(x) for x in ad.split(",")],
                          [int(x) for x in adf.split(",")], [int(x) for x in adr.split(",")])
    return rows

tables = []
for s in samples:
    v = pd.read_csv(f"{s}.variants.tsv", sep="\t")
    pile = parse_pileup(f"{s}.pileup_counts.tsv", set(v.position))
    rec = []
    for _, r in v.iterrows():
        is_snv = len(r.ref) == 1 and len(r.alt) == 1
        out = dict(af_pileup=None, alt_fwd=None, alt_rev=None, ref_fwd=None, ref_rev=None, strand_p=None)
        if is_snv and r.position in pile:
            alts, ad, adf, adr = pile[r.position]
            if r.alt in alts:
                i = alts.index(r.alt) + 1
                tot = ad[0] + ad[i]
                out.update(af_pileup=ad[i] / tot if tot else None, alt_fwd=adf[i], alt_rev=adr[i],
                           ref_fwd=adf[0], ref_rev=adr[0])
                out["strand_p"] = fisher_exact([[adf[0], adr[0]], [adf[i], adr[i]]])[1]
        rec.append(out)
    v = pd.concat([v, pd.DataFrame(rec)], axis=1)

    def verdict(r):
        if len(r.ref) != len(r.alt):
            return "indel (not assessed)"
        if pd.isna(r.af_pileup):
            return "not found in pileup"
        if (r.alt_fwd + r.alt_rev) < MIN_ALT_PILEUP:
            return "unsupported"
        if r.strand_p is not None and r.strand_p < STRAND_P:
            return "strand bias"
        return "supported"
    v["validation"] = v.apply(verdict, axis=1)
    tables.append(v)

allv = pd.concat(tables, ignore_index=True)
allv.round(4).to_csv("variants_validated.tsv", sep="\t", index=False)

# ---- Sample comparison -------------------------------------------------------
key = allv.assign(variant=allv.position.astype(str) + ":" + allv.ref + ">" + allv.alt)
by_sample = {s: set(key[key["sample"] == s].variant) for s in samples}
with open("sample_comparison.tsv", "w") as fh:
    fh.write("variant\t" + "\t".join(samples) + "\n")
    for var in sorted(set().union(*by_sample.values()), key=lambda x: int(x.split(":")[0])):
        fh.write(var + "\t" + "\t".join("yes" if var in by_sample[s] else "no" for s in samples) + "\n")

# ---- Genes from GFF ------------------------------------------------------------
genes = []
for line in open("NC_045512.2.gff3"):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    if len(f) > 8 and f[2] == "gene":
        name = [a.split("=")[1] for a in f[8].split(";") if a.startswith("Name=")]
        genes.append((int(f[3]), int(f[4]), name[0] if name else ""))

colors = {"supported": "#1b9e77", "strand bias": "#d95f02", "unsupported": "#7570b3",
          "not found in pileup": "#7570b3", "indel (not assessed)": "#666666"}

# ---- Figure 1: coverage profile + variants per sample ------------------------
fig, axes = plt.subplots(len(samples), 1, figsize=(11, 3.6 * len(samples)), sharex=True)
axes = [axes] if len(samples) == 1 else axes
for ax, s in zip(axes, samples):
    d = pd.read_csv(f"{s}.per-base.bed.gz", sep="\t", header=None, names=["c", "start", "end", "depth"])
    xs, ys = [], []
    for _, r in d.iterrows():
        xs += [r.start, r.end]; ys += [r.depth, r.depth]
    ax.fill_between(xs, ys, color="#bbbbbb", step=None)
    ax.set_yscale("log"); ax.set_ylim(1, max(ys) * 3)
    ax.axhline(10, color="k", lw=0.7, ls="--")
    ax.set_ylabel("depth (log)"); ax.set_title(s, loc="left", fontsize=10)
    ax2 = ax.twinx()
    sv = allv[allv["sample"] == s]
    for _, r in sv.iterrows():
        ax2.vlines(r.position, 0, r.allele_freq, color=colors[r.validation], lw=1.6)
        ax2.plot(r.position, r.allele_freq, "o", color=colors[r.validation], ms=4)
    ax2.set_ylim(0, 1); ax2.set_ylabel("allele frequency")
    for i, (a, b, n) in enumerate(genes):
        ax.axvspan(a, b, ymin=0, ymax=0.04, color="#2c7fb8" if i % 2 else "#7fcdbb")
        if b - a > 1500: ax.text((a + b) / 2, 1.3, n, ha="center", fontsize=7)
axes[-1].set_xlabel("genome position (NC_045512.2)")
handles = [plt.Line2D([0], [0], marker="o", color=c, ls="", label=l) for l, c in colors.items()
           if l in set(allv.validation)]
axes[0].legend(handles=handles, loc="upper right", fontsize=7, title="variant validation", title_fontsize=7)
fig.tight_layout(); fig.savefig("coverage_and_variants.png", dpi=150); plt.close(fig)

# ---- Figure 2: LoFreq vs pileup allele frequency ------------------------------
fig, ax = plt.subplots(figsize=(4.6, 4.4))
snv = allv.dropna(subset=["af_pileup"])
for s, m in zip(samples, ["o", "s", "^", "D"]):
    t = snv[snv["sample"] == s]
    ax.scatter(t.allele_freq, t.af_pileup, marker=m, label=s, alpha=0.8)
ax.plot([0, 0.5], [0, 0.5], "k--", lw=0.7)
ax.set_xlabel("allele frequency (LoFreq)"); ax.set_ylabel("allele frequency (bcftools pileup)")
ax.legend(fontsize=7); fig.tight_layout(); fig.savefig("af_concordance.png", dpi=150); plt.close(fig)
print(allv.groupby(["sample", "validation"]).size())
