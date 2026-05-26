"""Merge per-sample STAR ReadsPerGene.out.tab files into one count matrix.

Usage:
    python merge_star_counts.py <strandedness> <out_tsv> <sample=path> [<sample=path> ...]

STAR --quantMode GeneCounts writes ReadsPerGene.out.tab with 4 columns:
    gene_id | unstranded | forward(yes) | reverse
and 4 leading summary rows (N_unmapped, N_multimapping, N_noFeature,
N_ambiguous) which are skipped here.
"""
import sys
import pandas as pd

strand = sys.argv[1]
out_tsv = sys.argv[2]
pairs = sys.argv[3:]

# strandedness -> column index in ReadsPerGene.out.tab (0-based)
STRAND_COL = {"unstranded": 1, "forward": 2, "reverse": 3}
if strand not in STRAND_COL:
    sys.exit(f"Unknown strandedness '{strand}'; expected one of {list(STRAND_COL)}")
col = STRAND_COL[strand]

merged = None
for pair in pairs:
    sample, path = pair.split("=", 1)
    df = pd.read_csv(
        path, sep="\t", header=None, skiprows=4,
        usecols=[0, col], names=["gene_id", sample],
    ).set_index("gene_id")
    merged = df if merged is None else merged.join(df, how="outer")

merged = merged.fillna(0).astype(int)
merged.to_csv(out_tsv, sep="\t")
print(f"Wrote {merged.shape[0]} genes x {merged.shape[1]} samples to {out_tsv}")
