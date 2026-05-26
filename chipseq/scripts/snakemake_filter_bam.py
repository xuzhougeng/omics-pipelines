"""Filter a ChIP-seq BAM (no Tn5 shift).

Keeps reads that are: not duplicates, on a nuclear chromosome, above the MAPQ
threshold, and in a proper pair. Writes the filtered BAM and its index.

Usage: python snakemake_filter_bam.py <in_bam> <filtered_bam> <chroms> <mapq>
"""
import sys
import pysam

in_bam = sys.argv[1]
filtered_bam = sys.argv[2]
nuclear = set(sys.argv[3].split(","))
mapq_threshold = int(sys.argv[4])

STATS = {"total": 0, "dedup": 0, "nuclear": 0, "mapq": 0, "proper_pair": 0}

with pysam.AlignmentFile(in_bam, "rb", threads=2) as in_f:
    header = in_f.header
    with pysam.AlignmentFile(filtered_bam, "wb", header=header, threads=2) as filt_f:
        for read in in_f:
            STATS["total"] += 1
            if read.is_duplicate:
                continue
            STATS["dedup"] += 1
            if read.reference_name not in nuclear:
                continue
            STATS["nuclear"] += 1
            if read.mapping_quality <= mapq_threshold:
                continue
            STATS["mapq"] += 1
            if not read.is_proper_pair:
                continue
            STATS["proper_pair"] += 1
            filt_f.write(read)

pysam.index(filtered_bam)

for k, v in STATS.items():
    print(f"{k}\t{v}")
