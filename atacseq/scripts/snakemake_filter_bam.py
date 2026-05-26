#!/usr/bin/env python3
"""Filter BAM and apply Tn5 shift.
Usage: python snakemake_filter_bam.py <in_bam> <filtered_bam> <shifted_bam> <chroms> <mapq>
"""
import os
import sys
import pysam

in_bam = sys.argv[1]
filtered_bam = sys.argv[2]
shifted_bam = sys.argv[3]
nuclear = set(sys.argv[4].split(","))
mapq_threshold = int(sys.argv[5])
shifted_unsorted = shifted_bam.replace(".bam", ".unsorted.bam")

STATS = {"total": 0, "dedup": 0, "nuclear": 0, "mapq": 0, "proper_pair": 0}

with pysam.AlignmentFile(in_bam, "rb", threads=2) as in_f:
    header = in_f.header
    with (
        pysam.AlignmentFile(filtered_bam, "wb", header=header, threads=2) as filt_f,
        pysam.AlignmentFile(shifted_unsorted, "wb", header=header, threads=2) as shift_f,
    ):
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
            # Tn5 shift: +4 for +strand, -5 for -strand
            if read.is_reverse:
                read.pos = max(0, read.pos - 5)
            else:
                read.pos = read.pos + 4
            shift_f.write(read)

pysam.index(filtered_bam)
pysam.sort("-o", shifted_bam, shifted_unsorted)
pysam.index(shifted_bam)
os.remove(shifted_unsorted)
