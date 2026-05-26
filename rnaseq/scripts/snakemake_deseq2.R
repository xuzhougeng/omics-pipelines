#!/usr/bin/env Rscript
# DESeq2 differential expression for the RNA-seq pipeline.
# Called by Snakemake with:
#   Rscript snakemake_deseq2.R <counts_tsv> <samples_tsv> <out_dir> \
#                              <ref_group> <padj> <lfc>
#
# Grouping logic: the FIRST group in samples.tsv is the reference (WT). Every
# other group is contrasted against it (group2 vs group1, group3 vs group1, ...).

suppressMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
counts_tsv  <- args[1]
samples_tsv <- args[2]
out_dir     <- args[3]
ref_group   <- args[4]
padj_thr    <- as.numeric(args[5])
lfc_thr     <- as.numeric(args[6])

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load counts ----
counts <- read.delim(counts_tsv, row.names = 1, check.names = FALSE)

# ---- Load sample sheet ----
samples <- read.delim(samples_tsv, comment.char = "#", check.names = FALSE)
rownames(samples) <- samples$sample

# Keep only samples present in both, in sample-sheet order.
common <- samples$sample[samples$sample %in% colnames(counts)]
if (length(common) == 0) stop("No overlapping samples between counts and samples.tsv")
counts  <- counts[, common, drop = FALSE]
samples <- samples[common, , drop = FALSE]

# Group factor with the reference level (WT) first.
group_levels <- unique(samples$group)
if (!ref_group %in% group_levels) stop("ref_group not found in samples.tsv: ", ref_group)
group_levels <- c(ref_group, setdiff(group_levels, ref_group))
samples$group <- factor(samples$group, levels = group_levels)

coldata <- data.frame(group = samples$group, row.names = rownames(samples))

message(sprintf("DESeq2: %d genes x %d samples, reference group = '%s'",
                nrow(counts), ncol(counts), ref_group))

dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts)),
                              colData   = coldata,
                              design    = ~ group)

# Drop all-zero genes for cleaner dispersion estimation.
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- DESeq(dds)

# ---- Normalized counts (size-factor normalized) ----
norm <- counts(dds, normalized = TRUE)
norm_out <- file.path(out_dir, "normalized_counts.tsv")
write.table(data.frame(gene_id = rownames(norm), norm, check.names = FALSE),
            norm_out, sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote ", norm_out)

# ---- Per-contrast results: each non-reference group vs reference ----
contrast_groups <- setdiff(levels(samples$group), ref_group)
if (length(contrast_groups) == 0) {
  message("Only one group present; no contrasts to compute.")
}
for (g in contrast_groups) {
  res <- results(dds, contrast = c("group", g, ref_group))
  res <- res[order(res$padj), ]
  df <- data.frame(gene_id = rownames(res), as.data.frame(res), check.names = FALSE)
  df$significant <- !is.na(df$padj) & df$padj < padj_thr & abs(df$log2FoldChange) >= lfc_thr
  out_file <- file.path(out_dir, paste0(g, "_vs_", ref_group, ".deseq2.tsv"))
  write.table(df, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message(sprintf("Wrote %s (%d significant of %d)",
                  out_file, sum(df$significant), nrow(df)))
}
message("Done.")
