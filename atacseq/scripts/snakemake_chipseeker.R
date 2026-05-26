#!/usr/bin/env Rscript
# ChIPseeker annotation of consensus peaks
# Called by Snakemake with:
#   Rscript this_script.R <bed_file> <gtf_file> <out_tsv> \
#                         <organism> <anno_db> <tss_low> <tss_high>

library(ChIPseeker)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)

args <- commandArgs(trailingOnly = TRUE)
bed_file <- args[1]
gtf_file <- args[2]
out_tsv  <- args[3]
organism <- args[4]
anno_db  <- args[5]
tss_low  <- as.integer(args[6])
tss_high <- as.integer(args[7])

# Annotation DB is configurable (e.g. org.Mm.eg.db for mouse) and loaded
# at runtime rather than hard-imported.
library(anno_db, character.only = TRUE)

txdb_cache <- sub("\\.gtf$", "_txdb.sqlite", gtf_file)

message(sprintf("Loading peaks from %s", bed_file))
peaks <- import.bed(bed_file)
message(sprintf("Loaded %d peaks", length(peaks)))

if (is.null(peaks$name)) {
  peaks$name <- paste0("peak_", seq_len(length(peaks)))
}

# Peaks already come from BAMs restricted to the nuclear chromosomes in
# config.yaml, so just normalise the seqlevels to those actually present.
seqlevels(peaks) <- seqlevelsInUse(peaks)
message(sprintf("Annotating %d peaks across %d chromosomes",
                length(peaks), length(seqlevels(peaks))))

if (file.exists(txdb_cache)) {
  message(sprintf("Loading cached TxDb: %s", txdb_cache))
  txdb <- AnnotationDbi::loadDb(txdb_cache)
} else {
  message(sprintf("Building TxDb from %s (one-time, slow)", gtf_file))
  suppressMessages(
    txdb <- txdbmaker::makeTxDbFromGFF(gtf_file, format = "gtf", organism = organism)
  )
  AnnotationDbi::saveDb(txdb, txdb_cache)
  message(sprintf("TxDb cached: %s", txdb_cache))
}

message("Running annotatePeak...")
peak_anno <- annotatePeak(
  peaks,
  TxDb      = txdb,
  annoDb    = anno_db,
  tssRegion = c(tss_low, tss_high),
  genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR",
                                "Exon", "Intron", "Downstream", "Intergenic")
)

anno_df <- as.data.frame(peak_anno)

output <- data.frame(
  peak_id        = peaks$name,
  chrom          = seqnames(peaks),
  start          = start(peaks) - 1L,
  end            = end(peaks),
  annotation     = anno_df$annotation,
  distanceToTSS  = anno_df$distanceToTSS,
  geneId         = anno_df$geneId,
  geneChr        = anno_df$geneChr,
  geneStart      = anno_df$geneStart,
  geneEnd        = anno_df$geneEnd,
  geneStrand     = anno_df$geneStrand,
  geneLength     = anno_df$geneLength,
  transcriptId   = anno_df$transcriptId,
  stringsAsFactors = FALSE
)

gene_ids <- as.character(output$geneId)
gene_ids <- gsub("\\.\\d+$", "", gene_ids)
anno_db_obj <- get(anno_db, envir = asNamespace(anno_db))
symbols <- mapIds(anno_db_obj, keys = gene_ids, column = "SYMBOL",
                  keytype = "ENSEMBL", multiVals = "first")
output$geneSymbol <- symbols

message(sprintf("Writing %d annotated peaks to %s", nrow(output), out_tsv))
write.table(output, file = out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
message("Done.")
