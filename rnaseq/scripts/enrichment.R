#!/usr/bin/env Rscript
# GO / KEGG over-representation enrichment for the RNA-seq pipeline (step 06).
#
# This step intentionally runs on the *system* R (not a pixi environment): the
# clusterProfiler + organism annotation stack is large and version-sensitive,
# so the user installs it themselves (see README). The script fails fast with a
# clear message if a required package is missing.
#
# It does NOT re-run DESeq2. It consumes one contrast's DESeq2 result table
# (results/05-deseq2/{contrast}.deseq2.tsv) produced by snakemake_deseq2.R,
# splits the significant genes into up-/down-regulated sets using the existing
# `significant` flag and the sign of log2FoldChange, and runs enrichGO (BP/CC/MF)
# and enrichKEGG against the full set of tested genes as the background universe.
#
# Only human, mouse and arabidopsis are supported.

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript --vanilla scripts/enrichment.R",
      "    --deseq2 <results/05-deseq2/{contrast}.deseq2.tsv>",
      "    --output-dir <results/06-enrichment/{contrast}>",
      "    --organism <human|mouse|arabidopsis>",
      "    --gene-id-type <ENSEMBL|SYMBOL|TAIR|ENTREZID>",
      "    --strip-version <true|false>",
      "    --run-kegg <true|false>",
      "    --pvalue-cutoff <0.05>",
      "    --qvalue-cutoff <0.2>",
      sep = "\n"
    ),
    "\n"
  )
}

parse_args <- function(args) {
  if ("--help" %in% args || "-h" %in% args) {
    usage()
    quit(save = "no", status = 0)
  }

  expected <- c(
    "--deseq2", "--output-dir", "--organism", "--gene-id-type",
    "--strip-version", "--run-kegg", "--pvalue-cutoff", "--qvalue-cutoff"
  )
  values <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!(key %in% expected)) {
      stop(sprintf("Unknown argument: %s", key))
    }
    if (i == length(args)) {
      stop(sprintf("Missing value for argument: %s", key))
    }
    values[[sub("^--", "", key)]] <- args[[i + 1]]
    i <- i + 2
  }

  missing <- setdiff(sub("^--", "", expected), names(values))
  if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")))
  }
  values
}

parse_logical <- function(value, name) {
  normalized <- tolower(trimws(value))
  if (normalized %in% c("true", "t", "yes", "1")) return(TRUE)
  if (normalized %in% c("false", "f", "no", "0")) return(FALSE)
  stop(sprintf("Expected a boolean for %s, got: %s", name, value))
}

log_message <- function(...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, paste(..., collapse = "")))
}

# --- Per-organism annotation map -------------------------------------------
# org_db        : Bioconductor OrgDb package providing the annotations.
# kegg_organism : KEGG 3-letter organism code.
# kegg_keytype  : the OrgDb column whose values are valid KEGG gene IDs, i.e.
#                 the ID type enrichKEGG expects for this organism. Human/mouse
#                 KEGG uses NCBI Entrez gene IDs; arabidopsis uses TAIR loci.
ORGANISMS <- list(
  human       = list(org_db = "org.Hs.eg.db", kegg_organism = "hsa", kegg_keytype = "ENTREZID"),
  mouse       = list(org_db = "org.Mm.eg.db", kegg_organism = "mmu", kegg_keytype = "ENTREZID"),
  arabidopsis = list(org_db = "org.At.tair.db", kegg_organism = "ath", kegg_keytype = "TAIR")
)

args <- parse_args(commandArgs(trailingOnly = TRUE))

organism <- tolower(trimws(args$organism))
if (!organism %in% names(ORGANISMS)) {
  stop(sprintf(
    "Unsupported organism: '%s'. Supported: %s",
    args$organism, paste(names(ORGANISMS), collapse = ", ")
  ))
}
org_info <- ORGANISMS[[organism]]
org_db_name <- org_info$org_db

gene_id_type <- toupper(trimws(args$`gene-id-type`))
strip_version <- parse_logical(args$`strip-version`, "--strip-version")
run_kegg_flag <- parse_logical(args$`run-kegg`, "--run-kegg")
pvalue_cutoff <- as.numeric(args$`pvalue-cutoff`)
qvalue_cutoff <- as.numeric(args$`qvalue-cutoff`)
if (is.na(pvalue_cutoff) || is.na(qvalue_cutoff)) {
  stop("--pvalue-cutoff and --qvalue-cutoff must be numeric")
}

# --- Fail fast on missing packages (user installs these on system R) --------
required_packages <- c("clusterProfiler", "AnnotationDbi", "ggplot2", org_db_name)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(sprintf(
    paste0(
      "Missing R packages on the system R interpreter: %s\n",
      "Install them once with, e.g.:\n",
      "  Rscript -e 'if (!requireNamespace(\"BiocManager\", quietly=TRUE)) install.packages(\"BiocManager\"); ",
      "BiocManager::install(c(%s))'"
    ),
    paste(missing_packages, collapse = ", "),
    paste(sprintf('"%s"', missing_packages), collapse = ", ")
  ))
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(AnnotationDbi)
  library(ggplot2)
  library(org_db_name, character.only = TRUE)
})
org_db <- get(org_db_name)

deseq2_path <- normalizePath(args$deseq2, mustWork = TRUE)
output_dir <- normalizePath(args$`output-dir`, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Helpers ----------------------------------------------------------------
empty_enrichment_table <- function() {
  data.frame(
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer(),
    stringsAsFactors = FALSE
  )
}

write_tsv <- function(data, path, row_names = FALSE) {
  write.table(data, file = path, sep = "\t", quote = FALSE, row.names = row_names, col.names = !row_names)
}

gene_ratio_to_numeric <- function(values) {
  vapply(
    strsplit(values, "/", fixed = TRUE),
    function(x) {
      if (length(x) != 2) {
        return(NA_real_)
      }
      as.numeric(x[[1]]) / as.numeric(x[[2]])
    },
    numeric(1)
  )
}

# Strip Ensembl version suffixes (e.g. ENSG00000223972.5 -> ENSG00000223972).
clean_ids <- function(ids) {
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (strip_version) {
    ids <- sub("\\.[0-9]+$", "", ids)
  }
  unique(ids)
}

# Convert IDs from gene_id_type to the ID type KEGG expects for this organism.
to_kegg_ids <- function(ids) {
  if (length(ids) == 0) {
    return(character(0))
  }
  if (identical(gene_id_type, org_info$kegg_keytype)) {
    return(unique(ids))
  }
  mapped <- tryCatch(
    suppressWarnings(suppressMessages(
      clusterProfiler::bitr(
        ids,
        fromType = gene_id_type,
        toType = org_info$kegg_keytype,
        OrgDb = org_db
      )
    )),
    error = function(e) NULL
  )
  if (is.null(mapped) || nrow(mapped) == 0) {
    return(character(0))
  }
  unique(mapped[[org_info$kegg_keytype]])
}

save_enrichment_plots <- function(result_df, title, output_prefix) {
  if (nrow(result_df) == 0) {
    return(invisible(NULL))
  }

  plot_df <- result_df[order(result_df$p.adjust, result_df$pvalue), , drop = FALSE]
  plot_df <- head(plot_df, 20)
  plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))
  plot_df$neg_log10_padj <- -log10(plot_df$p.adjust)
  plot_df$gene_ratio_numeric <- gene_ratio_to_numeric(plot_df$GeneRatio)

  bar_plot <- ggplot(plot_df, aes(x = Description, y = neg_log10_padj)) +
    geom_col(fill = "#3B7A57") +
    coord_flip() +
    labs(title = title, x = NULL, y = expression(-log[10](adjusted~p))) +
    theme_bw(base_size = 11)

  dot_plot <- ggplot(
    plot_df,
    aes(x = gene_ratio_numeric, y = Description, size = Count, color = neg_log10_padj)
  ) +
    geom_point() +
    labs(
      title = title,
      x = "Gene Ratio",
      y = NULL,
      color = expression(-log[10](adjusted~p))
    ) +
    theme_bw(base_size = 11)

  ggsave(paste0(output_prefix, "_barplot.pdf"), bar_plot, width = 9, height = 6)
  ggsave(paste0(output_prefix, "_dotplot.pdf"), dot_plot, width = 9, height = 6)
}

run_go <- function(genes, ontology, universe, out_dir, label, summary_lines) {
  result_path <- file.path(out_dir, sprintf("%s.tsv", label))

  if (length(genes) == 0) {
    write_tsv(empty_enrichment_table(), result_path)
    summary_lines[[length(summary_lines) + 1]] <- sprintf("%s terms: 0 (empty gene set)", label)
    return(summary_lines)
  }

  ego <- tryCatch(
    enrichGO(
      gene = genes,
      universe = universe,
      OrgDb = org_db,
      keyType = gene_id_type,
      ont = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff,
      readable = FALSE
    ),
    error = function(e) e
  )

  if (inherits(ego, "error")) {
    write_tsv(empty_enrichment_table(), result_path)
    summary_lines[[length(summary_lines) + 1]] <- sprintf("%s failed: %s", label, conditionMessage(ego))
    return(summary_lines)
  }

  result_df <- as.data.frame(ego)
  if (nrow(result_df) == 0) {
    write_tsv(empty_enrichment_table(), result_path)
    summary_lines[[length(summary_lines) + 1]] <- sprintf("%s terms: 0", label)
    return(summary_lines)
  }

  write_tsv(result_df, result_path)
  save_enrichment_plots(result_df, label, file.path(out_dir, label))
  summary_lines[[length(summary_lines) + 1]] <- sprintf("%s terms: %d", label, nrow(result_df))
  summary_lines
}

run_kegg <- function(genes, universe, out_dir, label, summary_lines) {
  result_path <- file.path(out_dir, sprintf("%s.tsv", label))

  kegg_genes <- to_kegg_ids(genes)
  kegg_universe <- to_kegg_ids(universe)
  summary_lines[[length(summary_lines) + 1]] <- sprintf(
    "%s mapped genes: %d/%d", label, length(kegg_genes), length(genes)
  )

  if (length(kegg_genes) == 0) {
    write_tsv(empty_enrichment_table(), result_path)
    summary_lines[[length(summary_lines) + 1]] <- sprintf("%s terms: 0 (no KEGG-mapped genes)", label)
    return(summary_lines)
  }

  ekegg <- tryCatch(
    enrichKEGG(
      gene = kegg_genes,
      universe = if (length(kegg_universe) > 0) kegg_universe else NULL,
      organism = org_info$kegg_organism,
      keyType = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff
    ),
    error = function(e) e
  )

  if (inherits(ekegg, "error")) {
    write_tsv(empty_enrichment_table(), result_path)
    summary_lines[[length(summary_lines) + 1]] <- sprintf("%s failed: %s", label, conditionMessage(ekegg))
    return(summary_lines)
  }

  result_df <- as.data.frame(ekegg)
  if (nrow(result_df) == 0) {
    write_tsv(empty_enrichment_table(), result_path)
    summary_lines[[length(summary_lines) + 1]] <- sprintf("%s terms: 0", label)
    return(summary_lines)
  }

  write_tsv(result_df, result_path)
  save_enrichment_plots(result_df, label, file.path(out_dir, label))
  summary_lines[[length(summary_lines) + 1]] <- sprintf("%s terms: %d", label, nrow(result_df))
  summary_lines
}

# --- Load the DESeq2 result table -------------------------------------------
log_message("Reading DESeq2 results: ", deseq2_path)
res <- read.delim(deseq2_path, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("gene_id", "log2FoldChange", "significant")
missing_cols <- setdiff(required_cols, colnames(res))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "DESeq2 table is missing expected column(s): %s",
    paste(missing_cols, collapse = ", ")
  ))
}

sig <- as.logical(res$significant)
sig[is.na(sig)] <- FALSE
lfc <- suppressWarnings(as.numeric(res$log2FoldChange))

up_genes <- clean_ids(res$gene_id[sig & !is.na(lfc) & lfc > 0])
down_genes <- clean_ids(res$gene_id[sig & !is.na(lfc) & lfc < 0])
universe_genes <- clean_ids(res$gene_id)

log_message(sprintf(
  "Universe: %d genes | up: %d | down: %d (organism=%s, keytype=%s)",
  length(universe_genes), length(up_genes), length(down_genes), organism, gene_id_type
))

summary_lines <- c(
  sprintf("Organism: %s (OrgDb=%s, KEGG=%s)", organism, org_db_name, org_info$kegg_organism),
  sprintf("Gene ID type: %s (strip_version=%s)", gene_id_type, strip_version),
  sprintf("Background universe genes: %d", length(universe_genes)),
  sprintf("Significant up genes: %d", length(up_genes)),
  sprintf("Significant down genes: %d", length(down_genes))
)

# --- GO over-representation (BP/CC/MF, up & down) ---------------------------
summary_lines <- run_go(up_genes, "BP", universe_genes, output_dir, "go_bp_up", summary_lines)
summary_lines <- run_go(up_genes, "CC", universe_genes, output_dir, "go_cc_up", summary_lines)
summary_lines <- run_go(up_genes, "MF", universe_genes, output_dir, "go_mf_up", summary_lines)
summary_lines <- run_go(down_genes, "BP", universe_genes, output_dir, "go_bp_down", summary_lines)
summary_lines <- run_go(down_genes, "CC", universe_genes, output_dir, "go_cc_down", summary_lines)
summary_lines <- run_go(down_genes, "MF", universe_genes, output_dir, "go_mf_down", summary_lines)

# --- KEGG over-representation (up & down) -----------------------------------
if (run_kegg_flag) {
  summary_lines <- run_kegg(up_genes, universe_genes, output_dir, "kegg_up", summary_lines)
  summary_lines <- run_kegg(down_genes, universe_genes, output_dir, "kegg_down", summary_lines)
} else {
  summary_lines[[length(summary_lines) + 1]] <- "KEGG: skipped (run_kegg=false)"
}

writeLines(summary_lines, file.path(output_dir, "summary.txt"))
log_message("Enrichment finished successfully")
