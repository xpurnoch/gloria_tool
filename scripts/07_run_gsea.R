#!/usr/bin/env Rscript
# =============================================================================
# 07_run_gsea.R
#
# Step 7/9 — Gene Set Enrichment Analysis (GSEA) on TF families
#
# All parameters are read from config.sh via load_config.R.
# Edit config.sh to change thresholds — do not edit this script.
#
# Inputs (produced by step 06):
#   gsea_ranked_tfs.rnk      Ranked TF list (tf, rank_score; no header)
#   jaspar_tf_families.gmt   GMT file with TF family gene sets
#
# Outputs:
#   GSEA_TF_families.tsv     fgsea results for all tested TF families
#   GSEA_objects.RData       R objects for downstream steps (08, 09)
# =============================================================================
suppressPackageStartupMessages({
  library(fgsea)
  library(dplyr)
})

.this_script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(trailingOnly=FALSE)[startsWith(commandArgs(trailingOnly=FALSE), "--file=")])))
source(file.path(.this_script_dir, "..", "config.R"))

cat("[STEP 7/9] GSEA analysis\n")
cat(sprintf("  Parameters from config.sh:\n"))
cat(sprintf("    GSEA_MIN_SIZE  = %d\n",   GSEA_MIN_SIZE))
cat(sprintf("    GSEA_MAX_SIZE  = %d\n",   GSEA_MAX_SIZE))
cat(sprintf("    GSEA_N_PERM    = %d\n",   GSEA_N_PERM))
cat(sprintf("    GSEA_SEED      = %d\n",   GSEA_SEED))
cat(sprintf("    BH_FDR_THRESHOLD = %.3f\n", BH_FDR_THRESHOLD))

# =============================================================================
# write_empty_outputs
#
# Written when GSEA is skipped so steps 08 and 09 always find readable files.
# =============================================================================
write_empty_outputs <- function() {
  write.table(
    data.frame(pathway = character(), pval = numeric(), padj = numeric(),
               log2err = numeric(), ES = numeric(), NES = numeric(),
               size = integer(), leadingEdge = character()),
    "GSEA_TF_families.tsv", sep = "\t", quote = FALSE, row.names = FALSE
  )
  save(fg = data.frame(), file = "GSEA_objects.RData")
}

# =============================================================================
# load_gsea_inputs
# =============================================================================
load_gsea_inputs <- function() {
  ranks <- read.table("gsea_ranked_tfs.rnk", stringsAsFactors = FALSE)
  stats <- sort(setNames(ranks$V2, ranks$V1), decreasing = TRUE)
  pathways <- gmtPathways("jaspar_tf_families.gmt")

  overlap <- sum(names(stats) %in% unlist(pathways))
  cat(sprintf("  Ranked TFs: %d | TFs matching GMT pathways: %d\n", length(stats), overlap))
  if (overlap == 0) stop("[ERROR] No overlap between ranked TFs and GMT pathways")

  list(stats = stats, pathways = pathways)
}

# =============================================================================
# run_fgsea
#
# scoreType = "std" is valid because the rank list contains both positive and
# negative values. All parameters from config.sh (GSEA_MIN_SIZE, GSEA_MAX_SIZE,
# GSEA_N_PERM, GSEA_SEED).
# =============================================================================
run_fgsea <- function(stats, pathways) {
  set.seed(GSEA_SEED)

  fg <- suppressWarnings(fgsea(
    pathways    = pathways,
    stats       = stats,
    minSize     = GSEA_MIN_SIZE,
    maxSize     = GSEA_MAX_SIZE,
    scoreType   = "std",
    nPermSimple = GSEA_N_PERM
  )) %>% arrange(padj)

  fg$leadingEdge <- sapply(fg$leadingEdge, paste, collapse = ",")

  cat(sprintf("  %d TF families tested | %d significant (padj < %.3f)\n",
              nrow(fg), sum(fg$padj < BH_FDR_THRESHOLD, na.rm = TRUE), BH_FDR_THRESHOLD))
  fg
}

# =============================================================================
# save_results
# =============================================================================
save_results <- function(fg, stats, pathways) {
  write.table(fg, "GSEA_TF_families.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  save(fg, stats, pathways, file = "GSEA_objects.RData")
}

# =============================================================================
# MAIN
# =============================================================================
if (!file.exists("gsea_ranked_tfs.rnk") ||
    file.info("gsea_ranked_tfs.rnk")$size == 0) {
  cat("  WARNING: gsea_ranked_tfs.rnk is empty — skipping GSEA.\n")
  write_empty_outputs()
  cat("[STEP 7/9] Skipped\n")
  quit(save = "no", status = 0)
}

inputs <- load_gsea_inputs()
fg     <- run_fgsea(inputs$stats, inputs$pathways)
save_results(fg, inputs$stats, inputs$pathways)

cat("[STEP 7/9] Done\n")
