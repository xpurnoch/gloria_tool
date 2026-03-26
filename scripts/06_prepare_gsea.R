#!/usr/bin/env Rscript
# =============================================================================
# 06_prepare_gsea.R
#
# Step 6/9 — Prepare ranked TF list and GMT gene set file for GSEA
#
# All parameters are read from config.sh via load_config.R.
# Edit config.sh to change thresholds — do not edit this script.
#
# Inputs (produced by step 03):
#   TFBS_LTR_enrichment_results.tsv   Fisher test results
#
# Additional inputs:
#   jaspar_tf_families.csv   TF-to-family mapping (semicolon-separated;
#                            required columns: Name, Family)
#
# Outputs:
#   gsea_ranked_tfs.rnk      Two-column TSV: tf, rank_score (no header)
#   jaspar_tf_families.gmt   GMT file with one gene set per TF family
# =============================================================================
suppressPackageStartupMessages(library(dplyr))

.this_script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(trailingOnly=FALSE)[startsWith(commandArgs(trailingOnly=FALSE), "--file=")])))
source(file.path(.this_script_dir, "..", "config.R"))

cat("[STEP 6/9] GSEA preparation\n")
cat(sprintf("  BH_FDR_THRESHOLD = %.3f\n", BH_FDR_THRESHOLD))

# =============================================================================
# load_inputs
# =============================================================================
load_inputs <- function() {
  fisher <- read.table("TFBS_LTR_enrichment_results.tsv",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  tf_map <- read.table(JASPAR_FILE, sep = ";", header = TRUE, stringsAsFactors = FALSE)
  stopifnot("Family" %in% colnames(tf_map), "Name" %in% colnames(tf_map))

  list(fisher = fisher, tf_map = tf_map)
}

# =============================================================================
# compute_rank_scores
#
# Rank score = mean log2(OR) across all tested LTR families (finite OR only).
# All testable TFs are included regardless of significance so fgsea receives
# a complete ranked list and scoreType = "std" is valid (both directions).
# =============================================================================
compute_rank_scores <- function(fisher, tf_map) {
  rank_df <- fisher %>%
    filter(is.finite(odds_ratio), odds_ratio > 0) %>%
    mutate(log2OR = log2(odds_ratio)) %>%
    group_by(tf) %>%
    summarise(rank_score = mean(log2OR, na.rm = TRUE),
              n_families = n(),
              n_sig      = sum(q_value < BH_FDR_THRESHOLD, na.rm = TRUE),
              .groups    = "drop") %>%
    arrange(desc(rank_score))

  tf_in_map <- sum(rank_df$tf %in% tf_map$Name)
  cat(sprintf("  TFs in rank file: %d | matching JASPAR map: %d\n", nrow(rank_df), tf_in_map))
  if (tf_in_map == 0)
    warning("[WARN] No overlap between ranked TFs and JASPAR tf_map — GSEA will find nothing")

  rank_df
}

# =============================================================================
# write_rnk_file
# =============================================================================
write_rnk_file <- function(rank_df, output_path) {
  write.table(rank_df[, c("tf","rank_score")], output_path,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# =============================================================================
# build_gmt_file
#
# Format: FamilyName <TAB> NA <TAB> TF1 <TAB> TF2 ...
# =============================================================================
build_gmt_file <- function(tf_map, output_path) {
  gene_sets <- tf_map %>%
    group_by(Family) %>%
    summarise(genes = list(unique(Name)), .groups = "drop")

  gmt <- sapply(seq_len(nrow(gene_sets)), function(i)
    paste(gene_sets$Family[i], "NA", paste(gene_sets$genes[[i]], collapse = "\t"), sep = "\t"))

  writeLines(gmt, output_path)
  nrow(gene_sets)
}

# =============================================================================
# MAIN
# =============================================================================
inputs     <- load_inputs()
rank_df    <- compute_rank_scores(inputs$fisher, inputs$tf_map)

write_rnk_file(rank_df, "gsea_ranked_tfs.rnk")
n_families <- build_gmt_file(inputs$tf_map, "jaspar_tf_families.gmt")

cat(sprintf("  %d TFs ranked | %d TF families in GMT\n", nrow(rank_df), n_families))
cat("[STEP 6/9] Done\n")
