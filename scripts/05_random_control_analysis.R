#!/usr/bin/env Rscript
# =============================================================================
# 05_random_control_analysis.R
#
# Step 5/9 — Fisher exact tests on genomic shuffle controls and
#            LTR vs control comparison
#
# All parameters are read from config.sh via load_config.R.
# Edit config.sh to change thresholds — do not edit this script.
#
# Inputs (produced by steps 03 and 04):
#   TFBS_LTR_part1_objects.RData      LTR Fisher results and helper objects
#   all_ltr_combined.bed              Combined LTR BED (family labels)
#   random_genomic.bed                Genomic shuffle control coordinates
#   FIMO_RANDOM_GENOMIC/fimo.tsv      FIMO hits in control sequences
#
# Outputs:
#   ctrl_genomic_enrichment_results.tsv   Fisher results for genomic controls
#   ltr_vs_ctrl_comparison.tsv            LTR vs control comparison table
#   TFBS_random_controls_objects.RData    R objects for downstream steps (08)
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
})

.this_script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(trailingOnly=FALSE)[startsWith(commandArgs(trailingOnly=FALSE), "--file=")])))
source(file.path(.this_script_dir, "..", "config.R"))

cat("[STEP 5/9] Random control Fisher exact tests\n")
cat(sprintf("  Using %d cores\n", NCPUS))
cat(sprintf("  Parameters from config.sh:\n"))
cat(sprintf("    FISHER_MIN_FAMILY_SIZE   = %d\n",   FISHER_MIN_FAMILY_SIZE))
cat(sprintf("    FISHER_MIN_HITS_IN_FAMILY= %d\n",   FISHER_MIN_HITS_IN_FAMILY))
cat(sprintf("    LTR_SPECIFIC_LOG2OR      = %.2f\n",  LTR_SPECIFIC_LOG2OR))
cat(sprintf("    LTR_DEPLETED_LOG2OR      = %.2f\n",  LTR_DEPLETED_LOG2OR))
cat(sprintf("    BH_FDR_THRESHOLD         = %.3f\n",  BH_FDR_THRESHOLD))

# Load RData from step 03 — brings in fisher_df, ltr, N_total, all_tfs,
# run_fisher_sparse, and the Fisher threshold variables
load("TFBS_LTR_part1_objects.RData")
stopifnot(exists("N_total"), exists("all_tfs"), exists("fisher_df"),
          exists("run_fisher_sparse"))

# =============================================================================
# parse_family
# =============================================================================
parse_family <- function(x) sub("::.*", "", x)

# =============================================================================
# load_control_beds
#
# Loads combined LTR BED and genomic shuffle BED.
# Family labels are inherited from combined LTR BED by row position
# (bedtools shuffle preserves input row order).
# =============================================================================
load_control_beds <- function() {
  cat("  Loading genomic shuffle control...\n")

  read_bed6 <- function(path) {
    df <- read.table(path, sep = "\t", stringsAsFactors = FALSE,
                     colClasses = c("character","integer","integer",
                                    "character","character","character"))
    colnames(df) <- c("chr","start","end","name","dot","strand")
    df
  }

  all_ltr_bed <- read_bed6("all_ltr_combined.bed")
  all_ltr_bed$family <- parse_family(all_ltr_bed$name)

  genomic_bed <- read_bed6("random_genomic.bed")

  if (nrow(genomic_bed) != nrow(all_ltr_bed)) {
    cat(sprintf(
      "  WARNING: bedtools shuffle dropped %d / %d regions (%.1f%%) — continuing\n",
      nrow(all_ltr_bed) - nrow(genomic_bed), nrow(all_ltr_bed),
      100 * (nrow(all_ltr_bed) - nrow(genomic_bed)) / nrow(all_ltr_bed)
    ))
  }

  genomic_bed$family    <- all_ltr_bed$family[seq_len(nrow(genomic_bed))]
  genomic_bed$region_id <- genomic_bed$name

  cat(sprintf("  Genomic control regions: %d\n", nrow(genomic_bed)))
  list(all_ltr_bed = all_ltr_bed, genomic_bed = genomic_bed)
}

# =============================================================================
# load_fimo_control
#
# Loads FIMO hits for control sequences. Deduplicates to unique
# (region_id, tf) pairs and checks ID overlap.
# =============================================================================
load_fimo_control <- function(path, genomic_bed) {
  if (!file.exists(path)) {
    cat(sprintf("  WARNING: %s not found\n", path))
    return(data.frame(region_id = character(), tf = character()))
  }
  df <- read.delim(path, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(data.frame(region_id = character(), tf = character()))

  fimo_ctrl <- df %>%
    transmute(region_id = sub("\\)[0-9]+$", ")", sequence_name), tf = motif_alt_id) %>%
    distinct(region_id, tf)

  cat(sprintf("  FIMO genomic hits: %d unique (region, tf) pairs\n", nrow(fimo_ctrl)))

  n_overlap <- sum(unique(fimo_ctrl$region_id) %in% unique(genomic_bed$region_id))
  cat(sprintf("  Genomic ID overlap: %d / %d FIMO IDs matched\n",
              n_overlap, length(unique(fimo_ctrl$region_id))))

  if (n_overlap == 0 && nrow(fimo_ctrl) > 0) {
    cat("  [ERROR] Zero region_id matches — format comparison:\n")
    cat(paste0("    FIMO: ", head(unique(fimo_ctrl$region_id),   3), "\n"))
    cat(paste0("    BED : ", head(unique(genomic_bed$region_id), 3), "\n"))
    stop("region_id mismatch in genomic control")
  }
  fimo_ctrl
}

# =============================================================================
# run_ctrl_fisher
#
# Runs Fisher exact tests on control regions using the same logic as step 03.
# Thresholds come from config.sh via load_config.R (already in global env).
# =============================================================================
run_ctrl_fisher <- function(fimo_hits, regions_df, N_ctrl, label) {
  family_size_ctrl <- regions_df %>% count(family, name = "n_family")

  tf_total_ctrl <- fimo_hits %>%
    inner_join(regions_df %>% select(region_id, family),
               by = "region_id", relationship = "many-to-one") %>%
    count(tf, name = "n_tf_total")

  family_tf_hits_ctrl <- fimo_hits %>%
    inner_join(regions_df %>% select(region_id, family),
               by = "region_id", relationship = "many-to-one") %>%
    count(family, tf, name = "a")

  fisher_input_ctrl <- family_tf_hits_ctrl %>%
    inner_join(family_size_ctrl, by = "family") %>%
    filter(n_family >= FISHER_MIN_FAMILY_SIZE) %>%
    inner_join(tf_total_ctrl, by = "tf") %>%
    mutate(b = n_family - a, c = n_tf_total - a, d = N_ctrl - n_family - c) %>%
    filter(
      a >= FISHER_MIN_HITS_IN_FAMILY,
      a + b >= FISHER_MIN_FAMILY_SIZE,
      c + d >= FISHER_MIN_OUTSIDE_TOTAL,
      b + d > 0, a + c > 0
    )

  cat(sprintf("  [%s] Fisher tests to run: %d\n", label, nrow(fisher_input_ctrl)))

  if (nrow(fisher_input_ctrl) == 0) {
    cat(sprintf("  [%s] No testable pairs — returning empty\n", label))
    return(data.frame(
      family = character(), tf = character(),
      ctrl_in_family_with_tf = integer(), ctrl_in_family_without_tf = integer(),
      ctrl_outside_family_with_tf = integer(), ctrl_outside_family_without_tf = integer(),
      odds_ratio = numeric(), ci_low = numeric(), ci_high = numeric(),
      p_value = numeric(), q_value = numeric(), ctrl_label = character(),
      stringsAsFactors = FALSE
    ))
  }

  results_list <- mclapply(unique(fisher_input_ctrl$family), function(fam) {
    fam_data <- fisher_input_ctrl %>% filter(family == fam)
    res <- lapply(seq_len(nrow(fam_data)), function(i) {
      row <- fam_data[i, ]
      ft  <- run_fisher_sparse(row$a, row$b, row$c, row$d)
      if (is.null(ft)) return(NULL)
      data.frame(family = fam, tf = row$tf,
                 ctrl_in_family_with_tf         = row$a,
                 ctrl_in_family_without_tf      = row$b,
                 ctrl_outside_family_with_tf    = row$c,
                 ctrl_outside_family_without_tf = row$d,
                 odds_ratio = ft$odds_ratio, ci_low = ft$ci_low,
                 ci_high = ft$ci_high, p_value = ft$p_value,
                 stringsAsFactors = FALSE)
    })
    bind_rows(Filter(Negate(is.null), res))
  }, mc.cores = NCPUS)

  bind_rows(Filter(Negate(is.null), results_list)) %>%
    group_by(family) %>%
    mutate(q_value = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(ctrl_label = label)
}

# =============================================================================
# build_comparison_table
#
# Merges LTR and control Fisher results and classifies each pair.
# Category thresholds (LTR_SPECIFIC_LOG2OR, LTR_DEPLETED_LOG2OR,
# BH_FDR_THRESHOLD) come from config.sh via load_config.R.
#
# Categories:
#   "LTR-specific"               q_ltr < BH_FDR, OR_ltr > LTR_SPECIFIC, OR_ctrl < LTR_SPECIFIC
#   "Shared (possible artefact)" q_ltr < BH_FDR, OR_ltr > LTR_SPECIFIC, OR_ctrl >= LTR_SPECIFIC
#   "Depleted in LTR"            q_ltr < BH_FDR, OR_ltr < LTR_DEPLETED
#   "Not significant"            all other pairs
# =============================================================================
build_comparison_table <- function(fisher_df, genomic_fisher_df) {
  comparison_out <- merge(
    fisher_df[, c("family","tf","odds_ratio","q_value")],
    genomic_fisher_df[, c("family","tf","odds_ratio","q_value")],
    by = c("family","tf"), suffixes = c("_ltr","_ctrl")
  )

  comparison_out$log2OR_ltr  <- log2(pmax(comparison_out$odds_ratio_ltr,  1e-6))
  comparison_out$log2OR_ctrl <- log2(pmax(comparison_out$odds_ratio_ctrl, 1e-6))
  comparison_out <- comparison_out[
    is.finite(comparison_out$log2OR_ltr) & is.finite(comparison_out$log2OR_ctrl), ]

  comparison_out$category <- ifelse(
    comparison_out$q_value_ltr < BH_FDR_THRESHOLD &
      comparison_out$log2OR_ltr  >  LTR_SPECIFIC_LOG2OR &
      comparison_out$log2OR_ctrl <  LTR_SPECIFIC_LOG2OR,
    "LTR-specific",
    ifelse(
      comparison_out$q_value_ltr < BH_FDR_THRESHOLD &
        comparison_out$log2OR_ltr  >  LTR_SPECIFIC_LOG2OR &
        comparison_out$log2OR_ctrl >= LTR_SPECIFIC_LOG2OR,
      "Shared (possible artefact)",
      ifelse(
        comparison_out$q_value_ltr < BH_FDR_THRESHOLD &
          comparison_out$log2OR_ltr  <  LTR_DEPLETED_LOG2OR,
        "Depleted in LTR",
        "Not significant"
      )
    )
  )

  comparison_out[order(comparison_out$log2OR_ltr - comparison_out$log2OR_ctrl,
                       decreasing = TRUE), ]
}

# =============================================================================
# save_results
# =============================================================================
save_results <- function(genomic_fisher_df, comparison_out) {
  write.table(genomic_fisher_df, "ctrl_genomic_enrichment_results.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(comparison_out, "ltr_vs_ctrl_comparison.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)
  save(genomic_fisher_df, comparison_out,
       file = "TFBS_random_controls_objects.RData")
}

# =============================================================================
# MAIN
# =============================================================================
beds         <- load_control_beds()
fimo_genomic <- load_fimo_control("FIMO_RANDOM_GENOMIC/fimo.tsv", beds$genomic_bed)

genomic_fisher_df <- run_ctrl_fisher(
  fimo_genomic, beds$genomic_bed, nrow(beds$genomic_bed), "RANDOM_GENOMIC"
)

comparison_out <- build_comparison_table(fisher_df, genomic_fisher_df)
save_results(genomic_fisher_df, comparison_out)

cat(sprintf("  Genomic control: %d (family,tf) pairs tested\n", nrow(genomic_fisher_df)))
cat("[STEP 5/9] Done\n")
