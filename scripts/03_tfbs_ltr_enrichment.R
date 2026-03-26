#!/usr/bin/env Rscript
# =============================================================================
# 03_tfbs_ltr_enrichment.R
#
# Step 3/9 — Fisher exact tests for TFBS enrichment in LTR families
#
# For each (LTR family, TF) pair, tests whether the TF binding motif is
# significantly enriched in that family compared to all other LTR regions
# using a Fisher exact test. Multiple testing correction is applied per
# family using the Benjamini-Hochberg method.
#
# All parameters are read from config.sh via load_config.R.
# Edit config.sh to change thresholds — do not edit this script.
#
# Inputs (produced by steps 01 and 02):
#   LTR_5prime.bed            Full LTR coordinates (BED6)
#   solo_LTR.bed              Solo LTR coordinates (BED6)
#   FIMO_LTR/fimo.tsv         FIMO motif hits in full LTR sequences
#   FIMO_SOLO_LTR/fimo.tsv    FIMO motif hits in solo LTR sequences
#
# Outputs:
#   TFBS_LTR_enrichment_results.tsv   Fisher test results for all (family, TF) pairs
#   TFBS_LTR_part1_objects.RData      R objects for downstream steps (05, 08)
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
})

.this_script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(trailingOnly=FALSE)[startsWith(commandArgs(trailingOnly=FALSE), "--file=")])))
source(file.path(.this_script_dir, "..", "config.R"))

cat("[STEP 3/9] Fisher exact tests (LTR enrichment)\n")
cat(sprintf("  Using %d cores\n", NCPUS))
cat(sprintf("  Parameters from config.sh:\n"))
cat(sprintf("    FISHER_MIN_FAMILY_SIZE    = %d\n",  FISHER_MIN_FAMILY_SIZE))
cat(sprintf("    FISHER_MIN_HITS_IN_FAMILY = %d\n",  FISHER_MIN_HITS_IN_FAMILY))
cat(sprintf("    FISHER_MIN_OUTSIDE_TOTAL  = %d\n",  FISHER_MIN_OUTSIDE_TOTAL))
cat(sprintf("    BH_FDR_THRESHOLD          = %.3f\n", BH_FDR_THRESHOLD))

# =============================================================================
# load_ltr_bed
#
# Loads a BED6 file of LTR regions and extracts family names and ltr_ids.
# ltr_id format: FamilyName::chr:start-end(strand)
# =============================================================================
load_ltr_bed <- function(path, source_label) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  df <- read.table(path, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = c("character", "integer", "integer",
                                  "character", "character", "character"))
  colnames(df) <- c("chr", "start", "end", "name", "dot", "strand")
  df$family <- sub("::.*", "", df$name)
  df$source <- source_label
  df
}

# =============================================================================
# build_ltr_table
# =============================================================================
build_ltr_table <- function() {
  ltr_full <- load_ltr_bed("LTR_5prime.bed", "full")
  ltr_solo <- load_ltr_bed("solo_LTR.bed",   "solo")
  ltr      <- bind_rows(ltr_full, ltr_solo)

  if (is.null(ltr) || nrow(ltr) == 0) stop("[ERROR] No LTR regions found")

  ltr <- ltr %>%
    mutate(
      ltr_id    = paste0(name, "(", strand, ")"),
      length_bp = end - start + 1,
      length_kb = length_bp / 1000,
      genome    = sub("_.*", "", chr)
    )

  cat(sprintf("  Genomes detected : %s\n", paste(sort(unique(ltr$genome)), collapse = ", ")))
  cat(sprintf("  LTR families     : %d unique\n", length(unique(ltr$family))))
  cat(sprintf("  Full LTRs        : %d\n", sum(ltr$source == "full", na.rm = TRUE)))
  cat(sprintf("  Solo LTRs        : %d\n", sum(ltr$source == "solo", na.rm = TRUE)))
  cat(sprintf("  Total LTR regions: %d\n", nrow(ltr)))
  ltr
}

# =============================================================================
# load_fimo_results
# =============================================================================
load_fimo_results <- function() {
  load_fimo <- function(path, label) {
    if (!file.exists(path)) return(NULL)
    df <- tryCatch(read.delim(path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$fimo_source <- label
    df
  }
  fimo_full <- load_fimo("FIMO_LTR/fimo.tsv",      "full")
  fimo_solo <- load_fimo("FIMO_SOLO_LTR/fimo.tsv", "solo")
  fimo_raw  <- bind_rows(fimo_full, fimo_solo)
  if (is.null(fimo_raw) || nrow(fimo_raw) == 0) stop("[ERROR] FIMO found no motifs!")
  cat(sprintf("  FIMO hits (full) : %d\n", if (!is.null(fimo_full)) nrow(fimo_full) else 0))
  cat(sprintf("  FIMO hits (solo) : %d\n", if (!is.null(fimo_solo)) nrow(fimo_solo) else 0))
  fimo_raw
}

# =============================================================================
# parse_fimo_to_ltr_pairs
#
# Parses raw FIMO hits into unique (ltr_id, tf) pairs.
# Strips trailing FIMO position suffix after ")" as a safeguard.
# =============================================================================
parse_fimo_to_ltr_pairs <- function(fimo_raw, ltr) {
  fimo_parsed <- fimo_raw %>%
    mutate(ltr_id = sub("\\)[0-9]+$", ")", sequence_name), tf = motif_alt_id) %>%
    filter(nchar(ltr_id) > 0) %>%
    distinct(ltr_id, tf)

  cat(sprintf("  Parsed FIMO hits : %d unique (ltr_id, tf) pairs\n", nrow(fimo_parsed)))

  n_matched  <- sum(unique(fimo_parsed$ltr_id) %in% unique(ltr$ltr_id))
  n_fimo_ids <- length(unique(fimo_parsed$ltr_id))
  cat(sprintf("  ltr_id join check: %d / %d FIMO ltr_ids matched LTR table\n",
              n_matched, n_fimo_ids))

  if (n_matched == 0) {
    cat("\n  [ERROR] Zero ltr_id matches — format comparison:\n")
    cat(paste0("    FIMO : ", head(unique(fimo_parsed$ltr_id), 3), "\n"))
    cat(paste0("    BED  : ", head(unique(ltr$ltr_id), 3), "\n"))
    stop("ltr_id mismatch — check BED col4 format and bedtools getfasta -s output")
  }
  fimo_parsed
}

# =============================================================================
# build_fisher_input
#
# Constructs 2x2 contingency table entries for each (family, TF) pair.
# Pre-filtering thresholds come from config.sh via load_config.R.
#
# Contingency table:
#                    Has TF motif    No TF motif
#   In family X:          a               b
#   Outside family X:     c               d
# =============================================================================
build_fisher_input <- function(fimo_parsed, ltr, N_total) {
  ltr_family <- ltr %>% select(ltr_id, family)
  n_dup <- sum(duplicated(ltr_family$ltr_id))
  if (n_dup > 0) {
    warning(sprintf("[WARN] %d duplicate ltr_ids — dropping.", n_dup))
    ltr_family <- ltr_family %>% distinct(ltr_id, .keep_all = TRUE)
  }

  family_size    <- ltr %>% count(family, name = "n_family")
  all_families   <- unique(ltr$family)
  all_tfs        <- sort(unique(fimo_parsed$tf))

  cat(sprintf("  Families: %d | TFs: %d | Total LTRs: %d\n",
              length(all_families), length(all_tfs), N_total))

  tf_total <- fimo_parsed %>%
    semi_join(ltr_family, by = "ltr_id") %>%
    count(tf, name = "n_tf_total")

  family_tf_hits <- fimo_parsed %>%
    inner_join(ltr_family, by = "ltr_id", relationship = "many-to-one") %>%
    count(family, tf, name = "a")

  cat(sprintf("  Non-zero (family, tf) pairs: %d / %d possible\n",
              nrow(family_tf_hits), length(all_families) * length(all_tfs)))

  fisher_input <- family_tf_hits %>%
    inner_join(family_size, by = "family") %>%
    filter(n_family >= FISHER_MIN_FAMILY_SIZE) %>%
    inner_join(tf_total, by = "tf") %>%
    mutate(b = n_family - a, c = n_tf_total - a, d = N_total - n_family - c) %>%
    filter(
      a >= FISHER_MIN_HITS_IN_FAMILY,
      a + b >= FISHER_MIN_FAMILY_SIZE,
      c + d >= FISHER_MIN_OUTSIDE_TOTAL,
      b + d > 0, a + c > 0
    )

  cat(sprintf("  Fisher tests to run: %d\n", nrow(fisher_input)))
  fisher_input
}

# =============================================================================
# run_fisher_sparse
#
# Runs a single Fisher exact test. Minimum count checks use thresholds
# from config.sh (available as globals after load_config.R).
# Returns NULL if minimum requirements are not met.
# =============================================================================
run_fisher_sparse <- function(a, b, c, d) {
  if (a + b < FISHER_MIN_FAMILY_SIZE)  return(NULL)
  if (c + d < FISHER_MIN_OUTSIDE_TOTAL) return(NULL)
  if (a < FISHER_MIN_HITS_IN_FAMILY)    return(NULL)
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
  list(odds_ratio = unname(ft$estimate),
       ci_low     = ft$conf.int[1],
       ci_high    = ft$conf.int[2],
       p_value    = ft$p.value)
}

# =============================================================================
# run_fisher_tests
#
# Runs Fisher exact tests for all (family, TF) pairs in parallel.
# BH correction applied per family at BH_FDR_THRESHOLD from config.sh.
# =============================================================================
run_fisher_tests <- function(fisher_input, all_families) {
  cat("  Running Fisher tests...\n")

  results_list <- mclapply(all_families, function(fam) {
    fam_data <- fisher_input %>% filter(family == fam)
    if (nrow(fam_data) == 0) return(NULL)
    res <- lapply(seq_len(nrow(fam_data)), function(i) {
      row <- fam_data[i, ]
      ft  <- run_fisher_sparse(row$a, row$b, row$c, row$d)
      if (is.null(ft)) return(NULL)
      data.frame(family = fam, tf = row$tf,
                 ltr_in_family_with_tf         = row$a,
                 ltr_in_family_without_tf      = row$b,
                 ltr_outside_family_with_tf    = row$c,
                 ltr_outside_family_without_tf = row$d,
                 odds_ratio = ft$odds_ratio, ci_low = ft$ci_low,
                 ci_high = ft$ci_high, p_value = ft$p_value,
                 stringsAsFactors = FALSE)
    })
    bind_rows(Filter(Negate(is.null), res))
  }, mc.cores = NCPUS)

  fisher_df <- bind_rows(Filter(Negate(is.null), results_list)) %>%
    group_by(family) %>%
    mutate(q_value = p.adjust(p_value, method = "BH")) %>%
    ungroup()

  cat(sprintf("  %d Fisher tests completed\n", nrow(fisher_df)))
  fisher_df
}

# =============================================================================
# save_results
#
# Bundles run_fisher_sparse and all config thresholds into the RData so
# downstream scripts (05, 08) use identical parameters without re-reading
# config.sh themselves.
# =============================================================================
save_results <- function(fisher_df, ltr, N_total, all_tfs) {
  write.table(fisher_df, "TFBS_LTR_enrichment_results.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)
  save(ltr, fisher_df, N_total, all_tfs,
       run_fisher_sparse,
       FISHER_MIN_FAMILY_SIZE, FISHER_MIN_HITS_IN_FAMILY, 
       FISHER_MIN_OUTSIDE_TOTAL, BH_FDR_THRESHOLD,
       file = "TFBS_LTR_part1_objects.RData")
}

# =============================================================================
# MAIN
# =============================================================================
ltr          <- build_ltr_table()
fimo_raw     <- load_fimo_results()
fimo_parsed  <- parse_fimo_to_ltr_pairs(fimo_raw, ltr)
N_total      <- nrow(ltr)
all_families <- unique(ltr$family)
all_tfs      <- sort(unique(fimo_parsed$tf))
fisher_input <- build_fisher_input(fimo_parsed, ltr, N_total)
fisher_df    <- run_fisher_tests(fisher_input, all_families)
save_results(fisher_df, ltr, N_total, all_tfs)

cat("[STEP 3/9] Done\n")
