#!/usr/bin/env Rscript
# =============================================================================
# 09_report.R  (GLORIA — summary report)
#
# Step 9/9 — Generate a plain-text summary report of all analysis results
#
# All parameters are read from config.sh via load_config.R.
# Edit config.sh to change thresholds — do not edit this script.
#
# Usage:
#   Rscript 09_report.R
#   Rscript 09_report.R results.tsv
#   Rscript 09_report.R results.tsv GSEA_TF_families.tsv
#   Rscript 09_report.R results.tsv GSEA_TF_families.tsv report.txt
# =============================================================================

# =============================================================================
# setup_output_sink  (called before source() so it captures all output)
# =============================================================================
setup_output_sink <- function(output_path) {
  con <- file(output_path, open = "wt")
  sink(con, type = "output")
  sink(con, type = "message")
  con
}

suppressPackageStartupMessages(library(dplyr))

.this_script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(trailingOnly=FALSE)[startsWith(commandArgs(trailingOnly=FALSE), "--file=")])))
source(file.path(.this_script_dir, "..", "config.R"))

# Map config.sh names to local aliases used throughout report functions
Q_THRESH                    <- REPORT_Q_THRESH
OR_THRESH                   <- REPORT_OR_THRESH

# =============================================================================
# FORMATTING HELPERS
# =============================================================================
section_sep <- function(label = "") {
  cat("\n", strrep("-", 60), "\n", sep = "")
  if (nchar(label) > 0) cat(" ", label, "\n", sep = "")
  cat(strrep("-", 60), "\n", sep = "")
}

fmt_pct <- function(x, n) sprintf("%d (%.1f%%)", x, 100 * x / n)

print_table <- function(df) {
  if (nrow(df) == 0) { cat("  (no records)\n"); return(invisible()) }
  col_widths <- sapply(seq_along(df), function(i)
    max(nchar(colnames(df)[i]), max(nchar(as.character(df[[i]])), na.rm = TRUE)))
  header <- paste(mapply(formatC, colnames(df), width = col_widths, flag = "-"), collapse = "  ")
  cat("  ", header, "\n", sep = "")
  cat("  ", strrep("-", nchar(header)), "\n", sep = "")
  for (i in seq_len(nrow(df))) {
    row <- paste(mapply(formatC, as.character(unlist(df[i,])),
                        width = col_widths, flag = "-"), collapse = "  ")
    cat("  ", row, "\n", sep = "")
  }
}

fmt_qvalue <- function(x)
  ifelse(x == 0, "< 2.22e-308", formatC(x, format = "e", digits = 2))

# =============================================================================
# DATA LOADING
# =============================================================================
load_fisher_results <- function(path) {
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
  df <- read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  required <- c("family","tf","odds_ratio","q_value","ltr_in_family_with_tf",
                "ltr_in_family_without_tf","ltr_outside_family_with_tf",
                "ltr_outside_family_without_tf")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop(sprintf("Missing columns: %s", paste(missing, collapse = ", ")))
  df %>% mutate(log2OR = ifelse(is.finite(odds_ratio) & odds_ratio > 0,
                                log2(odds_ratio), NA_real_))
}

load_ltr_bed <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  tryCatch({
    bed <- read.table(path, sep = "\t", stringsAsFactors = FALSE,
                      colClasses = c("character","integer","integer",
                                     "character","character","character"))
    colnames(bed) <- c("chr","start","end","name","dot","strand")
    bed$family    <- sub("::.*", "", bed$name)
    bed$length_bp <- bed$end - bed$start
    bed
  }, error = function(e) NULL)
}

detect_genomes_from_bed <- function() {
  genomes <- character(0)
  for (bed_path in c("LTR_5prime.bed","solo_LTR.bed")) {
    if (!file.exists(bed_path) || file.info(bed_path)$size == 0) next
    bed_tmp <- read.table(bed_path, sep = "\t", stringsAsFactors = FALSE,
                          colClasses = c("character", rep("NULL", 5)))
    genomes <- c(genomes, sub("_.*", "", bed_tmp[[1]]))
  }
  sort(unique(genomes[genomes != ""]))
}

# =============================================================================
# REPORT SECTIONS
# =============================================================================
report_header <- function() {
  cat(strrep("=", 60), "\n")
  cat(" GLORIA — SUMMARY REPORT\n")
  cat(strrep("=", 60), "\n")
  genomes_detected <- detect_genomes_from_bed()
  if (length(genomes_detected) > 0)
    cat(sprintf("  Genomes             : %s\n", paste(genomes_detected, collapse = ", ")))
  cat(sprintf("  q-value threshold   : %.2f  (REPORT_Q_THRESH in config.sh)\n",  Q_THRESH))
  cat(sprintf("  OR threshold        : %.1f  (REPORT_OR_THRESH in config.sh)\n", OR_THRESH))
}

report_ltr_counts <- function() {
  section_sep("LTR REGION COUNTS")
  bed_full <- load_ltr_bed("LTR_5prime.bed")
  bed_solo <- load_ltr_bed("solo_LTR.bed")
  has_full <- !is.null(bed_full) && nrow(bed_full) > 0
  has_solo <- !is.null(bed_solo) && nrow(bed_solo) > 0

  if (!has_full && !has_solo) {
    cat("  BED files not found — run from the output run directory.\n")
    return(invisible())
  }

  n_full <- if (has_full) nrow(bed_full) else 0
  n_solo <- if (has_solo) nrow(bed_solo) else 0
  cat(sprintf("  Full LTR (LTR_5prime.bed) : %d\n", n_full))
  cat(sprintf("  Solo LTR (solo_LTR.bed)   : %d\n", n_solo))
  cat(sprintf("  Total                     : %d\n", n_full + n_solo))

  if (has_full)
    cat(sprintf("\n  Full LTR length [bp] — median: %d  mean: %d  min: %d  max: %d\n",
                as.integer(median(bed_full$length_bp)), as.integer(mean(bed_full$length_bp)),
                min(bed_full$length_bp), max(bed_full$length_bp)))
  if (has_solo)
    cat(sprintf("  Solo LTR length [bp] — median: %d  mean: %d  min: %d  max: %d\n",
                as.integer(median(bed_solo$length_bp)), as.integer(mean(bed_solo$length_bp)),
                min(bed_solo$length_bp), max(bed_solo$length_bp)))

  cat("\n  LTR region counts per family:\n\n")
  all_bed <- bind_rows(
    if (has_full) bed_full %>% mutate(source = "full") else NULL,
    if (has_solo) bed_solo %>% mutate(source = "solo") else NULL
  )
  family_counts <- all_bed %>%
    group_by(family) %>%
    summarise(n_full = sum(source == "full"), n_solo = sum(source == "solo"),
              n_total = n(), med_len = as.integer(median(length_bp)),
              min_len = min(length_bp), max_len = max(length_bp), .groups = "drop") %>%
    arrange(desc(n_total))
  print_table(family_counts)

  # Dominant family warning (threshold from config.sh)
  dominant_fam <- family_counts %>%
    mutate(pct = n_total / sum(n_total)) %>% filter(pct > REPORT_DOMINANT_FAMILY_FRAC)
  if (nrow(dominant_fam) > 0) {
    cat(sprintf("\n  !!! WARNING: %s accounts for %.1f%% of all LTR regions.\n",
                dominant_fam$family[1], dominant_fam$pct[1] * 100))
    cat("  Fisher tests for this family may be unreliable.\n")
    cat("  Use genomic control results (steps 4-5) instead.\n")
  }

  # Per-genome breakdown
  genomes_bed <- all_bed %>% mutate(genome = sub("_.*", "", chr)) %>% filter(genome != chr)
  if (n_distinct(genomes_bed$genome) > 1) {
    cat("\n  LTR counts per genome:\n\n")
    genome_counts <- genomes_bed %>%
      group_by(genome) %>%
      summarise(n_full = sum(source == "full"), n_solo = sum(source == "solo"),
                n_total = n(), .groups = "drop") %>% arrange(desc(n_total))
    print_table(genome_counts)
  }

  # Small family warning (threshold from config.sh)
  small_families <- family_counts %>% filter(n_total < REPORT_SMALL_FAMILY_WARN)
  if (nrow(small_families) > 0) {
    cat(sprintf("\n  WARNING: %d familie(s) have fewer than %d LTR regions:\n",
                nrow(small_families), REPORT_SMALL_FAMILY_WARN))
    cat(sprintf("    %s\n",
                paste(small_families$family, sprintf("(n=%d)", small_families$n_total),
                      collapse = ", ")))
  }
}

report_basic_overview <- function(df) {
  section_sep("BASIC OVERVIEW")
  n_total  <- nrow(df)
  n_fam    <- n_distinct(df$family)
  n_tfs    <- n_distinct(df$tf)
  n_enrich <- sum(df$q_value < Q_THRESH & df$odds_ratio > 1,         na.rm = TRUE)
  n_deplet <- sum(df$q_value < Q_THRESH & df$odds_ratio < 1,         na.rm = TRUE)
  n_strong <- sum(df$q_value < Q_THRESH & df$odds_ratio > OR_THRESH, na.rm = TRUE)

  cat(sprintf("  Tested pairs (family x tf)               : %d\n", n_total))
  cat(sprintf("  Tested LTR families                      : %d\n", n_fam))
  cat(sprintf("  Unique TFs                               : %d\n", n_tfs))
  cat(sprintf("\n"))
  cat(sprintf("  Enriched          (q < %.2f, OR > 1)     : %s\n", Q_THRESH, fmt_pct(n_enrich, n_total)))
  cat(sprintf("  Depleted          (q < %.2f, OR < 1)     : %s\n", Q_THRESH, fmt_pct(n_deplet, n_total)))
  cat(sprintf("  Strongly enriched (q < %.2f, OR > %.0f)     : %s\n",
              Q_THRESH, OR_THRESH, fmt_pct(n_strong, n_total)))

  list(n_total = n_total, n_families = n_fam, n_tfs = n_tfs,
       n_enrich = n_enrich, n_deplet = n_deplet, n_strong = n_strong)
}

report_top20 <- function(df) {
  section_sep("TOP 20 STRONGEST ASSOCIATIONS  (q < threshold, ranked by OR)")
  top20 <- df %>%
    filter(q_value < Q_THRESH, !is.na(log2OR)) %>%
    arrange(desc(log2OR)) %>% slice_head(n = 20) %>%
    transmute(family, tf, log2OR = round(log2OR, 2), q_value = fmt_qvalue(q_value),
              n_hits = ltr_in_family_with_tf,
              n_family = ltr_in_family_with_tf + ltr_in_family_without_tf)
  print_table(top20)
}

report_per_family <- function(df) {
  section_sep("PER-FAMILY OVERVIEW")
  family_summary <- df %>%
    group_by(family) %>%
    summarise(n_tested   = n(),
              n_sig      = sum(q_value < Q_THRESH,                            na.rm = TRUE),
              n_strong   = sum(q_value < Q_THRESH & odds_ratio > OR_THRESH,  na.rm = TRUE),
              n_deplete  = sum(q_value < Q_THRESH & odds_ratio < 1,          na.rm = TRUE),
              med_log2OR = round(median(log2OR, na.rm = TRUE), 2),
              max_log2OR = round(max(log2OR,    na.rm = TRUE), 2),
              top_tf     = tf[which.min(q_value)], .groups = "drop") %>%
    arrange(desc(n_sig))
  print_table(family_summary %>% select(-top_tf))
  cat("\n  Top TF per family (lowest q-value):\n")
  for (i in seq_len(nrow(family_summary)))
    cat(sprintf("    %-50s %s\n", family_summary$family[i], family_summary$top_tf[i]))
  family_summary
}

report_consistent_tfs <- function(df) {
  section_sep("TOP TFs — CONSISTENCY ACROSS FAMILIES")
  tf_summary <- df %>%
    filter(!is.na(log2OR)) %>%
    group_by(tf) %>%
    summarise(n_tested        = n(),
              n_sig           = sum(q_value < Q_THRESH, na.rm = TRUE),
              mean_log2OR     = round(mean(log2OR, na.rm = TRUE), 2),
              sd_log2OR       = round(sd(log2OR,   na.rm = TRUE), 2),
              max_log2OR      = round(max(log2OR,  na.rm = TRUE), 2),
              consistency_pct = round(100 * sum(q_value < Q_THRESH, na.rm = TRUE) / n(), 1),
              .groups         = "drop") %>%
    filter(n_tested >= REPORT_MIN_TESTED_FAMILIES) %>%
    arrange(desc(n_sig), desc(mean_log2OR)) %>% slice_head(n = 20)

  cat(sprintf("  TFs enriched in the most LTR families (min. %d tested families):\n\n",
              REPORT_MIN_TESTED_FAMILIES))
  print_table(tf_summary)
  tf_summary
}

report_no_signal_families <- function(family_summary) {
  section_sep("LTR FAMILIES WITH NO SIGNIFICANT ENRICHMENT")
  no_signal <- family_summary %>% filter(n_sig == 0) %>%
    select(family, n_tested, med_log2OR, max_log2OR)
  if (nrow(no_signal) > 0) {
    cat(sprintf("  %d familie(s) have no significant TF associations:\n\n", nrow(no_signal)))
    print_table(no_signal)
    cat("\n  Possible reasons: low LTR count, low sequence diversity,\n")
    cat("  biologically silent family, insufficient motif coverage.\n")
  } else {
    cat("  All LTR families have at least one significant association.\n")
  }
}

report_gsea <- function(gsea_file) {
  section_sep("GSEA — TF FAMILY ENRICHMENT")
  if (!file.exists(gsea_file)) {
    cat(sprintf("  GSEA file not found (%s) — section skipped.\n", gsea_file))
    return(invisible())
  }
  gsea <- read.table(gsea_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  if (nrow(gsea) == 0 || !"padj" %in% colnames(gsea)) {
    cat("  GSEA file is empty — no significant TF families.\n")
    return(invisible())
  }
  gsea_sig <- gsea %>% filter(padj < Q_THRESH) %>% arrange(desc(NES))
  cat(sprintf("  TF families tested            : %d\n",  nrow(gsea)))
  cat(sprintf("  Significant (padj < %.2f)     : %d\n", Q_THRESH, nrow(gsea_sig)))
  cat(sprintf("  of which enriched  (NES > 0)  : %d\n", sum(gsea_sig$NES > 0)))
  cat(sprintf("  of which depleted  (NES < 0)  : %d\n", sum(gsea_sig$NES < 0)))
  if (nrow(gsea_sig) > 0) {
    cat("\n  Significant TF families:\n\n")
    print_table(gsea_sig %>%
      transmute(TF_family = pathway, NES = round(NES, 3),
                padj = formatC(padj, format = "e", digits = 2), n_TF = size))
    cat("\n  Leading edge TFs (top contributors to enrichment):\n\n")
    for (i in seq_len(min(5, nrow(gsea_sig))))
      cat(sprintf("  %s:\n    %s\n\n", gsea_sig$pathway[i], gsea_sig$leadingEdge[i]))
  }
}

report_ltr_vs_ctrl <- function() {
  comp_file <- "ltr_vs_ctrl_comparison.tsv"
  if (!file.exists(comp_file)) return(invisible())
  comp <- read.table(comp_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cat("\n------------------------------------------------------------\n")
  cat(" LTR vs GENOMIC CONTROL — CATEGORY BREAKDOWN\n")
  cat("------------------------------------------------------------\n")
  counts <- table(comp$category)
  n_comp <- nrow(comp)
  categories <- list(
    list(key = "LTR-specific",
         lbl = "LTR-specific                (q < thresh, log2OR_ltr > 1, log2OR_ctrl < 1) "),
    list(key = "Shared (possible artefact)",
         lbl = "Shared (possible artefact)  (q < thresh, log2OR_ltr > 1, log2OR_ctrl >= 1)"),
    list(key = "Depleted in LTR",
         lbl = "Depleted in LTR             (q < thresh, log2OR_ltr < -1)                 "),
    list(key = "Not significant",
         lbl = "Not significant             (q > thresh)                                  ")
  )
  for (cat_item in categories) {
    n <- if (cat_item$key %in% names(counts)) counts[[cat_item$key]] else 0
    cat(sprintf("  %s : %d (%.1f%%)\n", cat_item$lbl, n, 100 * n / n_comp))
  }

  cat("\n  Top 10 LTR-specific associations:\n\n")
  top_sp <- head(comp[comp$category == "LTR-specific" & is.finite(comp$log2OR_ltr),
                       ][order(-comp[comp$category == "LTR-specific" &
                                     is.finite(comp$log2OR_ltr), "log2OR_ltr"]), ], 10)
  if (nrow(top_sp) > 0) {
    fw <- max(nchar(top_sp$family), nchar("family")) + 2
    tw <- max(nchar(top_sp$tf),     nchar("tf"))     + 2
    cat(sprintf("  %-*s  %-*s  %10s  %10s  %11s\n", fw,"family", tw,"tf",
                "ltr_log2OR","ctrl_log2OR","q_ltr"))
    cat(sprintf("  %s\n", strrep("-", fw + tw + 37)))
    for (i in seq_len(nrow(top_sp)))
      cat(sprintf("  %-*s  %-*s  %10.2f  %10.2f  %11s\n",
                  fw, top_sp$family[i], tw, top_sp$tf[i],
                  top_sp$log2OR_ltr[i], top_sp$log2OR_ctrl[i],
                  fmt_qvalue(top_sp$q_value_ltr[i])))
  }
}

report_summary <- function(df, counts, family_summary, tf_summary, gsea_file) {
  section_sep("SUMMARY")
  top_pair <- df %>% filter(q_value < Q_THRESH, !is.na(log2OR)) %>%
    arrange(desc(log2OR)) %>% slice(1)
  cat(sprintf("  Analysis tested %d (LTR family x TF) pairs across %d LTR families\n",
              counts$n_total, counts$n_families))
  cat(sprintf("  and %d unique transcription factors.\n\n", counts$n_tfs))
  cat(sprintf("  In total %d (%.1f%%) pairs show significant enrichment (q < %.2f),\n",
              counts$n_enrich, 100 * counts$n_enrich / counts$n_total, Q_THRESH))
  cat(sprintf("  of which %d pairs have a strong effect size (OR > %.0f).\n\n",
              counts$n_strong, OR_THRESH))
  cat(sprintf("  Strongest single association      : %s in family %s (log2OR = %.2f)\n",
              top_pair$tf, top_pair$family, top_pair$log2OR))
  cat(sprintf("  Family with most enrichments      : %s (%d significant TFs)\n",
              family_summary$family[1], family_summary$n_sig[1]))
  cat(sprintf("  Most consistent TF                : %s (%d families, mean log2OR = %.2f)\n",
              tf_summary$tf[1], tf_summary$n_sig[1], tf_summary$mean_log2OR[1]))
  if (file.exists(gsea_file)) {
    gsea_check <- tryCatch(read.table(gsea_file, sep="\t", header=TRUE, stringsAsFactors=FALSE),
                           error = function(e) NULL)
    if (!is.null(gsea_check) && "padj" %in% colnames(gsea_check))
      cat(sprintf("  Significant TF families (GSEA)    : %d\n",
                  sum(gsea_check$padj < Q_THRESH, na.rm = TRUE)))
  }
}

# =============================================================================
# MAIN
# =============================================================================
args_pre   <- commandArgs(trailingOnly = TRUE)
report_out <- if (length(args_pre) >= 3) args_pre[3] else "pipeline_report.txt"

sink_con <- setup_output_sink(report_out)
on.exit({ sink(type = "message"); sink(type = "output"); close(sink_con) }, add = TRUE)

args       <- commandArgs(trailingOnly = TRUE)
INPUT_FILE <- if (length(args) > 0) args[1] else "TFBS_LTR_enrichment_results.tsv"
GSEA_FILE  <- if (length(args) > 1) args[2] else "GSEA_TF_families.tsv"

df     <- load_fisher_results(INPUT_FILE)
has_ci <- all(c("ci_low","ci_high") %in% colnames(df))

report_header()
report_ltr_counts()
counts         <- report_basic_overview(df)
report_top20(df)
family_summary <- report_per_family(df)
tf_summary     <- report_consistent_tfs(df)
report_no_signal_families(family_summary)
report_gsea(GSEA_FILE)
report_ltr_vs_ctrl()
report_summary(df, counts, family_summary, tf_summary, GSEA_FILE)
