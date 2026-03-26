# =============================================================================
# config.R — GLORIA pipeline parameters for R scripts
#
# This file mirrors the numerical parameters from config.sh.
# If you change a value in config.sh, change it here too.
#
# Sourced at the top of every R script:
#   source(file.path(.this_script_dir, "config.R"))
# =============================================================================

# ---------------------------------------------------------------------------
# PATHS
# .this_script_dir is set by each R script before sourcing this file,
# pointing to the scripts/ directory — config.R lives one level up.
# ---------------------------------------------------------------------------
JASPAR_FILE <- normalizePath(
  file.path(.this_script_dir, "..", "jaspar_tf_families.csv"),
  mustWork = FALSE
)

# ---------------------------------------------------------------------------
# COMPUTE RESOURCES
# ---------------------------------------------------------------------------
NCPUS <- as.integer(Sys.getenv("PBS_NCPUS", unset = "8"))

# ---------------------------------------------------------------------------
# STEP 03 — Fisher exact tests
# ---------------------------------------------------------------------------
FISHER_MIN_FAMILY_SIZE   <- 20     # min LTR regions per family to test
FISHER_MIN_HITS_IN_FAMILY <- 10    # min hits within family to test
FISHER_MIN_OUTSIDE_TOTAL <- 10     # min total regions outside family
BH_FDR_THRESHOLD         <- 0.05   # Benjamini-Hochberg FDR threshold

# ---------------------------------------------------------------------------
# STEP 05 — LTR vs control comparison
# ---------------------------------------------------------------------------
LTR_SPECIFIC_LOG2OR  <-  1    # log2(OR) threshold for LTR-specific category
LTR_DEPLETED_LOG2OR  <- -1    # log2(OR) threshold for depleted category

# ---------------------------------------------------------------------------
# STEP 07 — fgsea
# ---------------------------------------------------------------------------
GSEA_MIN_SIZE  <- 5        # min TFs per family to test
GSEA_MAX_SIZE  <- 200      # max TFs per family to test
GSEA_N_PERM    <- 100000   # number of permutations
GSEA_SEED      <- 42       # random seed for reproducibility

# ---------------------------------------------------------------------------
# STEP 08 — Plots
# ---------------------------------------------------------------------------
PLOT_Q_THRESH         <- 0.05     # significance threshold
PLOT_OR_THRESH        <- 2.0      # strong enrichment threshold
PLOT_Q_NETWORK        <- 0.0001   # stricter threshold for network plot
PLOT_TOP_N_TFS        <- 20L      # top TFs in heatmap and network
PLOT_MIN_HITS_TF      <- 10L      # min hits for TF to appear in heatmap
PLOT_TOP_N_GSEA_BAR   <- 15L      # max families in GSEA barplot
PLOT_TOP_N_GSEA_CURVES <- 4L      # enrichment curves shown

# ---------------------------------------------------------------------------
# STEP 09 — Report
# ---------------------------------------------------------------------------
REPORT_Q_THRESH              <- 0.05   # significance threshold
REPORT_OR_THRESH             <- 2.0    # strong enrichment threshold
REPORT_MIN_TESTED_FAMILIES   <- 3L     # min families for consistent TF table
REPORT_SMALL_FAMILY_WARN     <- 50L    # family size warning threshold
REPORT_DOMINANT_FAMILY_FRAC  <- 0.80   # dominant family warning threshold
