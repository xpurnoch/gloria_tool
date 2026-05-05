# =============================================================================
# config.R — GLORIA pipeline parameters for R scripts
# =============================================================================

JASPAR_FILE <- normalizePath(
  file.path(.this_script_dir, "..", "jaspar_tf_families.csv"),
  mustWork = FALSE
)

# =============================================================================
# COMPUTE RESOURCES
# =============================================================================
NCPUS <- as.integer(Sys.getenv("PBS_NCPUS", unset = "8"))

# =============================================================================
# STEP 03 — Fisher exact tests
# =============================================================================
FISHER_MIN_FAMILY_SIZE   <- 20     # min LTR regions per family to test
FISHER_MIN_HITS_IN_FAMILY <- 10    # min hits within family to test
FISHER_MIN_OUTSIDE_TOTAL <- 10     # min total regions outside family
BH_FDR_THRESHOLD         <- 0.05   # Benjamini-Hochberg FDR threshold

# =============================================================================
# STEP 05 — LTR vs control comparison
# =============================================================================
LTR_SPECIFIC_LOG2OR  <-  1.5    # log2(OR) threshold for LTR-specific category
LTR_DEPLETED_LOG2OR  <- -1.5    # log2(OR) threshold for depleted category

# =============================================================================
# STEP 07 — fgsea
# =============================================================================
GSEA_MIN_SIZE  <- 5        # min TFs per family to test
GSEA_MAX_SIZE  <- 200      # max TFs per family to test
GSEA_N_PERM    <- 100000   # number of permutations
GSEA_SEED      <- 42       # random seed for reproducibility


# =============================================================================
# PLOT AND REPORT PARAMETERS — visual settings only, no need to change these
# =============================================================================
PLOT_Q_THRESH          <- 0.05     # significance threshold for all plots
PLOT_OR_THRESH         <- 2.0      # strong enrichment threshold (OR > this)
PLOT_Q_NETWORK         <- 0.0001   # stricter threshold for network plot
PLOT_TOP_N_TFS         <- 20L      # top TFs shown in heatmap and network
PLOT_MIN_HITS_TF       <- 10L      # min hits for TF to appear in heatmap
PLOT_TOP_N_GSEA_BAR    <- 15L      # max families shown in GSEA barplot
PLOT_TOP_N_GSEA_CURVES <- 4L       # enrichment curves shown in GSEA panel
PLOT_DPI               <- 300L     # resolution for all saved plots
PLOT_SCATTER_XLIM      <- c(-6, 8)    # x-axis limits for LTR vs ctrl scatter
PLOT_SCATTER_YLIM      <- c(-8, 22)   # y-axis limits for LTR vs ctrl scatter
PLOT_GSEA_NCOLS        <- 2L       # columns in GSEA enrichment curves panel

REPORT_Q_THRESH             <- 0.05   # significance threshold used in report
REPORT_OR_THRESH            <- 2.0    # strong enrichment threshold (OR > this)
REPORT_MIN_TESTED_FAMILIES  <- 3L     # min families tested for consistent TF table
REPORT_SMALL_FAMILY_WARN    <- 50L    # warn if family has fewer LTR regions than this
REPORT_DOMINANT_FAMILY_FRAC <- 0.80   # warn if one family exceeds this fraction of all LTRs
REPORT_TOP_N_ASSOCIATIONS   <- 20L    # rows in "Top N strongest associations" table
REPORT_TOP_N_CONSISTENT_TFS <- 20L    # rows in "Top N consistent TFs" table
REPORT_TOP_N_LTR_SPECIFIC   <- 10L    # rows in "Top N LTR-specific" ctrl section
REPORT_TOP_N_SHARED         <- 5L     # rows in "Top N shared artefacts" ctrl section
REPORT_TOP_N_GSEA_LEADING   <- 5L     # leading edge families printed in GSEA section
