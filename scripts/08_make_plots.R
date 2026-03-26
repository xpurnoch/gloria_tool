#!/usr/bin/env Rscript
# =============================================================================
# 08_make_plots.R
#
# Step 8/9 — Generate all analysis plots
#
# All parameters are read from config.sh via load_config.R.
# Edit config.sh to change thresholds — do not edit this script.
#
# Plots requiring control data (2, 7, 8, 9) are skipped when
# SKIP_CONTROL=1 or TFBS_random_controls_objects.RData is missing.
# Plot 3 is only generated when more than one genome is present.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(Biostrings)
  library(fgsea)
  library(cowplot)
  library(igraph)
  library(ggraph)
  library(tidyr)
})

.this_script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(trailingOnly=FALSE)[startsWith(commandArgs(trailingOnly=FALSE), "--file=")])))
source(file.path(.this_script_dir, "..", "config.R"))

cat("[STEP 8/9] Generating plots\n")
cat(sprintf("  Plot parameters from config.sh:\n"))
cat(sprintf("    PLOT_Q_THRESH=%s  PLOT_OR_THRESH=%s  PLOT_Q_NETWORK=%s\n",
            PLOT_Q_THRESH, PLOT_OR_THRESH, PLOT_Q_NETWORK))
cat(sprintf("    PLOT_TOP_N_TFS=%d  PLOT_MIN_HITS_TF=%d\n",
            PLOT_TOP_N_TFS, PLOT_MIN_HITS_TF))

ECDF_COLORS <- c("LTR" = "#3e2b6f", "RANDOM_GENOMIC" = "#f27070")

# =============================================================================
# make_clean_theme
# =============================================================================
make_clean_theme <- function() {
  theme_minimal(base_size = 14) +
    theme(plot.title   = element_text(face = "bold", hjust = 0.5),
          axis.title   = element_text(face = "bold"),
          axis.text    = element_text(color = "black"),
          panel.grid   = element_blank(),
          axis.line    = element_line(color = "black", linewidth = 0.5),
          legend.title = element_text(face = "bold"))
}

# =============================================================================
# load_data
# =============================================================================
load_data <- function() {
  load("TFBS_LTR_part1_objects.RData", envir = .GlobalEnv)

  skip_ctrl <- identical(Sys.getenv("SKIP_CONTROL"), "1")

  if (!skip_ctrl) {
    if (file.exists("TFBS_random_controls_objects.RData")) {
      load("TFBS_random_controls_objects.RData", envir = .GlobalEnv)
    } else {
      cat("  WARNING: TFBS_random_controls_objects.RData not found — control plots skipped\n")
      skip_ctrl <- TRUE
    }
  } else {
    cat("  INFO: SKIP_CONTROL=1 — plots depending on random control will be skipped\n")
  }
  skip_ctrl
}

# =============================================================================
# select_top_tfs
# =============================================================================
select_top_tfs <- function(fisher_df) {
  suppressWarnings({
    top_tfs <- fisher_df %>%
      filter(is.finite(odds_ratio), odds_ratio > 0,
             ltr_in_family_with_tf >= PLOT_MIN_HITS_TF,
             q_value < PLOT_Q_THRESH) %>%
      mutate(log_or = log2(odds_ratio)) %>%
      filter(is.finite(log_or)) %>%
      group_by(tf) %>%
      summarise(max_log2OR = max(log_or, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(max_log2OR)) %>%
      slice_head(n = PLOT_TOP_N_TFS) %>%
      pull(tf)
  })
  cat(sprintf("  Top TFs selected: %d\n", length(top_tfs)))
  top_tfs
}

# =============================================================================
# plot_heatmap
# =============================================================================
plot_heatmap <- function(fisher_df, top_tfs, clean_theme) {
  if (length(top_tfs) == 0) return(invisible(NULL))

  family_order <- fisher_df %>%
    group_by(family) %>%
    summarise(n_sig = sum(q_value < PLOT_Q_THRESH, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(n_sig)) %>% pull(family)

  heat_df <- fisher_df %>%
    filter(tf %in% top_tfs) %>%
    mutate(score = -log10(q_value + 1e-300))

  heat_df <- tidyr::complete(heat_df, tf, family) %>%
    mutate(family = factor(family, levels = family_order),
           family_lbl = gsub(".*_", "", as.character(family)))

  lbl_map <- heat_df %>% distinct(family, family_lbl) %>% { setNames(.$family_lbl, .$family) }

  ggsave("plots/heatmap_topTFs.png",
    ggplot(heat_df, aes(family, tf, fill = score)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_x_discrete(labels = lbl_map) +
      scale_fill_gradientn(colours = c("#3e2b6f","#bdb3d4","#f27070","#f1a04b"),
                           name = "-log10(q)", na.value = "grey92") +
      labs(x = "LTR family", y = "TF",
           title = sprintf("TF Enrichment Heatmap  (top %d TFs, q < %.2f)",
                           PLOT_TOP_N_TFS, PLOT_Q_THRESH)) +
      clean_theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)),
    width = 11, height = 8, dpi = 300)
  cat("  Plot saved: heatmap_topTFs.png\n")
}

# =============================================================================
# plot_violin_or_comparison
# =============================================================================
plot_violin_or_comparison <- function(fisher_df, genomic_fisher_df, clean_theme) {
  extract_log2or <- function(df, label) {
    if (is.null(df) || nrow(df) == 0) return(data.frame(dataset = label, log2OR = numeric(0)))
    df %>% filter(is.finite(odds_ratio), odds_ratio > 0) %>%
      transmute(dataset = label, log2OR = log2(odds_ratio))
  }
  or_comparison <- bind_rows(extract_log2or(fisher_df, "LTR"),
                              extract_log2or(genomic_fisher_df, "RANDOM_GENOMIC")) %>%
    mutate(dataset = factor(dataset, levels = c("LTR","RANDOM_GENOMIC")))
  if (nrow(or_comparison) == 0) return(invisible(NULL))

  ltr_vals  <- or_comparison$log2OR[or_comparison$dataset == "LTR"]
  geo_vals  <- or_comparison$log2OR[or_comparison$dataset == "RANDOM_GENOMIC"]
  p_geo     <- if (length(geo_vals) > 0) wilcox.test(ltr_vals, geo_vals)$p.value else NA
  cat(sprintf("  Wilcoxon LTR vs Genomic: p = %.2e\n", p_geo))

  ggsave("plots/violin_OR_comparison.png",
    ggplot(or_comparison, aes(dataset, log2OR, fill = dataset)) +
      geom_violin(alpha = 0.6, trim = FALSE, color = NA) +
      geom_boxplot(width = 0.2, outlier.alpha = 0.3, alpha = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      scale_fill_manual(values = ECDF_COLORS, labels = c("LTR","Random genomic")) +
      scale_x_discrete(labels = c("LTR","Genomic\nshuffle")) +
      labs(x = NULL, y = "Family Enrichment\nLog2(OR)",
           title = "Family Enrichment: LTR vs Genomic Control") +
      clean_theme + theme(legend.position = "none"),
    width = 7, height = 6, dpi = 300)
  cat("  Plot saved: violin_OR_comparison.png\n")
}

# =============================================================================
# plot_per_genome_ecdf
# =============================================================================
plot_per_genome_ecdf <- function(fisher_df, ltr, genomes, clean_theme) {
  family_genome_map <- ltr %>% distinct(family, genome) %>%
    group_by(family) %>% slice_head(n = 1) %>% ungroup()

  fisher_genome <- fisher_df %>%
    left_join(family_genome_map, by = "family") %>% filter(!is.na(genome)) %>%
    mutate(sig_plot = pmin(-log10(q_value + 1e-300),
                           quantile(-log10(q_value + 1e-300), 0.99, na.rm = TRUE)))

  genome_pal <- setNames(scales::hue_pal()(length(genomes)), genomes)

  ggsave("plots/per_genome_ecdf.png",
    ggplot(fisher_genome, aes(sig_plot, color = genome)) +
      stat_ecdf(linewidth = 1.2) +
      scale_color_manual(values = genome_pal) +
      labs(x = "-log10(q-value)", y = "Cumulative fraction",
           color = "Genome", title = "TFBS enrichment significance per genome") +
      clean_theme +
      theme(legend.position = "top",
            legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)),
    width = 10, height = 7, dpi = 300)
  cat("  Plot saved: per_genome_ecdf.png\n")
}

# =============================================================================
# plot_family_enrichment_summary
# =============================================================================
plot_family_enrichment_summary <- function(fisher_df, clean_theme) {
  family_bar <- fisher_df %>%
    group_by(family) %>%
    summarise(enriched = sum(q_value < PLOT_Q_THRESH & odds_ratio > PLOT_OR_THRESH, na.rm = TRUE),
              depleted = sum(q_value < PLOT_Q_THRESH & odds_ratio < 1,              na.rm = TRUE),
              .groups  = "drop") %>%
    mutate(family = factor(family, levels = family[order(enriched + depleted)])) %>%
    tidyr::pivot_longer(cols = c(enriched, depleted), names_to = "direction", values_to = "n") %>%
    mutate(n = ifelse(direction == "depleted", -n, n),
           direction = factor(direction, levels = c("enriched","depleted")))

  if (nrow(family_bar) == 0) return(invisible(NULL))

  family_bar <- family_bar %>%
    mutate(family_lbl = gsub(".*_", "", as.character(family)))
  lbl_map_fam <- family_bar %>% distinct(family, family_lbl) %>% { setNames(.$family_lbl, .$family) }

  ggsave("plots/per_family_enrichment_summary.png",
    ggplot(family_bar, aes(family, n, fill = direction)) +
      geom_col() + coord_flip() +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
      scale_x_discrete(labels = lbl_map_fam) +
      scale_fill_manual(values = c("enriched" = "#3e2b6f", "depleted" = "#f27070"),
                        labels = c(sprintf("Enriched (OR > %.0f)", PLOT_OR_THRESH),
                                   "Depleted (OR < 1)")) +
      labs(x = NULL, y = sprintf("Significant TF Associations (q < %.2f)", PLOT_Q_THRESH),
           fill = NULL, title = "Enriched vs Depleted TF Associations") +
      clean_theme + theme(legend.position = "top", axis.text.y = element_text(size = 9)),
    width = 10, height = max(5, n_distinct(fisher_df$family) * 0.5), dpi = 300)
  cat("  Plot saved: per_family_enrichment_summary.png\n")
}

# =============================================================================
# plot_gsea_results
# =============================================================================
plot_gsea_results <- function(fg, stats, pathways, clean_theme) {
  gsea_sig <- fg %>% filter(padj < PLOT_Q_THRESH) %>% arrange(desc(NES))
  if (nrow(gsea_sig) == 0) return(invisible(NULL))

  top_bar <- gsea_sig %>% arrange(desc(abs(NES))) %>% slice_head(n = PLOT_TOP_N_GSEA_BAR)
  ggsave("plots/GSEA_barplot_top_families.png",
    ggplot(top_bar, aes(reorder(pathway, NES), NES, fill = NES)) +
      geom_col(width = 0.7) + coord_flip() +
      scale_fill_gradientn(colours = c("#3e2b6f","#bdb3d4","#f27070","#f1a04b"), guide = "none") +
      labs(x = NULL, y = "NES",
           title = sprintf("Top Enriched TF Families\n(padj < %.2f)", PLOT_Q_THRESH)) +
      clean_theme,
    width = 8, height = max(4, nrow(top_bar) * 0.4), dpi = 300)
  cat("  Plot saved: GSEA_barplot_top_families.png\n")

  top_sets <- gsea_sig %>% arrange(desc(abs(NES))) %>%
    slice_head(n = PLOT_TOP_N_GSEA_CURVES) %>% pull(pathway)

  curve_plots <- lapply(top_sets, function(set_name) {
    nes_val   <- round(fg$NES[fg$pathway == set_name], 2)
    padj_val  <- fg$padj[fg$pathway == set_name]
    padj_text <- ifelse(padj_val < 0.001, "< 0.001", sprintf("= %.3f", padj_val))
    suppressWarnings(plotEnrichment(pathways[[set_name]], stats)) +
      labs(title = set_name, subtitle = paste0("NES = ", nes_val, ", q ", padj_text),
           x = "Rank", y = "Enrichment score") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
            plot.subtitle = element_text(hjust = 0.5, size = 9),
            panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
            axis.line = element_line(color = "black", linewidth = 0.5))
  })

  n_cols <- min(length(curve_plots), 2)
  n_rows <- ceiling(length(curve_plots) / n_cols)
  ggsave("plots/GSEA_enrichment_curves_top4.png",
         plot_grid(plotlist = curve_plots, ncol = n_cols, labels = "AUTO"),
         width = 10, height = max(6, n_rows * 3.5), dpi = 300)
  cat("  Plot saved: GSEA_enrichment_curves_top4.png\n")
}

# =============================================================================
# plot_gc_qc
# =============================================================================
plot_gc_qc <- function(clean_theme) {
  calc_gc <- function(fasta, label) {
    if (!file.exists(fasta)) return(data.frame(gc = numeric(0), dataset = character(0)))
    seqs <- readDNAStringSet(fasta)
    if (length(seqs) == 0) return(data.frame(gc = numeric(0), dataset = character(0)))
    data.frame(gc = rowSums(letterFrequency(seqs, c("G","C"), as.prob = TRUE)), dataset = label)
  }
  ltr_fa <- if (file.exists("all_ltr_combined.fa")) "all_ltr_combined.fa" else "LTR_5prime.fa"
  cat(sprintf("  GC QC using: %s\n", ltr_fa))

  gc_df <- bind_rows(calc_gc(ltr_fa, "LTR"), calc_gc("random_genomic.fa", "RANDOM_GENOMIC")) %>%
    mutate(dataset = factor(dataset, levels = c("LTR","RANDOM_GENOMIC")))
  if (nrow(gc_df) == 0) return(invisible(NULL))

  gc_ltr  <- gc_df$gc[gc_df$dataset == "LTR"]
  gc_geno <- gc_df$gc[gc_df$dataset == "RANDOM_GENOMIC"]
  p_gc    <- if (length(gc_geno) > 0) wilcox.test(gc_ltr, gc_geno)$p.value else NA
  cat(sprintf("  GC Wilcoxon LTR vs Genomic: p = %.2e\n", p_gc))

  ggsave("plots/QC_GC_violin.png",
    ggplot(gc_df, aes(dataset, gc, fill = dataset)) +
      geom_violin(alpha = 0.7, trim = FALSE, color = NA) +
      geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
      stat_summary(fun = median, geom = "point", size = 2.5) +
      scale_fill_manual(values = ECDF_COLORS, labels = c("LTR","Genomic shuffle")) +
      scale_x_discrete(labels = c("LTR","Genomic\nshuffle")) +
      labs(x = NULL, y = "GC fraction", title = "GC Content: LTR vs Genomic Control") +
      clean_theme + theme(legend.position = "none"),
    width = 7, height = 6, dpi = 400)

  ggsave("plots/QC_GC_density.png",
    ggplot(gc_df, aes(gc, fill = dataset, color = dataset)) +
      geom_density(alpha = 0.3, linewidth = 1) +
      scale_fill_manual(values  = ECDF_COLORS, labels = c("LTR","Genomic shuffle")) +
      scale_color_manual(values = ECDF_COLORS, labels = c("LTR","Genomic shuffle")) +
      labs(x = "GC fraction", y = "Density", fill = "Dataset", color = "Dataset",
           title = "GC Content Distribution") +
      clean_theme +
      theme(legend.position = "top", legend.direction = "horizontal",
            legend.title = element_blank()),
    width = 9, height = 6.5, dpi = 400)
  cat("  Plots saved: QC_GC_violin.png, QC_GC_density.png\n")
}

# =============================================================================
# plot_ltr_vs_ctrl_scatter
# =============================================================================
plot_ltr_vs_ctrl_scatter <- function(clean_theme) {
  comp_file <- "ltr_vs_ctrl_comparison.tsv"
  if (!file.exists(comp_file)) {
    cat("  [SKIP] ltr_vs_ctrl_scatter.png (file not found)\n")
    return(invisible(NULL))
  }
  comparison <- read.table(comp_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  scatter_colors <- c("LTR-specific" = "#2166ac", "Shared (possible artefact)" = "#d6604d",
                      "Depleted in LTR" = "#4dac26", "Not significant" = "grey70")

  ggsave("plots/ltr_vs_ctrl_scatter.png",
    ggplot(comparison, aes(x = log2OR_ctrl, y = log2OR_ltr, color = category)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = 0, color = "grey80") +
      geom_vline(xintercept = 0, color = "grey80") +
      geom_point(alpha = 0.6, size = 1.5,
                 position = position_jitter(width = 0.15, height = 0.15)) +
      scale_color_manual(values = scatter_colors) +
      scale_x_continuous(limits = c(-6, 8)) +
      scale_y_continuous(limits = c(-8, 22)) +
      labs(title = "Enrichment LTR vs Genomic",
           x = "Genomic control\nLog2(OR)", y = "LTR families\nLog2(OR)", color = NULL) +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)),
    width = 7, height = 7, dpi = 300)
  cat("  Plot saved: ltr_vs_ctrl_scatter.png\n")
}

# =============================================================================
# plot_regulatory_network
# =============================================================================
plot_regulatory_network <- function(fisher_df, top_tfs, clean_theme) {
  suppressWarnings({
    network_df <- fisher_df %>%
      filter(tf %in% top_tfs, q_value < PLOT_Q_NETWORK,
             odds_ratio > PLOT_OR_THRESH, is.finite(odds_ratio)) %>%
      select(family, tf)

    if (nrow(network_df) == 0) {
      cat(sprintf("  [SKIP] Network plot — no pairs passed q < %.4f & OR > %.0f\n",
                  PLOT_Q_NETWORK, PLOT_OR_THRESH))
      return(invisible(NULL))
    }

    g         <- graph_from_data_frame(network_df, directed = FALSE)
    V(g)$type <- ifelse(V(g)$name %in% network_df$family, "LTR family", "TF")
    V(g)$deg  <- degree(g)

    set.seed(GSEA_SEED)
    p_net <- ggraph(g, layout = "fr") +
      geom_edge_link(colour = "grey70", width = 0.6, alpha = 0.4) +
      geom_node_point(aes(color = type, size = deg)) +
      geom_node_text(aes(label = name), repel = TRUE, size = 3) +
      scale_color_manual(name = "Type",
                         values = c("LTR family" = "#1b9e77", "TF" = "#d95f02")) +
      scale_size(name = "Degree", range = c(3, 8)) +
      labs(title = sprintf("Core Regulatory Network\n(q < %.4f, OR > %.0f, Top %d TFs)",
                           PLOT_Q_NETWORK, PLOT_OR_THRESH, PLOT_TOP_N_TFS)) +
      clean_theme +
      theme(axis.title = element_blank(), axis.text  = element_blank(),
            axis.ticks = element_blank(), axis.line  = element_blank())

    ggsave("plots/core_regulatory_network_LTR_TF.png",
           plot = p_net, width = 12, height = 8, dpi = 400, bg = "white")
    cat("  Plot saved: core_regulatory_network_LTR_TF.png\n")
  })
}

# =============================================================================
# MAIN
# =============================================================================
dir.create("plots", showWarnings = FALSE)

SKIP_CONTROL <- load_data()
clean_theme  <- make_clean_theme()
genomes      <- sort(unique(ltr$genome))
multi_genome <- length(genomes) > 1

cat(sprintf("  Genomes in data : %s\n", paste(genomes, collapse = ", ")))

top_tfs <- select_top_tfs(fisher_df)

plot_heatmap(fisher_df, top_tfs, clean_theme)

if (!SKIP_CONTROL) {
  plot_violin_or_comparison(fisher_df, genomic_fisher_df, clean_theme)
} else {
  cat("  [SKIP] violin_OR_comparison.png (SKIP_CONTROL=1)\n")
}

if (multi_genome) plot_per_genome_ecdf(fisher_df, ltr, genomes, clean_theme)

plot_family_enrichment_summary(fisher_df, clean_theme)

if (file.exists("GSEA_objects.RData")) {
  load("GSEA_objects.RData")
  if (exists("fg") && nrow(fg) > 0) plot_gsea_results(fg, stats, pathways, clean_theme)
}

if (!SKIP_CONTROL) {
  plot_gc_qc(clean_theme)
  plot_ltr_vs_ctrl_scatter(clean_theme)
} else {
  cat("  [SKIP] QC_GC plots (SKIP_CONTROL=1)\n")
  cat("  [SKIP] ltr_vs_ctrl_scatter.png (SKIP_CONTROL=1)\n")
}

plot_regulatory_network(fisher_df, top_tfs, clean_theme)

cat("[STEP 8/9] DONE\n")
