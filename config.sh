#!/bin/bash
# =============================================================================
# config.sh — GLORIA pipeline configuration (bash scripts only)
#
# Parameters used by R scripts (03, 05, 06, 07, 08, 09) are in scripts/config.R
#
# This file is sourced automatically by all bash scripts:
#   source "$SCRIPT_DIR/../config.sh"
# =============================================================================
GLORIA_ROOT="/storage/brno2/home/xpurnoch/work/gloria_tool"

# =============================================================================
# PATHS — derived from GLORIA_ROOT, no manual editing needed
# =============================================================================
SCRIPTS_DIR="$GLORIA_ROOT/scripts"
ENVS_DIR="$GLORIA_ROOT/envs"
GENOMES_DIR="$GLORIA_ROOT/genomes"
OUTPUT_DIR="$GLORIA_ROOT/output"
MOTIFS_FILE="$GLORIA_ROOT/motifs.meme"
JASPAR_FILE="$GLORIA_ROOT/jaspar_tf_families.csv"

# =============================================================================
# CONDA — resolved automatically, no manual editing needed
# Run `mamba info --base` to verify CONDA_BASE on your system
# =============================================================================
CONDA_BASE="$(mamba info --base 2>/dev/null)"
if [[ -z "$CONDA_BASE" ]]; then
  echo "[config.sh] ERROR: Cannot find conda base." >&2
  exit 1
fi

# Paths to conda environments — update if you installed them elsewhere
# Confirm with: mamba env list
CONDA_ENV_DANTE_LTR="/storage/brno2/home/xpurnoch/.conda/envs/dante_ltr"
CONDA_ENV_MEME="/storage/brno2/home/xpurnoch/.conda/envs/meme"

# =============================================================================
# COMPUTE RESOURCES
# =============================================================================

# Number of CPU cores — overridden automatically by PBS_NCPUS under PBS
NCPUS="${PBS_NCPUS:-8}"

# =============================================================================
# STEP 01 — DANTE + LTR extraction + FIMO
# =============================================================================

# p-value threshold for FIMO motif scanning of full LTR sequences
FIMO_LTR_PVALUE="1e-5"

# Order of the Markov background model built from LTR sequences
MARKOV_ORDER=2

# =============================================================================
# STEP 02 — HMMER Solo LTR detection
# =============================================================================

# Minimum Solo LTR hit length in bp
MIN_HIT_LEN=80

# Minimum nhmmer bit score
MIN_HIT_SCORE=50

# Minimum sequences per family to build an HMM profile
MIN_SEQS_FOR_PROFILE=1

# p-value threshold for FIMO on Solo LTR sequences
FIMO_SOLO_PVALUE="1e-5"

# nhmmer E-value thresholds
NHMMER_EVALUE="1e-10"
NHMMER_INCEVALUE="1e-10"

# Multiplier for per-family Solo LTR length cap (relative to max full LTR length)
SOLO_LTR_MAX_LEN_MULTIPLIER=1.5

# Overlap fraction above which a hit is excluded as part of a full element
SOLO_LTR_MASK_OVERLAP=0.5

# =============================================================================
# STEP 04 — Random genomic shuffle controls
# =============================================================================

# p-value threshold for FIMO on control sequences
FIMO_CTRL_PVALUE="1e-5"

# Random seed for bedtools shuffle (reproducibility)
BEDTOOLS_SHUFFLE_SEED=123
