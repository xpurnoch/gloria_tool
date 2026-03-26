#!/bin/bash
# =============================================================================
# setup.sh — GLORIA one-time setup script
#
# Run this once after unpacking gloria_tool:
#   bash setup.sh
#
# What it does:
#   1. Creates conda environments from envs/*.yml
#   2. Writes correct conda env paths into config.sh
# =============================================================================
set -euo pipefail

module add mambaforge
GLORIA_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "============================================"
echo " GLORIA setup"
echo " Root: $GLORIA_ROOT"
echo "============================================"

# =============================================================================
# 1. Find conda/mamba
# =============================================================================
echo ""
echo "[1/2] Locating conda/mamba..."

if command -v mamba &>/dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &>/dev/null; then
    CONDA_CMD="conda"
else
    echo "[ERROR] Neither mamba nor conda found in PATH."
    echo "        Load the module first, e.g.: module add mambaforge"
    exit 1
fi

CONDA_BASE="$($CONDA_CMD info --base 2>/dev/null)"
[[ -z "$CONDA_BASE" ]] && { echo "[ERROR] Cannot determine conda base path"; exit 1; }
echo "  conda base : $CONDA_BASE"
echo "  command    : $CONDA_CMD"

# =============================================================================
# 2. Create conda environments
# =============================================================================
echo ""
echo "[2/2] Creating conda environments..."

source "${CONDA_BASE}/etc/profile.d/conda.sh"

for ENV_YML in "$GLORIA_ROOT/envs/"*.yml; do
    ENV_NAME=$(grep '^name:' "$ENV_YML" | awk '{print $2}')
    ENV_PATH="$HOME/.conda/envs/$ENV_NAME"

    if conda env list | grep -q "^$ENV_NAME "; then
        echo "  [SKIP] $ENV_NAME already exists — $(conda env list | grep "^$ENV_NAME " | awk '{print $NF}')"
    else
        echo "  Creating $ENV_NAME (this may take 10-20 minutes)..."
        $CONDA_CMD env create -f "$ENV_YML" -q
        echo "  [OK] $ENV_NAME created"
    fi

    # Get actual path (may differ from default if user has custom conda config)
    ENV_PATH=$(conda env list | grep "^$ENV_NAME " | awk '{print $NF}')
    echo "  Path: $ENV_PATH"

    # Write path into config.sh
    if [[ "$ENV_NAME" == "dante_ltr" ]]; then
        sed -i "s|CONDA_ENV_DANTE_LTR=.*|CONDA_ENV_DANTE_LTR=\"$ENV_PATH\"|" "$GLORIA_ROOT/config.sh"
    elif [[ "$ENV_NAME" == "meme" ]]; then
        sed -i "s|CONDA_ENV_MEME=.*|CONDA_ENV_MEME=\"$ENV_PATH\"|" "$GLORIA_ROOT/config.sh"
    fi
done

# =============================================================================
# Done
# =============================================================================
echo ""
echo "============================================"
echo " Setup complete"
echo "============================================"
echo ""
echo " Next steps:"
echo "   1. Copy genome files into:  $GLORIA_ROOT/genomes/"
echo "   2. Copy motif file to:      $GLORIA_ROOT/motifs.meme"
echo "   3. Copy JASPAR mapping to:  $GLORIA_ROOT/jaspar_tf_families.csv"
echo "   4. Submit a job:"
echo "      qsub -v GENOMES=\"genome.fna\" $GLORIA_ROOT/scripts/run_all.pbs"
echo ""
