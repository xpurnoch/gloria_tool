#!/bin/bash
# =============================================================================
# setup.sh
#
# First-time setup for GLORIA pipeline on MetaCentrum.
# Run this script from the gloria_tool root directory:
#   bash setup.sh
#
# What this script does:
#   1. Detects GLORIA_ROOT from current directory
#   2. Creates conda environments (skips if already exist)
#   3. Detects conda env paths from mamba env list
#   4. Writes paths into pipeline_config.sh
#   5. Sets GLORIA_ROOT in scripts/run_all.pbs
# =============================================================================
set -euo pipefail

RESET='\033[0m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'

info()    { echo -e "${GREEN}[INFO]${RESET} $*"; }
warn()    { echo -e "${YELLOW}[WARN]${RESET} $*"; }
error()   { echo -e "${RED}[ERROR]${RESET} $*"; exit 1; }

# =============================================================================
# STEP 1: Detect GLORIA_ROOT
# =============================================================================
echo ""
echo "============================================================"
echo " GLORIA Pipeline Setup"
echo "============================================================"
echo ""

GLORIA_ROOT="$(pwd)"
info "GLORIA_ROOT detected: $GLORIA_ROOT"

# Sanity check — are we in the right directory?
[[ -f "$GLORIA_ROOT/config.sh" ]] || \
  error "config.sh not found. Run setup.sh from the gloria_tool root directory."

# =============================================================================
# STEP 2: Load mambaforge and create conda environments
# =============================================================================
info "Loading mambaforge..."
module add mambaforge 2>/dev/null || warn "Could not load mambaforge module — assuming mamba is already in PATH"

if ! command -v mamba &>/dev/null; then
  error "mamba not found. Load mambaforge first: module add mambaforge"
fi

info "Creating conda environments..."

for ENV_YML in envs/dante_ltr.yml envs/meme.yml; do
  [[ -f "$ENV_YML" ]] || error "Environment file not found: $ENV_YML"

  ENV_NAME=$(grep '^name:' "$ENV_YML" | awk '{print $2}')
  [[ -z "$ENV_NAME" ]] && error "Could not parse env name from $ENV_YML"

  if mamba env list | grep -q "^${ENV_NAME}\s"; then
    info "Environment '$ENV_NAME' already exists — skipping creation"
  else
    info "Creating environment: $ENV_NAME"
    mamba env create -f "$ENV_YML"
    info "Environment '$ENV_NAME' created"
  fi
done

# =============================================================================
# STEP 3: Detect conda env paths
# =============================================================================
info "Detecting conda environment paths..."

get_env_path() {
  local ENV_NAME="$1"
  local PATH_FOUND
  PATH_FOUND=$(mamba env list | grep "^${ENV_NAME}\s" | awk '{print $NF}')
  if [[ -z "$PATH_FOUND" ]]; then
    error "Could not find path for conda environment: $ENV_NAME"
  fi
  echo "$PATH_FOUND"
}

DANTE_LTR_PATH=$(get_env_path "dante_ltr")
MEME_PATH=$(get_env_path "meme")

info "dante_ltr env path: $DANTE_LTR_PATH"
info "meme env path:      $MEME_PATH"

# =============================================================================
# STEP 4: Write paths into config.sh
# =============================================================================
info "Writing paths into config.sh..."

CONFIG="$GLORIA_ROOT/config.sh"

# Add or update CONDA_ENV_DANTE_LTR
if grep -q "^CONDA_ENV_DANTE_LTR=" "$CONFIG"; then
  sed -i "s|^CONDA_ENV_DANTE_LTR=.*|CONDA_ENV_DANTE_LTR=\"${DANTE_LTR_PATH}\"|" "$CONFIG"
else
  echo "" >> "$CONFIG"
  echo "# Conda environment paths (set by setup.sh)" >> "$CONFIG"
  echo "CONDA_ENV_DANTE_LTR=\"${DANTE_LTR_PATH}\"" >> "$CONFIG"
fi

# Add or update CONDA_ENV_MEME
if grep -q "^CONDA_ENV_MEME=" "$CONFIG"; then
  sed -i "s|^CONDA_ENV_MEME=.*|CONDA_ENV_MEME=\"${MEME_PATH}\"|" "$CONFIG"
else
  echo "CONDA_ENV_MEME=\"${MEME_PATH}\"" >> "$CONFIG"
fi

# Add or update GLORIA_ROOT
if grep -q "^GLORIA_ROOT=" "$CONFIG"; then
    sed -i "s|^GLORIA_ROOT=.*|GLORIA_ROOT=\"${GLORIA_ROOT}\"|" "$CONFIG"
  else
  echo "GLORIA_ROOT="${GLORIA_ROOT}"" >> "$CONFIG"
fi

info "config.sh updated"

# =============================================================================
# STEP 5: Set GLORIA_ROOT in scripts/run_all.pbs
# =============================================================================
PBS_SCRIPT="$GLORIA_ROOT/scripts/run_all.pbs"

if [[ -f "$PBS_SCRIPT" ]]; then
  info "Setting GLORIA_ROOT in scripts/run_all.pbs..."
  if grep -q "^GLORIA_ROOT=" "$PBS_SCRIPT"; then
    sed -i "s|^GLORIA_ROOT=.*|GLORIA_ROOT=\"${GLORIA_ROOT}\"|" "$PBS_SCRIPT"
  else
    sed -i "s|^# =* MAIN|GLORIA_ROOT=\"${GLORIA_ROOT}\"\n\n# === MAIN|" "$PBS_SCRIPT"
  fi
  info "scripts/run_all.pbs updated"
else
  warn "scripts/run_all.pbs not found — skipping GLORIA_ROOT setup"
fi

# =============================================================================
# SUMMARY
# =============================================================================
echo ""
echo "============================================================"
echo " Setup complete"
echo "============================================================"
echo "  GLORIA_ROOT      : $GLORIA_ROOT"
echo "  dante_ltr env    : $DANTE_LTR_PATH"
echo "  meme env         : $MEME_PATH"
echo "  Config file      : $CONFIG"
echo ""
echo "  Next step: submit a job with"
echo "    qsub -v GENOMES=\"genome1.fna\" scripts/run_all.pbs"
echo "============================================================"
