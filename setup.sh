#!/bin/bash
# =============================================================================
# setup.sh
#
# First-time setup for GLORIA pipeline on MetaCentrum.
# Run this script from the gloria_tool root directory:
#   bash setup.sh
#
# This script submits setup.pbs to a compute node and waits for it to finish.
# =============================================================================

echo ""
echo "============================================================"
echo " GLORIA Pipeline Setup"
echo "============================================================"
echo ""

# Sanity check for the correct working directory
[[ -f "$(pwd)/config.sh" ]] || \
  { echo "[ERROR] config.sh not found. Run setup.sh from the gloria_tool root directory."; exit 1; }

[[ -f "$(pwd)/scripts/setup.pbs" ]] || \
  { echo "[ERROR] scripts/setup.pbs not found."; exit 1; }

echo "[INFO] Submitting setup job to compute node..."
echo "[INFO] Working directory: $(pwd)"
echo ""

spinner() {
  local FRAMES=('⠋' '⠙' '⠹' '⠸' '⠼' '⠴' '⠦' '⠧' '⠇' '⠏')
  local i=0
  while true; do
    printf "\r[INFO] Setup running... ${FRAMES[$i]} (this may take up to 40 minutes)"
    i=$(( (i + 1) % ${#FRAMES[@]} ))
    sleep 0.1
  done
}

# Start spinner in background
spinner &
SPINNER_PID=$!

qsub -W block=true -v GLORIA_ROOT="$(pwd)" scripts/setup.pbs >/dev/null 2>&1
EXIT_CODE=$?

# Stop spinner
kill $SPINNER_PID 2>/dev/null
printf "\r%*s\r" 60 ""

echo ""
if [[ $EXIT_CODE -eq 0 ]]; then
  echo "============================================================"
  echo " Setup complete"
  echo "============================================================"
  echo ""
  echo "  Next step: submit the pipeline with"
  echo "    qsub -v GENOMES=\"test.fna\" scripts/run_all.pbs"
  echo "============================================================"
else
  echo "[ERROR] Setup job failed (exit code $EXIT_CODE)"
  echo "  Check the log file: setup_gloria.o<job_id>"
  exit $EXIT_CODE
fi
