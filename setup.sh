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
echo "[INFO] This may take up to 40 minutes. Do not close this terminal."
echo ""

qsub -W block=true -v GLORIA_ROOT="$(pwd)" scripts/setup.pbs

EXIT_CODE=$?

echo ""
if [[ $EXIT_CODE -eq 0 ]]; then
  echo "============================================================"
  echo " Setup complete"
  echo "============================================================"
  echo ""
  echo "  Next step: submit the pipeline with"
  echo "    qsub -v GENOMES=\"genome1.fna\" scripts/run_all.pbs"
  echo "============================================================"
else
  echo "[ERROR] Setup job failed (exit code $EXIT_CODE)"
  echo "  Check the log file: setup_gloria.o<job_id>"
  exit $EXIT_CODE
fi
