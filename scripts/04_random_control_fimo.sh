#!/bin/bash
# =============================================================================
# 04_random_control_fimo.sh
#
# Step 4/9 — Generate genomic shuffle controls and run FIMO on them
#
# Creates random genomic regions matching the size and chromosome distribution
# of all LTR regions (full + solo). Controls are used in step 05 to
# distinguish LTR-specific TF enrichment from genomic background signal.
#
# Usage:
#   bash 04_random_control_fimo.sh
#
# All parameters are read from config.sh (sourced automatically).
#
# Inputs (produced by steps 01 and 02):
#   LTR_5prime.bed / .fa        Full LTR coordinates and sequences
#   solo_LTR.bed / .fa          Solo LTR coordinates and sequences
#   merged_genomes.fasta[.fai]  Merged genome + index
#   ltr_background.txt          Markov background model (from LTR seqs)
#   genome_*/dante_ltr_results.gff3
#
# Outputs:
#   all_ltr_combined.bed          Merged BED of all LTR regions (full + solo)
#   all_ltr_combined.fa           Merged FASTA of all LTR sequences
#   all_ltr_mask.bed              Exclusion mask for genomic shuffle
#   random_genomic.bed            Shuffled control region coordinates
#   random_genomic.fa             Shuffled control region sequences
#   FIMO_RANDOM_GENOMIC/fimo.tsv  FIMO motif hits in control sequences
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

MOTIFS="$MOTIFS_FILE"

echo "[STEP 4/9] Generate controls + FIMO on controls"

# ==========================================================
# INPUT CHECKS
# ==========================================================
[[ ! -s LTR_5prime.bed           ]] && { echo "[ERROR] LTR_5prime.bed missing";           exit 10; }
[[ ! -s LTR_5prime.fa            ]] && { echo "[ERROR] LTR_5prime.fa missing";            exit 11; }
[[ ! -s merged_genomes.fasta.fai ]] && { echo "[ERROR] merged_genomes.fasta.fai missing"; exit 12; }
[[ ! -s "$MOTIFS"                ]] && { echo "[ERROR] Motifs file missing: $MOTIFS";     exit 13; }

export PATH="$CONDA_ENV_DANTE_LTR/bin:$PATH"

# ==========================================================
# run_fimo_parallel  (same logic as steps 01 and 02)
# ==========================================================
run_fimo_parallel() {
  set +e

  local MOTIFS_F="$1"
  local SEQS="$2"
  local OUTDIR="$3"
  local N_CHUNKS="$4"
  shift 4
  local EXTRA_FLAGS=("$@")

  if [[ ! -s "$MOTIFS_F" ]]; then
    echo "[fimo_parallel] ERROR: motif file missing: $MOTIFS_F"
    set -e; return 1
  fi
  if [[ ! -s "$SEQS" ]]; then
    echo "[fimo_parallel] ERROR: sequence file missing: $SEQS"
    set -e; return 1
  fi
  [[ "$N_CHUNKS" -lt 1 ]] && N_CHUNKS=1

  mkdir -p "$OUTDIR"

  local TMPBASE
  TMPBASE=$(mktemp -d "${SCRATCHDIR}/fimo_split_XXXXXX")

  local HEADER_FILE="${TMPBASE}/meme_header.txt"
  awk '/^MOTIF /{exit} {print}' "$MOTIFS_F" > "$HEADER_FILE"

  awk -v base="${TMPBASE}/motif_" '
    /^MOTIF / {
      if (out) close(out)
      motif_idx++
      out = base motif_idx ".txt"
    }
    out { print > out }
  ' "$MOTIFS_F"

  local N_MOTIFS
  N_MOTIFS=$(ls "${TMPBASE}"/motif_*.txt 2>/dev/null | wc -l)

  if [[ "$N_MOTIFS" -eq 0 ]]; then
    echo "[fimo_parallel] ERROR: no motifs found in $MOTIFS_F"
    rm -rf "$TMPBASE"
    set -e; return 1
  fi

  [[ "$N_CHUNKS" -gt "$N_MOTIFS" ]] && N_CHUNKS="$N_MOTIFS"

  echo "[fimo_parallel] Motifs: $N_MOTIFS | Chunks: $N_CHUNKS | Sequences: $(basename "$SEQS")"

  local i chunk_id chunk_file
  for (( i=1; i<=N_MOTIFS; i++ )); do
    chunk_id=$(( (i - 1) % N_CHUNKS + 1 ))
    chunk_file="${TMPBASE}/chunk_${chunk_id}.meme"
    if [[ ! -f "$chunk_file" ]]; then
      cat "$HEADER_FILE" > "$chunk_file"
    fi
    cat "${TMPBASE}/motif_${i}.txt" >> "$chunk_file"
  done

  local PIDS=()
  local CHUNK_DIRS=()
  for (( chunk_id=1; chunk_id<=N_CHUNKS; chunk_id++ )); do
    chunk_file="${TMPBASE}/chunk_${chunk_id}.meme"
    [[ ! -s "$chunk_file" ]] && continue
    local CHUNK_OUT="${TMPBASE}/fimo_out_${chunk_id}"
    mkdir -p "$CHUNK_OUT"
    CHUNK_DIRS+=("$CHUNK_OUT")
    (
      fimo \
        --oc "$CHUNK_OUT" \
        "${EXTRA_FLAGS[@]}" \
        "$chunk_file" \
        "$SEQS" \
        >/dev/null 2>&1
    ) &
    PIDS+=($!)
  done

  local FAILED=0
  for PID in "${PIDS[@]}"; do
    wait "$PID"
    local RC=$?
    [[ $RC -ne 0 ]] && FAILED=$(( FAILED + 1 ))
  done

  if [[ "$FAILED" -eq "${#PIDS[@]}" ]]; then
    echo "[fimo_parallel] ERROR: all FIMO chunks failed"
    rm -rf "$TMPBASE"
    set -e; return 1
  fi
  [[ "$FAILED" -gt 0 ]] && \
    echo "[fimo_parallel] WARNING: $FAILED / ${#PIDS[@]} chunks failed — merging available results"

  local MERGED_TSV="${OUTDIR}/fimo.tsv"
  local HEADER_WRITTEN=0
  rm -f "$MERGED_TSV"

  for CHUNK_OUT in "${CHUNK_DIRS[@]}"; do
    local CHUNK_TSV="${CHUNK_OUT}/fimo.tsv"
    [[ ! -s "$CHUNK_TSV" ]] && continue
    if [[ "$HEADER_WRITTEN" -eq 0 ]]; then
      grep -v '^#' "$CHUNK_TSV" | head -1 >> "$MERGED_TSV"
      HEADER_WRITTEN=1
    fi
    grep -v '^#' "$CHUNK_TSV" | tail -n +2 >> "$MERGED_TSV"
  done

  if [[ ! -s "$MERGED_TSV" ]]; then
    echo "[fimo_parallel] WARNING: no motif hits found — writing empty TSV"
    echo -e "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence" \
      > "$MERGED_TSV"
    rm -rf "$TMPBASE"
    set -e; return 0
  fi

  echo "[fimo_parallel] Merged hits: $(tail -n +2 "$MERGED_TSV" | wc -l) → ${MERGED_TSV}"
  rm -rf "$TMPBASE"
  set -e; return 0
}

# ==========================================================
# combine_ltr_beds
#
# Merges full LTR and Solo LTR BED + FASTA files into single
# combined files used as the template for genomic shuffling.
# ==========================================================
combine_ltr_beds() {
  touch solo_LTR.bed solo_LTR.fa

  cat LTR_5prime.bed solo_LTR.bed \
    | sort -k1,1 -k2,2n \
    > all_ltr_combined.bed

  if [[ ! -s all_ltr_combined.bed ]]; then
    echo "[ERROR] Combined LTR BED is empty"
    return 1
  fi

  cat LTR_5prime.fa solo_LTR.fa > all_ltr_combined.fa

  if [[ ! -s all_ltr_combined.fa ]]; then
    echo "[ERROR] Combined LTR FASTA is empty"
    return 1
  fi

  N_FULL=$(wc -l < LTR_5prime.bed)
  N_SOLO=$(wc -l < solo_LTR.bed)
  N_ALL=$(wc -l  < all_ltr_combined.bed)

  echo "  Full LTRs   : $N_FULL"
  echo "  Solo LTRs   : $N_SOLO"
  echo "  Combined    : $N_ALL regions for controls"
  return 0
}

# ==========================================================
# build_ltr_exclusion_mask
#
# Builds a BED mask of all known LTR-related genomic features
# to be excluded when placing random control regions.
# ==========================================================
build_ltr_exclusion_mask() {
  local TMP_MASK
  TMP_MASK=$(mktemp)

  for GFF in genome_*/dante_ltr_results.gff3; do
    [[ ! -s "$GFF" ]] && continue
    awk 'BEGIN{OFS="\t"}
      !/^#/ && $3 ~ /LTR|target_site/ {
        print $1, $4-1, $5
      }' "$GFF"
  done >> "$TMP_MASK"

  [[ -s solo_LTR.bed ]] && \
    awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' solo_LTR.bed >> "$TMP_MASK"

  sort -k1,1 -k2,2n "$TMP_MASK" | bedtools merge -i stdin > all_ltr_mask.bed
  rm -f "$TMP_MASK"

  echo "  Exclusion mask: $(wc -l < all_ltr_mask.bed) merged regions"
  return 0
}

# ==========================================================
# shuffle_genomic_controls
#
# Creates random control regions by shuffling the combined LTR
# BED within each chromosome, excluding known LTR regions.
# Random seed from config: BEDTOOLS_SHUFFLE_SEED.
# ==========================================================
shuffle_genomic_controls() {
  echo "  Generating genomic shuffle control (seed=$BEDTOOLS_SHUFFLE_SEED)..."

  bedtools shuffle \
    -i all_ltr_combined.bed \
    -g merged_genomes.fasta.fai \
    -chrom \
    -seed "$BEDTOOLS_SHUFFLE_SEED" \
    -excl all_ltr_mask.bed \
    2>/dev/null | \
    awk 'BEGIN{OFS="\t"}{
      id = $1 "_" $2 "_" $3 "_" $6
      print $1, $2, $3, id, ".", $6
    }' > random_genomic.bed

  if [[ ! -s random_genomic.bed ]]; then
    echo "[ERROR] Genomic shuffle failed"
    return 1
  fi

  echo "  Controls: $(wc -l < random_genomic.bed) genomic regions"
  return 0
}

# ==========================================================
# extract_control_sequences
#
# Extracts FASTA sequences for shuffled control regions.
# -s is intentionally NOT used — control IDs do not follow
# the canonical ltr_id format and need no strand suffix.
# ==========================================================
extract_control_sequences() {
  bedtools getfasta \
    -fi merged_genomes.fasta \
    -bed random_genomic.bed \
    -fo random_genomic.fa \
    -nameOnly \
    >/dev/null 2>&1

  if [[ ! -s random_genomic.fa ]]; then
    echo "[ERROR] Genomic FASTA extraction failed"
    return 1
  fi

  echo "  Control sequences extracted: $(grep -c '^>' random_genomic.fa)"
  return 0
}

# ==========================================================
# run_fimo_controls
#
# Scans genomic shuffle control sequences for TF binding motifs.
# Intentionally uses NO background model (FIMO default) — controls
# represent random genomic sequence so the background should reflect
# general genomic composition, not LTR-specific composition.
# p-value threshold from config: FIMO_CTRL_PVALUE.
# ==========================================================
run_fimo_controls() {
  local MOTIFS_F="$1"
  local N_CHUNKS="$2"

  echo "  Running parallel FIMO on controls — $N_CHUNKS chunks..."

  export PATH="$CONDA_ENV_MEME/bin:$PATH"

  run_fimo_parallel \
    "$MOTIFS_F" \
    random_genomic.fa \
    FIMO_RANDOM_GENOMIC \
    "$N_CHUNKS" \
    --no-pgc --thresh "$FIMO_CTRL_PVALUE"

  if [[ ! -s FIMO_RANDOM_GENOMIC/fimo.tsv ]]; then
    echo "[ERROR] FIMO genomic output empty"
    return 1
  fi

  echo "  Genomic: $(tail -n +2 FIMO_RANDOM_GENOMIC/fimo.tsv | wc -l) hits"
  return 0
}

# =============================================================================
# MAIN
# =============================================================================
combine_ltr_beds
build_ltr_exclusion_mask
shuffle_genomic_controls
extract_control_sequences
run_fimo_controls "$MOTIFS" "${NCPUS}"

echo "[STEP 4/9] Done"

