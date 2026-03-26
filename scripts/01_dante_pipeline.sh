#!/bin/bash
# =============================================================================
# 01_dante_pipeline.sh
#
# Step 1/9 — DANTE annotation, LTR region extraction, and FIMO scanning
#
# Usage:
#   bash 01_dante_pipeline.sh
#
# All parameters are read from config.sh (sourced automatically).
# Genome files are read from GENOMES_DIR, motifs from MOTIFS_FILE.
#
# Outputs (written to RUN_DIR, set by run_all.pbs):
#   LTR_5prime.bed          Merged BED of all full LTR regions (all genomes)
#   LTR_5prime.fa           FASTA sequences of all full LTR regions
#   merged_genomes.fasta    Concatenated prefixed genome sequences
#   ltr_background.txt      2nd-order Markov background model (from LTR seqs)
#   FIMO_LTR/fimo.tsv       FIMO motif hits in full LTR sequences
#   genome_<PREFIX>/        Per-genome working directories
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

MOTIFS="$SCRATCHDIR/motifs.meme"

# Build list of genome paths from GENOMES variable (basenames only,
# files were copied to SCRATCHDIR by run_all.pbs)
GENOMES_ARRAY=()
for G in $GENOMES; do
  GENOMES_ARRAY+=("$SCRATCHDIR/$G")
done

[[ ${#GENOMES_ARRAY[@]} -eq 0 ]] && { echo "[ERROR] No genomes in GENOMES variable"; exit 1; }
[[ ! -s "$MOTIFS" ]]             && { echo "[ERROR] Motifs file not found: $MOTIFS"; exit 1; }

echo "[STEP 1/9] DANTE + LTR extraction + FIMO (${#GENOMES_ARRAY[@]} genome(s))"

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV_DANTE_LTR" >/dev/null 2>&1

export TMPDIR="${RUN_DIR}/tmp"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
mkdir -p "$TMPDIR"

rm -f LTR_5prime.bed LTR_5prime.fa

# ==========================================================
# derive_genome_prefix
#
# Derives a safe filesystem prefix from a genome FASTA filename.
# Strips the directory path and all extensions, then replaces
# spaces, dots, slashes, and hyphens with underscores.
# ==========================================================
derive_genome_prefix() {
  local GENOME_FASTA="$1"
  local BASENAME PREFIX
  BASENAME=$(basename "$GENOME_FASTA")
  PREFIX="${BASENAME%%.*}"
  echo "$PREFIX" | tr ' ./-' '_'
}

# ==========================================================
# prefix_fasta_headers
#
# Rewrites every FASTA header in a genome file to include the
# genome prefix, preventing chromosome name collisions when
# sequences from multiple genomes are concatenated later.
# For example, ">Chr1" becomes ">prefix_Chr1".
# ==========================================================
prefix_fasta_headers() {
  local INPUT="$1"
  local PREFIX="$2"
  local OUTPUT="$3"

  awk -v pfx="$PREFIX" '
    /^>/ { print ">" pfx "_" substr($0,2); next }
    { print }
  ' "$INPUT" > "$OUTPUT"
  return 0
}

# ==========================================================
# run_dante
#
# Runs DANTE to annotate retrotransposon protein domains in a
# genome FASTA. Output is a GFF3 file with domain annotations
# used as input for DANTE_LTR in the next step.
# ==========================================================
run_dante() {
  local INPUT="$1"
  local OUTPUT="$2"
  local N_CORES="$3"

  dante \
    -q "$INPUT" \
    -o "$OUTPUT" \
    -c "$N_CORES" \
    >/dev/null 2>&1

  if [[ ! -s "$OUTPUT" ]]; then
    echo "[ERROR] DANTE failed for $INPUT"
    return 1
  fi
  return 0
}

# ==========================================================
# run_dante_ltr
#
# Runs DANTE_LTR to predict complete LTR retrotransposon element
# boundaries from DANTE domain annotations.
# ==========================================================
run_dante_ltr() {
  local DANTE_GFF="$1"
  local GENOME_FA="$2"
  local OUTPUT_PREFIX="$3"

  dante_ltr \
    -g "$DANTE_GFF" \
    -s "$GENOME_FA" \
    -o "$OUTPUT_PREFIX" \
    -M 1 \
    >/dev/null 2>&1

  if [[ ! -s "${OUTPUT_PREFIX}.gff3" ]]; then
    echo "[ERROR] DANTE_LTR failed for $GENOME_FA"
    return 1
  fi
  return 0
}

# ==========================================================
# extract_ltr_bed
#
# Extracts full LTR regions from a DANTE_LTR GFF3 file and
# writes them in BED6 format.
#
# BED col4 (name) format: FamilyName::chr:start-end
# Strand is in col6 only — bedtools getfasta -s appends
# "(+)" or "(-)" to FASTA headers automatically.
# ==========================================================
extract_ltr_bed() {
  local INPUT_GFF="$1"
  local OUTPUT_BED="$2"

  awk '$3=="long_terminal_repeat" && $9 ~ /LTR=[53]LTR/' \
    "$INPUT_GFF" \
  | awk 'BEGIN{OFS="\t"}{
      match($9, /Final_Classification=([^;]+)/, a)
      class  = a[1]
      gsub(/[^A-Za-z0-9_-]/, "_", class)
      chr    = $1
      start  = $4 - 1
      end    = $5
      strand = $7
      name   = class "::" chr ":" start "-" end
      print chr, start, end, name, ".", strand
    }' > "$OUTPUT_BED"
  return 0
}

# ==========================================================
# process_genome
#
# Orchestrates all per-genome steps for a single genome FASTA.
# Designed to be called in a background subshell (&).
# ==========================================================
process_genome() {
  local GENOME_FASTA="$1"
  local N_CORES="${NCPUS:-8}"

  [[ ! -f "$GENOME_FASTA" ]] && { echo "[ERROR] Genome not found: $GENOME_FASTA"; return 1; }

  local PREFIX
  PREFIX=$(derive_genome_prefix "$GENOME_FASTA")

  echo "  Processing genome: $PREFIX (PID $$)"

  local WORKDIR="genome_${PREFIX}"
  mkdir -p "$WORKDIR"

  local PREFIXED_FA="${WORKDIR}/${PREFIX}_prefixed.fasta"
  prefix_fasta_headers "$GENOME_FASTA" "$PREFIX" "$PREFIXED_FA"
  run_dante     "$PREFIXED_FA" "${WORKDIR}/dante_output.gff3" "$N_CORES" || return 1
  run_dante_ltr "${WORKDIR}/dante_output.gff3" "$PREFIXED_FA" "${WORKDIR}/dante_ltr_results" || return 1
  extract_ltr_bed "${WORKDIR}/dante_ltr_results.gff3" "${WORKDIR}/LTR_5prime.bed"

  echo "    LTRs found in $PREFIX: $(wc -l < "${WORKDIR}/LTR_5prime.bed")"
  return 0
}

export -f derive_genome_prefix
export -f prefix_fasta_headers
export -f run_dante
export -f run_dante_ltr
export -f extract_ltr_bed
export -f process_genome

# ==========================================================
# run_genomes_parallel
#
# Launches process_genome in the background for every genome
# in the provided list, then collects all exit codes.
# ==========================================================
run_genomes_parallel() {
  set +e

  local PIDS=()
  for GENOME_FASTA in "$@"; do
    process_genome "$GENOME_FASTA" &
    PIDS+=($!)
  done

  local FAILED=0
  for PID in "${PIDS[@]}"; do
    wait "$PID"
    local RC=$?
    [[ $RC -ne 0 ]] && { echo "[ERROR] Genome processing failed (PID $PID)"; FAILED=1; }
  done

  set -e
  [[ $FAILED -eq 1 ]] && return 1
  return 0
}

# ==========================================================
# merge_ltr_beds
#
# Concatenates per-genome LTR BED files into a single sorted,
# deduplicated LTR_5prime.bed.
# ==========================================================
merge_ltr_beds() {
  for GENOME_FASTA in "$@"; do
    local PREFIX
    PREFIX=$(derive_genome_prefix "$GENOME_FASTA")
    local BED="genome_${PREFIX}/LTR_5prime.bed"
    [[ -s "$BED" ]] && cat "$BED" >> LTR_5prime.bed
  done

  if [[ ! -s LTR_5prime.bed ]]; then
    echo "[ERROR] No LTR regions found across all genomes"
    return 1
  fi

  sort -u LTR_5prime.bed -o LTR_5prime.bed

  echo "  Total LTR regions: $(wc -l < LTR_5prime.bed)"
  echo "  Unique families:   $(awk '{if(index($4,"::")>0) print substr($4,1,index($4,"::")-1); \
    else print $4}' LTR_5prime.bed | sort -u | wc -l)"
  return 0
}

# ==========================================================
# build_merged_genome
#
# Concatenates all per-genome prefixed FASTA files into a single
# merged_genomes.fasta and builds a samtools faidx index.
# ==========================================================
build_merged_genome() {
  local MERGED="merged_genomes.fasta"
  rm -f "$MERGED"

  for GENOME_FASTA in "$@"; do
    local PREFIX
    PREFIX=$(derive_genome_prefix "$GENOME_FASTA")
    cat "genome_${PREFIX}/${PREFIX}_prefixed.fasta" >> "$MERGED"
  done

  samtools faidx "$MERGED" 2>/dev/null

  echo "  Merged genome: $MERGED ($(wc -l < "${MERGED}.fai") sequences)"
  return 0
}

# ==========================================================
# extract_ltr_sequences
#
# Extracts LTR sequences from the merged genome FASTA using
# bedtools getfasta in strand-aware mode (-s).
#
# Header format: FamilyName::chr:start-end(strand)
# This is the canonical ltr_id used throughout the pipeline.
# ==========================================================
extract_ltr_sequences() {
  bedtools getfasta \
    -fi merged_genomes.fasta \
    -bed LTR_5prime.bed \
    -fo LTR_5prime.fa \
    -s \
    -nameOnly \
    >/dev/null 2>&1

  if [[ ! -s LTR_5prime.fa ]]; then
    echo "[ERROR] FASTA extraction failed"
    return 1
  fi

  echo "  LTR sequences extracted: $(grep -c '^>' LTR_5prime.fa)"
  return 0
}

# ==========================================================
# build_background_model
#
# Computes a Markov background model from the FULL LTR sequences
# (LTR_5prime.fa) using fasta-get-markov (MEME suite).
#
# Using LTR sequences rather than the whole genome ensures the
# background reflects the nucleotide composition of the sequences
# actually scanned by FIMO, giving a fair per-family comparison.
#
# Model order is set by MARKOV_ORDER in config.sh (default: 2).
# ==========================================================
build_background_model() {
  source "${CONDA_BASE}/etc/profile.d/conda.sh"
  conda activate "$CONDA_ENV_MEME" >/dev/null 2>&1

  fasta-get-markov -m "$MARKOV_ORDER" LTR_5prime.fa > ltr_background.txt 2>/dev/null

  if [[ -s ltr_background.txt ]]; then
    echo "  Background model (order ${MARKOV_ORDER}, from LTR seqs): ltr_background.txt"
  else
    echo "  WARNING: background model empty — FIMO will use default"
  fi
  return 0
}

# ==========================================================
# split_motif_file  /  launch_fimo_chunks  /  wait_for_fimo_chunks
# merge_fimo_chunks  /  run_fimo_parallel
#
# Parallel FIMO orchestration — splits the motif file into
# N_CHUNKS chunks (round-robin), runs one FIMO process per chunk,
# then merges all per-chunk fimo.tsv files into a single TSV.
# set -e is disabled during background-job sections to allow
# proper exit-code collection via wait().
# ==========================================================
split_motif_file() {
  local MOTIFS_F="$1"
  local TMPBASE="$2"
  local N_CHUNKS="$3"

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

  FIMO_N_MOTIFS=$(ls "${TMPBASE}"/motif_*.txt 2>/dev/null | wc -l)

  if [[ "$FIMO_N_MOTIFS" -eq 0 ]]; then
    echo "[fimo_parallel] ERROR: no motifs found in $MOTIFS_F"
    return 1
  fi

  [[ "$N_CHUNKS" -gt "$FIMO_N_MOTIFS" ]] && N_CHUNKS="$FIMO_N_MOTIFS"

  local i chunk_id chunk_file
  for (( i=1; i<=FIMO_N_MOTIFS; i++ )); do
    chunk_id=$(( (i - 1) % N_CHUNKS + 1 ))
    chunk_file="${TMPBASE}/chunk_${chunk_id}.meme"
    if [[ ! -f "$chunk_file" ]]; then
      cat "$HEADER_FILE" > "$chunk_file"
    fi
    cat "${TMPBASE}/motif_${i}.txt" >> "$chunk_file"
  done

  FIMO_N_CHUNKS=$(ls "${TMPBASE}"/chunk_*.meme 2>/dev/null | wc -l)
  return 0
}

launch_fimo_chunks() {
  local N_CHUNKS="$1"
  local TMPBASE="$2"
  local SEQS="$3"
  shift 3
  local EXTRA_FLAGS=("$@")

  local chunk_id chunk_file CHUNK_OUT
  for (( chunk_id=1; chunk_id<=N_CHUNKS; chunk_id++ )); do
    chunk_file="${TMPBASE}/chunk_${chunk_id}.meme"
    [[ ! -s "$chunk_file" ]] && continue
    CHUNK_OUT="${TMPBASE}/fimo_out_${chunk_id}"
    mkdir -p "$CHUNK_OUT"
    (
      fimo \
        --oc "$CHUNK_OUT" \
        "${EXTRA_FLAGS[@]}" \
        "$chunk_file" \
        "$SEQS" \
        >/dev/null 2>&1
    ) &
    FIMO_PIDS+=($!)
  done
  return 0
}

wait_for_fimo_chunks() {
  FIMO_FAILED=0
  for PID in "${FIMO_PIDS[@]}"; do
    wait "$PID"
    local RC=$?
    [[ $RC -ne 0 ]] && FIMO_FAILED=$(( FIMO_FAILED + 1 ))
  done
  return 0
}

merge_fimo_chunks() {
  local N_CHUNKS="$1"
  local TMPBASE="$2"
  local MERGED_TSV="$3"

  local HEADER_WRITTEN=0
  rm -f "$MERGED_TSV"

  local chunk_id CHUNK_TSV
  for (( chunk_id=1; chunk_id<=N_CHUNKS; chunk_id++ )); do
    CHUNK_TSV="${TMPBASE}/fimo_out_${chunk_id}/fimo.tsv"
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
    return 0
  fi

  echo "[fimo_parallel] Merged hits: $(tail -n +2 "$MERGED_TSV" | wc -l) → ${MERGED_TSV}"
  return 0
}

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

  local FIMO_N_MOTIFS=0
  local FIMO_N_CHUNKS=0
  local FIMO_PIDS=()
  local FIMO_FAILED=0

  split_motif_file "$MOTIFS_F" "$TMPBASE" "$N_CHUNKS"
  if [[ $? -ne 0 ]]; then
    rm -rf "$TMPBASE"
    set -e; return 1
  fi

  echo "[fimo_parallel] Motifs: $FIMO_N_MOTIFS | Chunks: $FIMO_N_CHUNKS | Sequences: $(basename "$SEQS")"

  launch_fimo_chunks "$FIMO_N_CHUNKS" "$TMPBASE" "$SEQS" "${EXTRA_FLAGS[@]}"
  wait_for_fimo_chunks

  if [[ "$FIMO_FAILED" -eq "${#FIMO_PIDS[@]}" ]]; then
    echo "[fimo_parallel] ERROR: all FIMO chunks failed"
    rm -rf "$TMPBASE"
    set -e; return 1
  fi
  [[ "$FIMO_FAILED" -gt 0 ]] && \
    echo "[fimo_parallel] WARNING: $FIMO_FAILED / ${#FIMO_PIDS[@]} chunks failed — merging available results"

  merge_fimo_chunks "$FIMO_N_CHUNKS" "$TMPBASE" "${OUTDIR}/fimo.tsv"

  rm -rf "$TMPBASE"
  set -e; return 0
}

# ==========================================================
# run_fimo_ltr
#
# Scans full LTR sequences for TF binding motif occurrences.
# Uses ltr_background.txt (built from LTR seqs) if available.
# p-value threshold from config: FIMO_LTR_PVALUE.
# ==========================================================
run_fimo_ltr() {
  local MOTIFS_F="$1"
  local N_CHUNKS="$2"

  local BGFLAG=()
  [[ -s ltr_background.txt ]] && BGFLAG=("--bgfile" "ltr_background.txt")

  echo "  Running parallel FIMO (full LTR) — $N_CHUNKS chunks..."

  run_fimo_parallel \
    "$MOTIFS_F" \
    LTR_5prime.fa \
    FIMO_LTR \
    "$N_CHUNKS" \
    --no-pgc --thresh "$FIMO_LTR_PVALUE" "${BGFLAG[@]}"

  if [[ ! -s FIMO_LTR/fimo.tsv ]]; then
    echo "[ERROR] FIMO failed"
    return 1
  fi

  echo "  FIMO (full LTR): $(tail -n +2 FIMO_LTR/fimo.tsv | wc -l) motif hits"
  return 0
}

# =============================================================================
# MAIN
# =============================================================================
run_genomes_parallel "${GENOMES_ARRAY[@]}"
merge_ltr_beds       "${GENOMES_ARRAY[@]}"
build_merged_genome  "${GENOMES_ARRAY[@]}"
extract_ltr_sequences
build_background_model
run_fimo_ltr "$MOTIFS" "${NCPUS}"

echo "[STEP 1/9] Done"
