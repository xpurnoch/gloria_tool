#!/bin/bash
# =============================================================================
# 02_hmmer_solo_ltr.sh
#
# Step 2/9 — HMMER-based Solo LTR detection and FIMO scanning
#
# Solo LTRs are remnants of complete LTR retrotransposons that lost their
# internal sequence via homologous recombination between the 5' and 3' LTR.
#
# Inputs (from step 01):
#   LTR_5prime.fa             Full LTR sequences (used to build HMM profiles)
#   LTR_5prime.bed            Full LTR coordinates (used to compute length caps)
#   merged_genomes.fasta      Merged genome sequences
#   merged_genomes.fasta.fai  Samtools index
#   ltr_background.txt        Markov background model (built from LTR seqs)
#   genome_*/dante_ltr_results.gff3  Per-genome DANTE_LTR annotations
#
# Outputs:
#   solo_LTR.bed              Solo LTR coordinates (BED6)
#   solo_LTR.fa               Solo LTR sequences (FASTA)
#   FIMO_SOLO_LTR/fimo.tsv    FIMO motif hits in Solo LTR sequences
#   family_fastas/            Per-family FASTA files
#   hmmer_profiles/           HMM profiles and nhmmer result tables
#   full_rte_mask.bed         Merged retrotransposon body mask
#   family_max_len.txt        Per-family maximum Solo LTR length caps
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

MOTIFS="$MOTIFS_FILE"

echo "[STEP 2/9] HMMER Solo LTR detection + FIMO"

# ==========================================================
# INPUT CHECKS
# ==========================================================
[[ ! -s LTR_5prime.fa            ]] && { echo "[ERROR] LTR_5prime.fa missing";            exit 20; }
[[ ! -s LTR_5prime.bed           ]] && { echo "[ERROR] LTR_5prime.bed missing";           exit 21; }
[[ ! -s merged_genomes.fasta     ]] && { echo "[ERROR] merged_genomes.fasta missing";     exit 22; }
[[ ! -s merged_genomes.fasta.fai ]] && { echo "[ERROR] merged_genomes.fasta.fai missing"; exit 23; }
[[ ! -s "$MOTIFS"                ]] && { echo "[ERROR] Motifs file missing: $MOTIFS";     exit 24; }

N_GFF=$(ls genome_*/dante_ltr_results.gff3 2>/dev/null | wc -l)
[[ $N_GFF -eq 0 ]] && { echo "[ERROR] No dante_ltr_results.gff3 found in genome_* dirs"; exit 25; }

export PATH="$CONDA_ENV_DANTE_LTR/bin:$PATH"

mkdir -p hmmer_profiles family_fastas

N_CORES="${NCPUS}"
echo "  Using $N_CORES cores"
echo "  Parameters:"
echo "    MIN_HIT_LEN=$MIN_HIT_LEN"
echo "    MIN_HIT_SCORE=$MIN_HIT_SCORE"
echo "    MIN_SEQS_FOR_PROFILE=$MIN_SEQS_FOR_PROFILE"
echo "    NHMMER_EVALUE=$NHMMER_EVALUE"
echo "    SOLO_LTR_MAX_LEN_MULTIPLIER=$SOLO_LTR_MAX_LEN_MULTIPLIER"
echo "    SOLO_LTR_MASK_OVERLAP=$SOLO_LTR_MASK_OVERLAP"

# ==========================================================
# run_fimo_parallel
#
# Splits a MEME motif file into N_CHUNKS chunks,
# runs one FIMO process per chunk in parallel, then merges all
# per-chunk fimo.tsv files into a single output TSV.
# Returns 0 even when no hits are found (empty-but-valid TSV).
# Returns 1 only on hard failures.
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

  for (( chunk_id=1; chunk_id<=N_CHUNKS; chunk_id++ )); do
    local CHUNK_TSV="${TMPBASE}/fimo_out_${chunk_id}/fimo.tsv"
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
# build_rte_mask
#
# Builds a merged BED mask of all known retrotransposon bodies
# from per-genome DANTE_LTR GFF3 files. Used to exclude full
# elements from Solo LTR hit candidates.
# ==========================================================
build_rte_mask() {
  echo "  Building full retrotransposon body mask from DANTE_LTR GFF3..."

  local TMP_MASK
  TMP_MASK=$(mktemp)

  for GFF in genome_*/dante_ltr_results.gff3; do
    [[ ! -s "$GFF" ]] && continue
    awk 'BEGIN{OFS="\t"}
      !/^#/ && $3 ~ /LTR|transposable_element|target_site/ {
        print $1, $4-1, $5
      }' "$GFF"
  done >> "$TMP_MASK"

  sort -k1,1 -k2,2n "$TMP_MASK" | bedtools merge -i stdin > full_rte_mask.bed
  rm -f "$TMP_MASK"

  N_MASKED=$(wc -l < full_rte_mask.bed)
  echo "  Mask regions (merged element bodies): $N_MASKED"

  if [[ $N_MASKED -eq 0 ]]; then
    echo "[ERROR] Full retrotransposon mask is empty"
    return 1
  fi
  return 0
}

# ==========================================================
# split_ltr_by_family
#
# Splits LTR_5prime.fa into one FASTA file per LTR family,
# written to family_fastas/. Family name = prefix before "::".
# ==========================================================
split_ltr_by_family() {
  echo "  Splitting sequences per family..."

  rm -f family_fastas/*.fa

  awk '
    /^>/ {
      header = substr($0, 2)
      if (index(header, "::") > 0)
        fam = substr(header, 1, index(header, "::") - 1)
      else
        fam = header
      gsub(/[^A-Za-z0-9_-]/, "_", fam)
      outfile = "family_fastas/" fam ".fa"
    }
    { print > outfile }
  ' LTR_5prime.fa

  N_FAMILIES=$(ls family_fastas/*.fa 2>/dev/null | wc -l)
  echo "  Total families: $N_FAMILIES"

  if [[ $N_FAMILIES -eq 0 ]]; then
    echo "[ERROR] No family FASTA files created"
    return 1
  fi
  return 0
}

# ==========================================================
# build_hmm_profile
#
# Builds a single HMM profile for one LTR family FASTA.
# Single-sequence families use --singlemx (no alignment needed).
# Multi-sequence families are aligned with MAFFT first.
# MIN_SEQS is read from MIN_SEQS_FOR_PROFILE in config.sh.
# ==========================================================
build_hmm_profile() {
  local FA="$1"
  local MIN_SEQS="$2"
  local FAMILY
  FAMILY=$(basename "$FA" .fa)
  local ALN="hmmer_profiles/${FAMILY}.aln"
  local HMM="hmmer_profiles/${FAMILY}.hmm"
  local N_SEQS
  N_SEQS=$(grep -c '^>' "$FA" 2>/dev/null || N_SEQS=0)

  if [[ $N_SEQS -lt $MIN_SEQS ]]; then
    echo "    [SKIP] $FAMILY: only $N_SEQS sequence(s) — below MIN_SEQS_FOR_PROFILE=$MIN_SEQS"
    return 0
  fi

  if [[ $N_SEQS -eq 1 ]]; then
    hmmbuild --singlemx --dna --cpu 1 -n "$FAMILY" "$HMM" "$FA" >/dev/null 2>&1
    if [[ -s "$HMM" ]]; then
      sed -i "s/^NAME.*/NAME  $FAMILY/" "$HMM"
    fi
  else
    if [[ $N_SEQS -gt 500 ]]; then
      mafft --auto --quiet --parttree --thread 1 "$FA" > "$ALN" 2>/dev/null
    else
      mafft --auto --quiet --thread 1 "$FA" > "$ALN" 2>/dev/null
    fi

    local ALN_CLEAN="${ALN}.clean"
    awk '
      /^>/ {
        if (header != "" && seq ~ /[ACGTNacgtn]/)
          { print header; print seq_raw }
        header=$0; seq=""; seq_raw=""
        next
      }
      { seq=seq $0; seq_raw=seq_raw"\n"$0 }
      END {
        if (header != "" && seq ~ /[ACGTNacgtn]/)
          { print header; print seq_raw }
      }
    ' "$ALN" > "$ALN_CLEAN" 2>/dev/null

    local N_CLEAN
    N_CLEAN=$(grep -c '^>' "$ALN_CLEAN" 2>/dev/null || N_CLEAN=0)

    if [[ "$N_CLEAN" -eq 0 ]]; then
        echo "    [WARN] All sequences gap-only in alignment for $FAMILY -- using --singlemx fallback"
        hmmbuild --singlemx --dna --cpu 1 -n "$FAMILY" "$HMM" "$FA" >/dev/null 2>&1
        if ! hmmstat "$HMM" >/dev/null 2>&1; then
            echo "    [INFO] HMM profile invalid for $FAMILY — skipping"
            rm -f "$HMM"
            return 0
        fi
    elif [[ "$N_CLEAN" -eq 1 ]]; then
      hmmbuild --singlemx --dna --cpu 1 -n "$FAMILY" "$HMM" "$ALN_CLEAN" >/dev/null 2>&1
    else
      hmmbuild --dna --cpu 1 -n "$FAMILY" "$HMM" "$ALN_CLEAN" >/dev/null 2>&1
    fi
  fi

  if [[ ! -s "$HMM" ]]; then
    echo "    [INFO] HMM profile skipped for $FAMILY — alignment insufficient"
    return 0
  fi
  echo "    Profile built: $FAMILY ($N_SEQS seqs)"
  return 0
}

export -f build_hmm_profile

# ==========================================================
# build_all_hmm_profiles
#
# Builds HMM profiles for all families in family_fastas/
# running build_hmm_profile in parallel using xargs -P.
# ==========================================================
build_all_hmm_profiles() {
  local MIN_SEQS="$1"
  local N_JOBS="$2"

  echo "  Building HMM profiles (min. $MIN_SEQS sequences per family)..."

  ls family_fastas/*.fa | \
    xargs -P "$N_JOBS" -I{} bash -c 'build_hmm_profile "$@"' _ {} "$MIN_SEQS"

  N_PROFILES=$(ls hmmer_profiles/*.hmm 2>/dev/null | wc -l)
  echo "  HMM profiles built: $N_PROFILES / $N_FAMILIES"

  if [[ $N_PROFILES -eq 0 ]]; then
    echo "[ERROR] No HMM profiles built"
    return 1
  fi
  return 0
}

# ==========================================================
# run_nhmmer_profile
#
# Runs nhmmer for a single HMM profile against the merged genome.
# E-value thresholds set by NHMMER_EVALUE / NHMMER_INCEVALUE in config.sh.
# ==========================================================
run_nhmmer_profile() {
  local HMM="$1"
  local GENOME="$2"
  local PROFILE_NUM="$3"
  local PROFILE_TOTAL="$4"
  local N_CORES="$5"

  local FAMILY
  FAMILY=$(basename "$HMM" .hmm)
  local OUTFILE="hmmer_profiles/${FAMILY}.tbl"

  echo "  [$PROFILE_NUM/$PROFILE_TOTAL] nhmmer: $FAMILY ..."

  nhmmer \
    --cpu "$N_CORES" \
    --tblout "$OUTFILE" \
    -E "$NHMMER_EVALUE" \
    --incE "$NHMMER_INCEVALUE" \
    --noali \
    "$HMM" \
    "$GENOME" \
    >/dev/null 2>&1
  local NHMMER_RC=$?

  if [[ $NHMMER_RC -ne 0 ]]; then
    echo "    [WARN] nhmmer failed for $FAMILY (exit $NHMMER_RC) — skipping"
    return 0
  fi
  if [[ ! -s "$OUTFILE" ]]; then
    echo "    [WARN] nhmmer produced no hits for $FAMILY"
    return 0
  fi

  local N_HITS
  N_HITS=$(grep -vc '^#' "$OUTFILE" 2>/dev/null || echo 0)
  echo "    Done: $FAMILY ($N_HITS hits)"

  grep -v '^#' "$OUTFILE" >> hmmer_raw_hits.tbl || true
  return 0
}

# ==========================================================
# run_all_nhmmer
#
# Iterates over all HMM profiles and calls run_nhmmer_profile
# sequentially (each profile uses all cores via --cpu).
# Fixes nhmmer alphabet detection on low-entropy sequences
# by replacing the first 4 bases of every sequence with GATC.
# ==========================================================
run_all_nhmmer() {
  local GENOME="$1"
  local N_CORES="$2"

  local GENOME_FIXED="${SCRATCHDIR}/merged_nhmmer_fixed.fa"
  awk '
    /^>/ { print; first=1; next }
    first { print "GATC" substr($0, 5); first=0; next }
    { print }
  ' "$GENOME" > "$GENOME_FIXED"
  echo "  Genome fixed for nhmmer alphabet detection"

  echo "  Running nhmmer sequentially, each profile using $N_CORES cores..."

  rm -f hmmer_raw_hits.tbl

  local N_PROFILES_TOTAL PROFILE_NUM
  N_PROFILES_TOTAL=$(ls hmmer_profiles/*.hmm 2>/dev/null | wc -l)
  PROFILE_NUM=0

  for HMM in hmmer_profiles/*.hmm; do
    PROFILE_NUM=$(( PROFILE_NUM + 1 ))
    run_nhmmer_profile "$HMM" "$GENOME_FIXED" "$PROFILE_NUM" "$N_PROFILES_TOTAL" "$N_CORES"
  done

  if [[ ! -s hmmer_raw_hits.tbl ]]; then
    echo "[ERROR] nhmmer produced no output"
    return 1
  fi

  N_RAW=$(wc -l < hmmer_raw_hits.tbl)
  echo "  Raw nhmmer hits: $N_RAW"
  return 0
}

# ==========================================================
# convert_nhmmer_to_bed
#
# Converts nhmmer tblout to BED6 with 0-based coordinates.
# BED col4: FamilyName::chr:start-end  (strand in col6 only).
# ==========================================================
convert_nhmmer_to_bed() {
  echo "  Converting nhmmer hits to BED..."

  awk 'BEGIN{OFS="\t"}{
    chr    = $1
    family = $3
    start  = $9  - 1
    end    = $10
    score  = $14
    strand = $12
    if (start > end) { tmp=start; start=end; end=tmp }
    if (start < 0) start = 0
    if (end > start) {
      gsub(/[^A-Za-z0-9_-]/, "_", family)
      name = family "::" chr ":" start "-" end
      print chr, start, end, name, score, strand
    }
  }' hmmer_raw_hits.tbl | sort -k1,1 -k2,2n > hmmer_hits_raw.bed

  if [[ ! -s hmmer_hits_raw.bed ]]; then
    echo "[ERROR] No valid nhmmer BED entries"
    return 1
  fi
  return 0
}

# ==========================================================
# filter_by_length_and_score
#
# Retains only hits meeting MIN_HIT_LEN and MIN_HIT_SCORE
# thresholds (set in config.sh).
# ==========================================================
filter_by_length_and_score() {
  local MIN_LEN="$1"
  local MIN_SCORE="$2"

  echo "  Filtering by length (>= ${MIN_LEN} bp) and score (>= ${MIN_SCORE})..."

  awk -v min_len="$MIN_LEN" -v min_score="$MIN_SCORE" \
    'BEGIN{OFS="\t"} ($3 - $2) >= min_len && $5 >= min_score' \
    hmmer_hits_raw.bed > hmmer_hits_filtered.bed

  N_FILTERED_SCORE=$(wc -l < hmmer_hits_filtered.bed)
  echo "  Hits after length+score filter: $N_FILTERED_SCORE (was $N_RAW raw)"
  return 0
}

# ==========================================================
# filter_by_family_max_length
#
# Removes hits exceeding SOLO_LTR_MAX_LEN_MULTIPLIER × the
# maximum full LTR length per family (set in config.sh).
# ==========================================================
filter_by_family_max_length() {
  echo "  Filtering by per-family maximum length..."

  awk -v mult="$SOLO_LTR_MAX_LEN_MULTIPLIER" 'BEGIN{OFS="\t"}{
    name = $4
    if (index(name, "::") > 0)
      fam = substr(name, 1, index(name, "::") - 1)
    else
      fam = name
    len = $3 - $2
    if (len > max[fam]) max[fam] = len
  } END {
    for (fam in max) print fam, int(max[fam] * mult)
  }' LTR_5prime.bed > family_max_len.txt

  echo "  Per-family max Solo LTR length (${SOLO_LTR_MAX_LEN_MULTIPLIER}x Full LTR max):"
  awk '{printf "    %-55s %d bp\n", $1, $2}' family_max_len.txt

  awk 'BEGIN{OFS="\t"}
    NR==FNR { maxlen[$1] = $2; next }
    {
      name = $4
      if (index(name, "::") > 0)
        fam = substr(name, 1, index(name, "::") - 1)
      else
        fam = name
      len = $3 - $2
      if (fam in maxlen && len <= maxlen[fam])
        print
      else if (!(fam in maxlen))
        print
    }
  ' family_max_len.txt hmmer_hits_filtered.bed > hmmer_hits_lenfiltered.bed

  N_LENFILTERED=$(wc -l < hmmer_hits_lenfiltered.bed)
  echo "  Hits after per-family length cap: $N_LENFILTERED (was $N_FILTERED_SCORE)"

  cp hmmer_hits_lenfiltered.bed hmmer_hits_filtered.bed
  return 0
}

# ==========================================================
# filter_against_rte_mask
#
# Removes nhmmer hits overlapping known retrotransposon bodies
# by more than SOLO_LTR_MASK_OVERLAP fraction (config.sh).
# ==========================================================
filter_against_rte_mask() {
  local MASK="$1"

  echo "  Filtering against full retrotransposon body mask..."

  bedtools intersect \
    -a hmmer_hits_filtered.bed \
    -b "$MASK" \
    -v \
    -f "$SOLO_LTR_MASK_OVERLAP" \
    > solo_LTR.bed

  if [[ ! -s solo_LTR.bed ]]; then
    echo "  WARNING: No Solo LTRs found after filtering"
    touch solo_LTR.bed solo_LTR.fa
    return 1
  fi

  N_SOLO=$(wc -l < solo_LTR.bed)
  echo "  Solo LTRs after filtering: $N_SOLO"
  return 0
}

# ==========================================================
# extract_solo_ltr_sequences
#
# Extracts Solo LTR sequences from the merged genome using
# bedtools getfasta with strand-aware mode (-s).
# Header format: FamilyName::chr:start-end(strand)
# ==========================================================
extract_solo_ltr_sequences() {
  bedtools getfasta \
    -fi merged_genomes.fasta \
    -bed solo_LTR.bed \
    -fo solo_LTR.fa \
    -nameOnly \
    -s \
    >/dev/null 2>&1

  if [[ ! -s solo_LTR.fa ]]; then
    echo "[ERROR] Solo LTR FASTA extraction failed"
    return 1
  fi

  echo "  Solo LTR sequences extracted: $(grep -c '^>' solo_LTR.fa)"
  return 0
}

# ==========================================================
# run_fimo_solo_ltr
#
# Scans Solo LTR sequences for TF binding motif occurrences.
# Uses ltr_background.txt (built from LTR seqs in step 01).
# p-value threshold from config: FIMO_SOLO_PVALUE.
# ==========================================================
run_fimo_solo_ltr() {
  local MOTIFS_F="$1"
  local N_CHUNKS="$2"

  echo "  Running parallel FIMO on Solo LTRs — $N_CHUNKS chunks..."

  local BGFLAG=()
  [[ -s ltr_background.txt ]] && BGFLAG=("--bgfile" "ltr_background.txt")

  run_fimo_parallel \
    "$MOTIFS_F" \
    solo_LTR.fa \
    FIMO_SOLO_LTR \
    "$N_CHUNKS" \
    --no-pgc --thresh "$FIMO_SOLO_PVALUE" "${BGFLAG[@]}"

  if [[ ! -f FIMO_SOLO_LTR/fimo.tsv ]]; then
    echo "  WARNING: FIMO_SOLO_LTR/fimo.tsv missing — creating fallback empty TSV"
    mkdir -p FIMO_SOLO_LTR
    head -1 FIMO_LTR/fimo.tsv > FIMO_SOLO_LTR/fimo.tsv 2>/dev/null || \
      echo -e "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence" \
      > FIMO_SOLO_LTR/fimo.tsv
  fi

  echo "  FIMO (Solo LTR): $(tail -n +2 FIMO_SOLO_LTR/fimo.tsv | wc -l) motif hits"
  return 0
}

# ==========================================================
# print_summary
# ==========================================================
print_summary() {
  echo ""
  echo "  === Solo LTR summary ==="
  echo "  Cores used              : $N_CORES"
  echo "  Min seqs for profile    : $MIN_SEQS_FOR_PROFILE"
  echo "  Profiles built          : $N_PROFILES"
  echo "  Mask regions            : $N_MASKED"
  echo "  Raw nhmmer hits         : $N_RAW"
  echo "  After length+score      : $N_FILTERED_SCORE"
  echo "  After per-family len cap: $N_LENFILTERED"
  echo "  Solo LTRs (final)       : $N_SOLO"

  awk '{
    name = $4
    if (index(name, "::") > 0)
      fam = substr(name, 1, index(name, "::") - 1)
    else
      fam = name
    print fam
  }' solo_LTR.bed | sort | uniq -c | sort -rn | \
    awk '{printf "    %-60s %d\n", $2, $1}'
  return 0
}

# =============================================================================
# MAIN
# =============================================================================
build_rte_mask
split_ltr_by_family
build_all_hmm_profiles "$MIN_SEQS_FOR_PROFILE" "$N_CORES"
run_all_nhmmer "$PWD/merged_genomes.fasta" "$N_CORES"
convert_nhmmer_to_bed
filter_by_length_and_score "$MIN_HIT_LEN" "$MIN_HIT_SCORE"

if [[ ! -s hmmer_hits_filtered.bed ]]; then
  touch solo_LTR.bed solo_LTR.fa
  echo "[STEP 2/9] Done (0 hits after filtering)"
  exit 0
fi

filter_by_family_max_length

filter_against_rte_mask full_rte_mask.bed || {
  echo "[STEP 2/9] Done (0 Solo LTRs)"
  exit 0
}

extract_solo_ltr_sequences

# Switch to meme environment for FIMO.
export PATH="$CONDA_ENV_MEME/bin:$PATH"
run_fimo_solo_ltr "$MOTIFS" "$N_CORES"
print_summary

echo "[STEP 2/9] Done"
