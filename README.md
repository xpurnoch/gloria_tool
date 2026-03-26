# GLORIA — Genome LTR Oriented Regulation Analysis

A bioinformatic pipeline for identifying transcription factor binding site (TFBS) enrichment
in LTR retrotransposon families across one or more plant genomes.

---

## Directory structure

```
gloria_tool/
├── config.sh                  ← All parameters and thresholds (edit here!)
├── motifs.meme                ← MEME-format TF binding motif file
├── jaspar_tf_families.csv     ← JASPAR TF-to-family mapping, semicolon-separated
├── genomes/                   ← Place genome FASTA files here (.fna / .fa / .fasta)
├── output/                    ← Run results are written here
├── envs/
│   ├── dante_ltr.yml          ← Conda environment for DANTE, HMMER, bedtools, R
│   └── meme.yml               ← Conda environment for FIMO (MEME suite)
└── scripts/
    ├── 01_dante_pipeline.sh
    ├── 02_hmmer_solo_ltr.sh
    ├── 03_tfbs_ltr_enrichment.R
    ├── 04_random_control_fimo.sh
    ├── 05_random_control_analysis.R
    ├── 06_prepare_gsea.R
    ├── 07_run_gsea.R
    ├── 08_make_plots.R
    ├── 09_report.R
    └── run_all.pbs               ← PBS job submission script
```

---

## First-time setup

### 1. Create conda environments on MetaCentrum

```bash
module add mambaforge
mamba env create -f envs/dante_ltr.yml
mamba env create -f envs/meme.yml
```

### 2. Place input files
Copy genome FASTA files into `genomes/` (`.fna`, `.fa`, or `.fasta` extension)

### 3. Configure

Open `config.sh` and set `CONDA_ENV_DANTE_LTR` and `CONDA_ENV_MEME` to the
absolute paths of the conda environments you created in step 1:

```bash
mamba env list
```
Copy the absolute paths of the environments and set them in `config.sh`:

```bash
CONDA_ENV_DANTE_LTR="$HOME/.conda/envs/dante_ltr"
CONDA_ENV_MEME="$HOME/.conda/envs/meme"
```
### 4. Setup root directory

In scripts/run_all.pbs, set the `GLORIA_ROOT` variable to the absolute path of the
gloria_tool directory:

```bash
GLORIA_ROOT="/storage/brno2/home/NAME/work/gloria_tool"
```

---

## Running the pipeline

```bash
# Submit with default settings
qsub -v "GENOMES=potato.fna" scripts/run_all.pbs

```

PBS resources are set in the `#PBS` header of `run_all.pbs`. Adjust `ncpus`,
`mem`, `scratch_local`, and `walltime` to match your cluster and data size.

---

## Outputs

Each run creates a subdirectory in `output/` named after the genome file(s):

```
output/<genome_name>/
├── LTR_5prime.bed / .fa              Full LTR regions
├── solo_LTR.bed / .fa                Solo LTR regions
├── FIMO_LTR/fimo.tsv                 FIMO hits in full LTRs
├── FIMO_SOLO_LTR/fimo.tsv            FIMO hits in Solo LTRs
├── FIMO_RANDOM_GENOMIC/fimo.tsv      FIMO hits in controls (if not skipped)
├── TFBS_LTR_enrichment_results.tsv   Fisher test results (main output)
├── ltr_vs_ctrl_comparison.tsv        LTR vs control classification
├── GSEA_TF_families.tsv              fgsea results
├── pipeline_report.txt               Plain-text summary report
└── plots/                            All analysis plots (PNG)
```
---

## Customising parameters

All numerical thresholds live in `config.sh`. You never need to edit the scripts
themselves. Key parameters:

| Parameter | Default | Controls |
|-----------|---------|----------|
| `FIMO_LTR_PVALUE` | `1e-5` | FIMO p-value for LTR scanning |
| `MARKOV_ORDER` | `2` | Order of Markov background model |
| `MIN_HIT_LEN` | `80` | Minimum Solo LTR hit length (bp) |
| `MIN_HIT_SCORE` | `50` | Minimum nhmmer bit score |
| `FISHER_MIN_FAMILY_SIZE` | `20` | Minimum LTR regions per family for Fisher test |
| `FISHER_MIN_HITS_IN_FAMILY` | `10` | Minimum hits within family for Fisher test |
| `BH_FDR_THRESHOLD` | `0.05` | Benjamini-Hochberg FDR threshold |
| `BEDTOOLS_SHUFFLE_SEED` | `123` | Reproducibility seed for genomic shuffle |
| `GSEA_N_PERM` | `100000` | Number of fgsea permutations |

See `config.sh` for the full list with descriptions.
