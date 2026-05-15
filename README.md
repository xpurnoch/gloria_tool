# GLORIA — Genome LTR Oriented Regulation Analysis
A bioinformatic pipeline for identifying transcription factor binding site (TFBS) enrichment
in LTR retrotransposon families across one or more plant genomes.

---

## Directory structure
```
gloria_tool/
├── setup.sh                   ← Setup conda environments and set root directory
├── config.sh                  ← All parameters and thresholds for Bash scripts
├── config.R                   ← All parameters and thresholds for R scripts
├── motifs.meme                ← MEME-format TF binding motif file
├── jaspar_tf_families.csv     ← JASPAR TF-to-family mapping file
├── genomes/                   ← Place genome FASTA files here (.fna / .fa / .fasta)
├── output/                    ← Results are written here
├── envs/
│   ├── dante_ltr.yml          ← Environment for DANTE, HMMER, bedtools, R
│   └── meme.yml               ← Environment for FIMO (MEME suite)
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
    ├── setup.pbs               ← PBS job script for creating conda environments
    └── run_all.pbs             ← PBS job script for running all steps
```

---

## First-time setup

Clone the repository and enter the project directory:

```bash
git clone https://github.com/xpurnoch/gloria_tool
cd gloria_tool
```

Run `setup.sh` from `gloria_tool/` directory to create conda environments and set root directory:

```bash
bash setup.sh
```

This step may take a while to complete (approximately 40 minutes). The terminal will be blocked
until setup is finished — do not close it.

> **Why does setup run on a compute node?**
> The `dante_ltr` environment includes the R package `GenomeInfoDbData`, which fails to install
> on MetaCentrum login nodes due to thread limits. Running setup on a compute node avoids this.

---

## Running the pipeline

### Place input files

Copy genome FASTA files into `genomes/` (`.fna`, `.fa`, or `.fasta` extension).

### Run the pipeline

Run `scripts/run_all.pbs` from `gloria_tool/` directory to submit job to the PBS cluster:

```bash
qsub -v "GENOMES=test.fna" scripts/run_all.pbs
```

For multiple genomes, use a space-separated list:

```bash
qsub -v "GENOMES=potato.fna tomato.fna" scripts/run_all.pbs
```

### Monitor the job

Use `qstat` to monitor the job status:

```bash
qstat -u $USER
```

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

All numerical thresholds are located in the config files. See `config.sh` and `config.R` for details.

Before submitting the job, you can open `scripts/run_all.pbs` and edit the `#PBS` header
depending on the number and size of input genomes:

- `ncpus` — number of available CPU cores (default: `8`)
- `mem` — total memory used (default: `200gb`)
- `walltime` — maximum runtime (default: `120:00:00`)
