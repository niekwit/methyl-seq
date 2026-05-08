# Snakemake workflow: `methyl-seq`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.25.5-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/niekwit/methyl-seq/workflows/Tests/badge.svg?branch=main)](https://github.com/niekwit/methyl-seq/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

A Snakemake workflow for EM-seq / bisulfite sequencing (WGBS) data analysis. Starting from FASTQ files, it performs read trimming, alignment with Bismark, methylation extraction, QC, and downstream analysis including CpG methylation visualisation and optional DMR calling.

## Installation

**1. Install Snakemake (≥ 8.25.5) via conda/mamba:**

```bash
mamba create -n snakemake -c conda-forge -c bioconda snakemake=8.25.2
conda activate snakemake
```

**2. Clone this repository:**

```bash
git clone https://github.com/niekwit/methyl-seq.git
cd methyl-seq
```

Conda environments for all tools are created automatically on first run — no manual software installation is needed.

## Usage

### 1. Set up your project directory

Create a new directory for your analysis and copy or symlink the workflow into it:

```
my_analysis/
├── reads/                  # FASTQ files (see naming convention below)
├── bed/                    # BED files for regions of interest
├── config/
│   ├── config.yaml
│   └── samples.csv
└── workflow -> /path/to/methyl-seq/workflow
```

Alternatively, run directly from the cloned repository by placing `reads/`, `bed/`, and `config/` inside it.

### 2. FASTQ naming convention

| Library type | Expected filename pattern                              |
| ------------ | ------------------------------------------------------ |
| Paired-end   | `{sample}_R1_001.fastq.gz`, `{sample}_R2_001.fastq.gz` |
| Single-end   | `{sample}.fastq.gz`                                    |

Paired-end vs. single-end is detected automatically from the filenames.

### 3. Configure `config/samples.csv`

Two required columns. Sample and condition names must contain only alphanumerics and underscores.

```csv
sample,condition
WT_1,WT
WT_2,WT
KO_1,KO
KO_2,KO
```

### 4. Configure `config/config.yaml`

```yaml
genome: hg38 # hg19 | hg38 | mm38 | mm39
ensembl_genome_build: 114 # Ensembl release number
temp_dir: /tmp # use "/local/" on Cambridge HPC

# Trim Galore clipping arguments (adjust for EM-seq or WGBS library prep)
trim_galore_args: "--clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10"

bismark:
  align: "" # extra arguments for bismark alignment
  deduplicate: ""
  extract: ""
  coverage: ""
  report: ""

deeptools:
  bigwig_summary:
    binSize: 10000
    extra: ""
  plotPCA:
    extra: ""

# CpG methylation boxplots
boxplot:
  plot: True
  cpg_n: 50 # CpGs per probe
  min_reads: 5 # minimum read coverage to keep a CpG
  regions:
    # Genomic regions (BED format) to plot. Whole genome is always included.
    CpG_islands: "config/annotations/hg38_CpG_islands.bed"
    promoters: "config/annotations/hg38_promoters.bed"

# DMR analysis (optional)
DMR:
  run: True
  reference_condition: WT # must match a condition in samples.csv
  tile_size: 1000
  step_size: 1000
  min_per_group: 2
  difference_threshold: 25
  qvalue_threshold: 0.01
```

### 5. Run the workflow

```bash
# Dry run (check what will be executed)
snakemake --snakefile workflow/Snakefile --use-conda -n

# Run locally
snakemake --snakefile workflow/Snakefile --use-conda -c <threads>

# Run with Singularity container
snakemake --snakefile workflow/Snakefile --use-conda --use-singularity -c <threads>

# Run on a SLURM cluster (adjust profile as needed)
snakemake --snakefile workflow/Snakefile --use-conda --profile slurm
```

---

## Re-running boxplots only

If alignment is already complete and you only want to regenerate boxplots (e.g. with different regions), use the standalone sub-workflow:

```bash
snakemake --snakefile workflow_boxplot_only/Snakefile --use-conda -c <threads>
```

This starts from `resources/filtered_cpg_probes.bed` and `results/bed/CpG_merged_{condition}.bed`, which must already exist from a previous full run.

---

## Expected output

```
results/
├── multiqc/
│   └── multiqc_report.html          # Trimming, alignment, and QC summary
├── bigwig/
│   ├── {condition}.bw               # Average CpG methylation BigWig per condition
│   └── coverage/
│       └── {condition}.bw           # Average read coverage BigWig per condition
├── deeptools/
│   └── PCA.tab
├── plots/
│   ├── PCA.pdf                      # PCA of methylation profiles
│   ├── scree.pdf
│   ├── methylation_conversion_rate.pdf
│   ├── methylation_conversion_rate.csv
│   └── boxplots.pdf                 # CpG methylation boxplots (if boxplot.plot: True)
├── preseq/
│   └── library_complexity_summary.txt
└── dmrs/                            # Only produced if DMR.run: True
    ├── hypermethylated_DMRs.bed
    ├── hypomethylated_DMRs.bed
    ├── hypermethylated_DMRs_annotated.tab
    ├── hypomethylated_DMRs_annotated.tab
    ├── all_methylation_tiles.rds
    ├── differential_methylation_tiles.rds
    ├── significant_differential_methylation_tiles.rds
    └── plots/dmrs/
        ├── DMR_distance_to_TSS.pdf
        ├── DMR_genomic_distribution.pdf
        └── DMR_volcano.pdf
```

---

## Authors

- Niek Wit
  - University of Cambridge
  - [ORCID profile](https://orcid.org/0009-0002-4330-5333)

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.
