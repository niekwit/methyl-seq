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
  # Standard regions (whole_genome, genic, exon, intron, intergenic, and
  # promoter/cpg_islands where available for the genome) are generated
  # automatically. Optionally add extra custom regions on top:
  # regions:
  #   my_custom_region: "config/annotations/my_custom_region.bed"
  #
  # Optionally add transposable-element (TE) class blocks, keyed by a
  # RepeatMasker repClass value, for a separate combined te_boxplots.pdf.
  # A class's mere presence enables it. Repeat elements overlapping gene
  # bodies are excluded first, then filtered by min_length. `family`
  # (optional, comma-separated repFamily values) and `subfamilies`
  # (optional, repName prefixes, requires `family`) add further
  # breakdowns plotted alongside the class total.
  LINE:
    min_length: 6000
    family: L1
    subfamilies:
      - L1MdA
      - L1MdF
    # Optional: 5' UTR vs. remainder methylation boxplot for this family
    # (requires family to be exactly "L1")
    # utr_analysis:
    #   genbank: config/annotations/L1_consensus.gb
  LTR:
    min_length: 6000

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

#### Boxplot regions

`boxplots.pdf` always includes the standard regions generated automatically for the configured genome (`whole_genome`, `genic`, `exon`, `intron`, `intergenic`, and `promoter`/`cpg_islands` where available) — no config needed for these.

To plot additional custom regions alongside them, add a `regions:` map under `boxplot` (name → BED file path, one boxplot facet per entry):
```yaml
boxplot:
  ...
  regions:
    my_custom_region: "config/annotations/my_custom_region.bed"
```
The BED file just needs `chrom`, `start`, `end` columns. A region name can't reuse a standard region name (`genic`, `exon`, ...) or a TE class/family/subfamily name (below) — the workflow will error out at startup if it does.

#### TE (transposable element) boxplots

Optionally add one or more TE class blocks directly under `boxplot`, keyed by a RepeatMasker `repClass` value (e.g. `LINE`, `LTR`, `SINE`, `DNA`) — a class block's mere presence turns it on. These produce a separate, combined `te_boxplots.pdf` (not mixed into `boxplots.pdf`), one facet per class/family/subfamily, in the order they're configured.

```yaml
boxplot:
  ...
  LINE:
    min_length: 6000   # required: minimum element length (bp)
    family: L1          # optional: comma-separated repFamily value(s), e.g. "L1,L2"
    subfamilies:          # optional: repName prefixes, requires `family`
      - L1MdA
      - L1MdF
  LTR:
    min_length: 6000       # class-only: no family/subfamily breakdown
```
For each class block, repeat elements overlapping a gene body are excluded first, then the remainder is filtered by `min_length`. This gives, in order: the class total (all filtered elements), then one facet per listed `family` (`repFamily` match), then one facet per listed `subfamilies` entry (`repName` prefix match, e.g. `"L1MdA"` matches `"L1MdA"`/`"L1MdA_I"`/`"L1MdA_II"`/... but not `"L1MdAxyz"`, and is restricted to the configured family/families). A class, family, or subfamily name that matches zero elements is a hard error (usually a typo or an overly strict `min_length`) rather than a silently empty facet.

Enabling any TE class block also makes the `genic`/`promoter`/`cpg_islands` regions in `boxplots.pdf` transposon-subtracted, so their CpG signal isn't contaminated by repeat-element methylation patterns.

TE analysis requires a genome with a RepeatMasker track (all supported genomes except `dm6`).

#### LINE1 5' UTR vs. remainder boxplots

Optionally add a `utr_analysis:` sub-block under a TE class block to split each near-full-length LINE1 (L1) element into its 5' UTR and "remainder" (ORF1 + linker + ORF2 + 3' UTR), and boxplot %CpG methylation of the two separately. This is only valid when that class block's `family` is exactly `L1` (LINE1 is the only TE with this ORF1/ORF2 UTR architecture) — the workflow errors out at startup otherwise.

```yaml
boxplot:
  ...
  LINE:
    min_length: 6000
    family: L1
    subfamilies:
      - L1MdA
      - L1MdF
    utr_analysis:
      genbank: config/annotations/L1_consensus.gb   # required
      min_identity: 0.6    # optional, default 0.6
      min_coverage: 0.5    # optional, default 0.5
      min_mapq: 20          # optional, default 20
      single_orf: false     # optional, default false
```

`genbank` must be a GenBank record for a consensus L1 with `CDS` features whose `/product` qualifiers are exactly `ORF1` and `ORF2` (e.g. mouse [M13002](https://www.ncbi.nlm.nih.gov/nuccore/M13002), human [AF148856](https://www.ncbi.nlm.nih.gov/nuccore/AF148856)). For each near-full-length L1 element (the configured `family`'s own filtered set), ORF1/ORF2 are located on that element's own sequence via minimap2; `min_identity`/`min_coverage` control how confident a hit must be, and `min_mapq` its minimum mapping quality. By default an element needs both ORFs confidently located to be annotated; `single_orf: true` also keeps elements where only one was found (the "UTR" on the missing ORF's side then also contains that ORF and the inter-ORF sequence).

This produces a separate `results/plots/te_5utr_boxplots.pdf`, faceted by TE (the family total, then its configured subfamilies) x region (5' UTR / L1 remainder) — each L1 element's UTR/remainder span is its own boxplot data point (not split into the CpG-probe-sized chunks the other boxplots use).

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
│   ├── boxplots.pdf                 # CpG methylation boxplots (if boxplot.plot: True)
│   ├── te_boxplots.pdf              # TE class/family/subfamily boxplots (if any boxplot TE class block is configured)
│   └── te_5utr_boxplots.pdf         # LINE1 5' UTR vs. remainder boxplots (if any boxplot TE class block configures utr_analysis)
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
