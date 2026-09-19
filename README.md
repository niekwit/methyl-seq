# Snakemake workflow: `methyl-seq`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.25.5-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/niekwit/methyl-seq/actions/workflows/main.yaml/badge.svg?branch=main)](https://github.com/niekwit/methyl-seq/actions/workflows/main.yaml?query=branch%3Amain)
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

Most final outputs (bigwig tracks, boxplots, the ICR heatmap) are per-`condition`, not per-`sample`: replicate samples sharing a condition are averaged/merged together (e.g. `results/bigwig/{condition}.bw`). Per-sample outputs (alignment, deduplication, coverage) still exist under `results/bismark/{sample}/` for inspection, but condition is the unit most downstream analyses report at.

### 4. Configure `config/config.yaml`

```yaml
genome: hg38 # hg19 | hg38 | mm38 | mm39 | dm6
ensembl_genome_build: 114 # Ensembl release number
temp_dir: /tmp # use "/local/" on Cambridge HPC

# Trim Galore clipping arguments (adjust for EM-seq or WGBS library prep)
trim_galore_args: "--clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10"

bismark:
  align: "" # extra arguments for bismark alignment
  deduplicate: "" # extra arguments for deduplicate_bismark
  extract: "" # extra arguments for bismark_methylation_extractor
  coverage: "" # extra arguments for the same command's --cytosine_report step

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

`temp_dir` is where Trim Galore and other tools write scratch files — point it at fast local/scratch storage on a cluster rather than a shared home directory. `trim_galore_args` should be adjusted for your library prep: the default hard-clips 10 bp from both ends of both mates, appropriate for EM-seq's end-repair bias; WGBS libraries may need different values. Each `bismark:` value is appended as extra command-line arguments to the corresponding step (`align` → `bismark`, `deduplicate` → `deduplicate_bismark`, `extract`/`coverage` → `bismark_methylation_extractor`, which handles both extraction and coverage/cytosine-report generation in one invocation — leave a value empty to use Bismark's own defaults).

#### QC: PCA and sample clustering

`deeptools:bigwig_summary` and `deeptools:plotPCA` control `multiBigwigSummary` and `plotPCA` (both from [deepTools](https://deeptools.readthedocs.io/)), run on the per-condition methylation bigwigs to check how similar/distinct your conditions are. `bigwig_summary:binSize` sets the genomic bin size for summarizing signal (smaller = finer-grained but slower); `extra` on either passes additional arguments straight through to the respective deepTools command. Produces `results/deeptools/PCA.tab` (the underlying values) and `results/plots/PCA.pdf`/`scree.pdf`.

#### DMR analysis

Optional (`DMR:run: True`); uses [methylKit](https://bioconductor.org/packages/methylKit/) to tile the genome into `tile_size`-bp windows (stepping by `step_size`, so overlapping tiles are possible if `step_size < tile_size`), keep only tiles with at least `min_per_group` samples covered in both groups, and test each tile for a methylation difference between `reference_condition`'s samples and every other sample pooled together as the comparison group (a single binary comparison, not one comparison per non-reference condition — relevant if `samples.csv` has more than two conditions). A tile is called a significant DMR if its absolute methylation difference exceeds `difference_threshold` (percentage points) and its q-value is below `qvalue_threshold`. Produces `hypermethylated_DMRs.bed`/`hypomethylated_DMRs.bed` (and `..._annotated.tab`, with nearest-gene/genomic-feature annotation via ChIPseeker) under `results/dmrs/`, plus `DMR_volcano.pdf`, `DMR_genomic_distribution.pdf`, and `DMR_distance_to_TSS.pdf` under `results/plots/dmrs/`.

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

#### Paternally imprinted region (ICR) heatmap

For `mm39` and `hg38`, a heatmap of average %CpG methylation across paternally imprinted control regions (ICRs) is generated automatically — no config needed. It's produced from `results/bigwig/{condition}.bw` via `bigWigAverageOverBed`, so it needs no extra alignment or extraction step beyond what the rest of the workflow already does.

- `mm39`: 4 hand-picked ICRs (`Gpr1-Zdbf2`, `Gtl2/Dlk1`, `H19/Igf2`, `Rasgrf1`) in `workflow/resources/icr_regions_mm39.bed`.
- `hg38`: the 25 canonical human ICRs first enumerated by [Skaar et al. 2012](https://pubmed.ncbi.nlm.nih.gov/23744971/), with hg38 coordinates and gene-symbol annotation from the [humanicr.org](https://humanicr.org/) database ([Sanchez-Delgado et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35786392/)), in `workflow/resources/icr_regions_hg38.bed` — full per-region citations in the companion `icr_regions_hg38.references.tsv`. Regenerate either file with:
  ```bash
  python workflow/scripts/get_icr_regions_hg38.py \
      --out workflow/resources/icr_regions_hg38.bed \
      --refs-out workflow/resources/icr_regions_hg38.references.tsv
  ```
  (a standalone utility, not part of the Snakemake DAG, since it queries humanicr.org and UCSC's REST API live).

Every other genome (`hg19`, `mm38`, `dm6`, `test`) has no ICR heatmap — `test` gets a single-region (`Rasgrf1`-only) subset in its own shifted mini-genome coordinates purely so this code path is exercised in CI.

This produces `results/plots/icr_heatmap.pdf` and `icr_heatmap_data.csv`.

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
│   ├── multiqc_report.html          # FastQC summary across all samples
│   └── multiqc_bismark.html         # Bismark alignment/dedup/methylation-extraction/nucleotide-coverage summary across all samples
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
│   ├── boxplots_data.csv
│   ├── te_boxplots.pdf              # TE class/family/subfamily boxplots (if any boxplot TE class block is configured)
│   ├── te_boxplots_data.csv
│   ├── te_5utr_boxplots.pdf         # LINE1 5' UTR vs. remainder boxplots (if any boxplot TE class block configures utr_analysis)
│   ├── te_5utr_boxplots_data.csv
│   ├── icr_heatmap.pdf              # ICR methylation heatmap (mm39/hg38 only)
│   └── icr_heatmap_data.csv
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

## Testing

CI (`.github/workflows/main.yaml`) runs on every push/PR to `main`: formatting (`snakefmt`, via super-linter), `snakemake --lint`, and a full pipeline run against the small, checked-in `.test/` fixture dataset (a real, subsetted mm39 EM-seq dataset with all annotation tracks carved to match — see `.test/make_test_data.py`). To reproduce any of these locally:

```bash
# Lint the workflow
snakemake --directory .test --snakefile workflow/Snakefile --lint

# Run the full pipeline against the test fixtures
snakemake --directory .test --snakefile workflow/Snakefile --use-conda -c <threads>
```

## Authors

- Niek Wit
  - University of Cambridge
  - [ORCID profile](https://orcid.org/0009-0002-4330-5333)

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.
