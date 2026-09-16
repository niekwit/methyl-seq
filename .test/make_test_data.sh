#!/usr/bin/env bash
#
# make_test_data.sh
# -----------------------------------------------------------------------------
# Build a tiny methyl-seq test dataset for CI from real Bismark output.
#
# It keeps every read (both mates) that Bismark aligned to a single small
# locus (+ flanks), plus a capped number of spike-in control reads, and writes
# them back out as FASTQ under .test/reads/. Optionally it also carves the
# matching mini genome FASTA and mini GTF out of the full reference.
#
# Nothing is subsampled by default: coverage depth over the locus is preserved,
# only the genomic breadth is reduced, so methylKit/boxplots stay meaningful.
#
# Requirements: samtools >= 1.12 (needs `view -N` / --qname-file).
#
# Usage:
#   .test/make_test_data.sh /path/to/results/bismark
#   BISMARK_DIR=/path/to/results/bismark .test/make_test_data.sh
#   REF_FASTA=/path/to/combined_genome.fa .test/make_test_data.sh /path/to/results/bismark
#   REF_FASTA=... REF_GTF=/path/to/annotation.gtf.gz .test/make_test_data.sh /path/to/results/bismark
#   FLANK=250000 THREADS=8 KEEP_COORDS=1 .test/make_test_data.sh /path/to/results/bismark
# -----------------------------------------------------------------------------
set -euo pipefail

case "${1:-}" in -h|--help) sed -n '2,38p' "$0"; exit 0;; esac

# ============================== CONFIG =======================================
# Directory holding <sample>/<sample>.deduplicated.bam
# (positional $1 takes precedence over the BISMARK_DIR env var)
BISMARK_DIR="${1:-${BISMARK_DIR:-}}"
[ -n "$BISMARK_DIR" ] || { echo "ERROR: BISMARK_DIR not set. Pass it as \$1 or the BISMARK_DIR env var." >&2; exit 1; }

# Where to write test data (defaults to the folder this script lives in = .test/)
OUT_DIR="${OUT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
READS_DIR="$OUT_DIR/reads"
RES_DIR="$OUT_DIR/resources"

# Target region (hg38). The "chr" prefix is optional and auto-matched to the BAM.
REGION_CHR="9"
REGION_START=89761645
REGION_END=89761980
REGION_NAME="Rasgrf1"

# bp added on each side of the target -> more reads + surrounding gene context
FLANK="${FLANK:-100000}"

THREADS="${THREADS:-4}"

# Spike-in / conversion-rate controls: keep some reads so the
# methylation_conversion_rate rule has input. Set INCLUDE_CONTROLS=0 to skip.
INCLUDE_CONTROLS="${INCLUDE_CONTROLS:-1}"
CONTROL_MAX="${CONTROL_MAX:-2000}"   # cap kept control reads per sample

# Optional extra downsampling of the *region* reads, passed to `samtools view -s`
# as seed.fraction (e.g. 42.25 keeps 25%). Empty = keep all. Mates stay paired.
SUBSAMPLE="${SUBSAMPLE:-}"

# Full reference to carve the mini genome from (e.g. resources/combined_genome.fa
# or the Ensembl primary assembly). Leave empty to skip mini-genome creation.
REF_FASTA="${REF_FASTA:-}"

# Optional GTF (plain or .gz) to carve a matching mini annotation from. Its
# coordinates get the same shift/rename as the mini genome. No workflow rule
# consumes the GTF today, but it is handy for regenerating annotation BEDs.
REF_GTF="${REF_GTF:-}"

# 1 = N-pad the 5' side so the mini genome keeps real hg38 coordinates (bigger
#     file, but TxDb annotation / real BED files line up).
# 0 = compact contig, region shifted to start near position 1 (fast; DMR
#     annotation still runs, just against the chromosome start).
KEEP_COORDS="${KEEP_COORDS:-0}"
# ============================================================================

command -v samtools >/dev/null || { echo "ERROR: samtools not found in PATH" >&2; exit 1; }

MINI_CONTIG="${REGION_CHR#chr}"            # header used in the test workflow (Ensembl style)
WIN_START=$(( REGION_START - FLANK )); (( WIN_START < 1 )) && WIN_START=1
WIN_END=$(( REGION_END + FLANK ))
# shift applied to every coordinate written into the test data
if [ "$KEEP_COORDS" = "1" ]; then OFFSET=0; else OFFSET=$(( WIN_START - 1 )); fi

mkdir -p "$READS_DIR" "$RES_DIR"
rm -f "$READS_DIR"/*.fastq.gz          # start clean so re-runs don't mix old output

mapfile -t BAMS < <(ls "$BISMARK_DIR"/*/*.deduplicated.bam 2>/dev/null || true)
[ "${#BAMS[@]}" -gt 0 ] || { echo "ERROR: no *.deduplicated.bam under $BISMARK_DIR" >&2; exit 1; }

# --- match the target contig name to whatever the BAM uses -------------------
contig_in() {  # contig_in <name> <header-text>
  awk -v c="SN:$1" '$1=="@SQ"{for(i=1;i<=NF;i++) if($i==c) ok=1} END{exit !ok}' <<<"$2"
}
HDR=$(samtools view -H "${BAMS[0]}")
if   contig_in "$REGION_CHR"        "$HDR"; then CONTIG="$REGION_CHR"
elif contig_in "${REGION_CHR#chr}"  "$HDR"; then CONTIG="${REGION_CHR#chr}"
elif contig_in "chr${REGION_CHR#chr}" "$HDR"; then CONTIG="chr${REGION_CHR#chr}"
else echo "ERROR: neither '$REGION_CHR' nor '${REGION_CHR#chr}' is a contig in ${BAMS[0]}" >&2; exit 1
fi

echo "Region      : ${CONTIG}:${REGION_START}-${REGION_END} (${REGION_NAME})"
echo "Window +flank: ${CONTIG}:${WIN_START}-${WIN_END}  (${FLANK} bp each side)"
echo "Output       : $READS_DIR"
echo

# ============================ per-sample =====================================
for bam in "${BAMS[@]}"; do
  sample=$(basename "$(dirname "$bam")")
  echo ">> $sample"
  names=$(mktemp); rnames=$(mktemp); cnames=$(mktemp)

  # Pass 1: single scan -> QNAMEs overlapping the flanked window (rnames) and
  # QNAMEs on the spike-in control contigs (cnames). Secondary/supp dropped.
  samtools view -@ "$THREADS" -F 0x900 "$bam" \
    | awk -v c="$CONTIG" -v s="$WIN_START" -v e="$WIN_END" -v rf="$rnames" -v cf="$cnames" '
        $3==c && $4<=e && ($4+length($10))>=s { print $1 > rf; next }
        $3 ~ /^(phage_lambda|plasmid_puc19c|Lambda|pUC19|J02459\.1)$/ { print $1 > cf }'
  sort -u "$rnames" -o "$rnames"
  sort -u "$cnames" -o "$cnames"
  region_n=$(wc -l < "$rnames")

  # cap kept control reads (head on a file -> no SIGPIPE)
  ctrl_n=0
  cp "$rnames" "$names"
  if [ "$INCLUDE_CONTROLS" = "1" ]; then
    head -n "$CONTROL_MAX" "$cnames" >> "$names"
    ctrl_n=$(head -n "$CONTROL_MAX" "$cnames" | wc -l)
    sort -u "$names" -o "$names"
  fi
  rm -f "$rnames" "$cnames"
  echo "   region read names: $region_n   control read names: $ctrl_n"
  [ "$region_n" -gt 0 ] || echo "   WARNING: no reads found in region for $sample"

  # paired- or single-end? (read one FLAG; bit 0x1 = paired. Avoid a pipe whose
  # exit status pipefail would poison via SIGPIPE.)
  first_flag=$(samtools view -F 0x900 "$bam" 2>/dev/null | awk 'NR==1{print $2; exit}') || true
  if [ $(( ${first_flag:-0} % 2 )) -eq 1 ]; then PE=1; else PE=0; fi

  sflag=(); [ -n "$SUBSAMPLE" ] && sflag=(-s "$SUBSAMPLE")

  # Pass 2: pull those reads (both mates via QNAME) and convert to FASTQ.
  # Bismark stores reads in aligned orientation; `samtools fastq` restores the
  # original (bisulfite-converted) read from the FLAG, so re-alignment is valid.
  if [ "$PE" = "1" ]; then
    samtools view -@ "$THREADS" -b -F 0x900 -N "$names" "${sflag[@]}" "$bam" \
      | samtools collate -@ "$THREADS" -u -O - \
      | samtools fastq -@ "$THREADS" -n \
          -1 "$READS_DIR/${sample}_R1_001.fastq.gz" \
          -2 "$READS_DIR/${sample}_R2_001.fastq.gz" \
          -0 /dev/null -s /dev/null -
  else
    samtools view -@ "$THREADS" -b -F 0x900 -N "$names" "${sflag[@]}" "$bam" \
      | samtools fastq -@ "$THREADS" -n -0 "$READS_DIR/${sample}.fastq.gz" -
  fi
  rm -f "$names"
done

# ============================ mini genome ====================================
echo
if [ -n "$REF_FASTA" ] && [ -f "$REF_FASTA" ]; then
  echo ">> mini genome from $REF_FASTA"
  [ -f "$REF_FASTA.fai" ] || samtools faidx "$REF_FASTA"
  if   cut -f1 "$REF_FASTA.fai" | grep -qx "$REGION_CHR";       then RCONTIG="$REGION_CHR"
  elif cut -f1 "$REF_FASTA.fai" | grep -qx "${REGION_CHR#chr}"; then RCONTIG="${REGION_CHR#chr}"
  else echo "   ERROR: '$REGION_CHR' not found in $REF_FASTA.fai" >&2; exit 1
  fi

  if [ "$KEEP_COORDS" = "1" ]; then
    { head -c "$(( WIN_START - 1 ))" /dev/zero | tr '\0' 'N'
      samtools faidx "$REF_FASTA" "${RCONTIG}:${WIN_START}-${WIN_END}" | tail -n +2 | tr -d '\n'
    } | fold -w 60 > "$RES_DIR/.seq"
  else
    samtools faidx "$REF_FASTA" "${RCONTIG}:${WIN_START}-${WIN_END}" \
      | tail -n +2 | tr -d '\n' | fold -w 60 > "$RES_DIR/.seq"
  fi
  { printf '>%s\n' "$MINI_CONTIG"; cat "$RES_DIR/.seq"; echo; } > "$RES_DIR/genome.fa"
  rm -f "$RES_DIR/.seq"
  samtools faidx "$RES_DIR/genome.fa"
  cut -f1,2 "$RES_DIR/genome.fa.fai" > "$RES_DIR/chrom.sizes"
  printf '%s\t%d\t%d\t%s\n' "$MINI_CONTIG" "$(( REGION_START - OFFSET - 1 ))" \
         "$(( REGION_END - OFFSET ))" "$REGION_NAME" > "$RES_DIR/target_region.bed"
  echo "   $RES_DIR/genome.fa  (contig '$MINI_CONTIG', coordinate offset applied to test data = $OFFSET)"
  echo "   $RES_DIR/target_region.bed  (target locus in mini-genome coordinates)"
else
  echo ">> REF_FASTA not set/found -> skipped mini genome."
  echo "   Re-run with e.g.  REF_FASTA=resources/combined_genome.fa .test/make_test_data.sh"
fi

# ============================ mini GTF =======================================
echo
gtf_cat() { case "$1" in *.gz) zcat "$1";; *) cat "$1";; esac; }

if [ -n "$REF_GTF" ] && [ -f "$REF_GTF" ]; then
  echo ">> mini GTF from $REF_GTF"

  # match the GTF's seqname style (Ensembl '9' vs UCSC 'chr9'). Capture awk's
  # output, not the pipeline exit status (awk's early `exit` SIGPIPEs gtf_cat,
  # which pipefail would otherwise report as failure).
  GCONTIG=""
  for cand in "${REGION_CHR#chr}" "chr${REGION_CHR#chr}" "$REGION_CHR"; do
    hit=$(gtf_cat "$REF_GTF" | awk -F'\t' -v c="$cand" \
            '!/^#/ && $1==c { print "y"; exit }') || true
    [ "$hit" = y ] && { GCONTIG="$cand"; break; }
  done
  [ -n "$GCONTIG" ] || { echo "   ERROR: no contig for '$REGION_CHR' in $REF_GTF" >&2; exit 1; }

  # keep features overlapping the window, clamp overhangs, shift + rename
  gtf_cat "$REF_GTF" \
    | awk -F'\t' -v OFS='\t' -v c="$GCONTIG" -v ws="$WIN_START" -v we="$WIN_END" \
          -v off="$OFFSET" -v mc="$MINI_CONTIG" '
        /^#/ { next }
        $1==c && $4<=we && $5>=ws {
          lo = ws - off; hi = we - off
          s = $4 - off; if (s < lo) s = lo
          e = $5 - off; if (e > hi) e = hi
          $1 = mc; $4 = s; $5 = e
          print
        }' > "$RES_DIR/genes.gtf"

  gzip -c "$RES_DIR/genes.gtf" > "$RES_DIR/genes.gtf.gz"
  n=$(grep -vc '^#' "$RES_DIR/genes.gtf" || true)
  echo "   $RES_DIR/genes.gtf(.gz)  ($n features; gene models overhanging the window are clamped)"
elif [ -n "$REF_GTF" ]; then
  echo ">> REF_GTF '$REF_GTF' not found -> skipped mini GTF."
fi

echo
echo "Done."
du -sh "$READS_DIR"/* 2>/dev/null || true
