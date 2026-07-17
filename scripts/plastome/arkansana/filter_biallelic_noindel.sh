#!/usr/bin/env bash
# filter_biallelic_noindel.sh
# Usage: ./filter_biallelic_noindel.sh input.vcf[.gz] output.vcf.gz
# Steps:
#  - Drop sample "Physaria_517"
#  - Mask all hets to missing
#  - Remove sites with any missing genotypes
#  - Keep only variable biallelic SNPs
#  - Strip all INFO, keep only GT
#  - Validate that all remaining sites are variable

set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <input.vcf[.gz]> <output.vcf.gz>"
  exit 1
fi

command -v bcftools >/dev/null || { echo "Error: bcftools not found in PATH"; exit 127; }

INPUT="$1"
OUTPUT="$2"
# DROP_SAMPLE="Physaria_517"

[ -f "$INPUT" ] || { echo "Error: input not found: $INPUT"; exit 2; }

# Ensure +setGT plugin is available
if ! bcftools +setGT -h >/dev/null 2>&1; then
  echo "Error: bcftools +setGT plugin not available (check BCFTOOLS_PLUGINS)."
  exit 127
fi

TMPDIR="$(mktemp -d 2>/dev/null || mktemp -dt filtervcf)"
trap 'rm -rf "$TMPDIR"' EXIT

STEP_IN="$INPUT"
STEP_EXCL="$TMPDIR/step_exclude.vcf.gz"
STEP0="$TMPDIR/step0.nohet.vcf.gz"       # hets masked to missing
STEP0NM="$TMPDIR/step0.nomissing.vcf.gz" # sites with any missing removed
STEP1="$TMPDIR/step1.biallelic_snps.vcf.gz"
CHECK_FILE="$TMPDIR/check_monomorphic.txt"

# # (-1) Drop the specified sample if present
# if bcftools query -l "$INPUT" | grep -qx "$DROP_SAMPLE"; then
#   bcftools view -s "^${DROP_SAMPLE}" -Oz -o "$STEP_EXCL" "$INPUT"
#   STEP_IN="$STEP_EXCL"
# fi

# 0) Mask all heterozygous genotypes to missing
bcftools +setGT -Oz -o "$STEP0" "$STEP_IN" -- -t q -n . -i 'GT="het"'

# 0b) Remove sites with ANY missing genotypes
bcftools view -e 'COUNT(GT="mis")>0' -Oz -o "$STEP0NM" "$STEP0"

# 1) Keep only variable biallelic SNPs (must have both ref and alt homozygotes)
#    Using COUNT on GT avoids depending on INFO/AC.
bcftools view -m2 -M2 -v snps -i 'COUNT(GT="alt")>0 && COUNT(GT="ref")>0' -Oz -o "$STEP1" "$STEP0NM"

# 2) Drop all INFO; keep only GT in FORMAT
bcftools annotate -x INFO,^FMT/GT -Oz -o "$OUTPUT" "$STEP1"

# 3) Index final VCF
bcftools index -t "$OUTPUT"

# 4) VALIDATION: ensure there are no monomorphic sites (all 0/0 or all 1/1)
#    Write any offenders (variant lines only) to a file; if non-empty, fail.
bcftools view -H -i 'COUNT(GT="alt")==0 || COUNT(GT="ref")==0' -Ov -o "$CHECK_FILE" "$OUTPUT"

if [ -s "$CHECK_FILE" ]; then
  BAD=$(wc -l < "$CHECK_FILE")
  echo "Validation FAILED: $BAD monomorphic site(s) found after filtering."
  echo "Example (first offending line):"
  head -n 1 "$CHECK_FILE"
  exit 3
else
  echo "Validation passed: all sites are variable (contain both ref and alt homozygotes)."
fi

echo "Done: $OUTPUT"
