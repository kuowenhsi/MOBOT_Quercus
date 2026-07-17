#!/usr/bin/env bash
# filter_biallelic_noindel.sh
#
# Usage:
#   bash filter_biallelic_noindel.sh input.vcf.gz output.vcf.gz
#
# Steps:
#   1. Remove specified samples when present
#   2. Convert heterozygous genotypes to missing
#   3. Remove sites containing any missing genotype
#   4. Retain variable biallelic SNPs
#   5. Remove INFO and retain only FORMAT/GT
#   6. Validate the final VCF

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input.vcf[.gz]> <output.vcf.gz>" >&2
    exit 1
fi

INPUT="$1"
OUTPUT="$2"

DROP_SAMPLES=(
    "Quercus_414"
    "Quercus_415"
    "Quercus_419"
    "Quercus_422"
    "Quercus_424"
    "Quercus_427"
    "Quercus_433"
    "Quercus_435"
)

command -v bcftools >/dev/null 2>&1 || {
    echo "Error: bcftools not found in PATH." >&2
    exit 127
}

[[ -f "$INPUT" ]] || {
    echo "Error: input file not found: $INPUT" >&2
    exit 2
}

if [[ "$OUTPUT" != *.vcf.gz ]]; then
    echo "Error: output filename must end in .vcf.gz" >&2
    exit 2
fi

if ! bcftools +setGT -h >/dev/null 2>&1; then
    echo "Error: bcftools +setGT plugin is unavailable." >&2
    echo "Check the BCFTOOLS_PLUGINS environment variable." >&2
    exit 127
fi

TMPDIR="$(mktemp -d 2>/dev/null || mktemp -dt filtervcf)"
trap 'rm -rf "$TMPDIR"' EXIT

SAMPLE_LIST="$TMPDIR/input.samples.txt"
STEP_EXCL="$TMPDIR/step_exclude.vcf.gz"
STEP_NOHET="$TMPDIR/step_nohet.vcf.gz"
STEP_NOMISSING="$TMPDIR/step_nomissing.vcf.gz"
STEP_SNPS="$TMPDIR/step_biallelic_snps.vcf.gz"
CHECK_FILE="$TMPDIR/check_monomorphic.vcf"

bcftools query -l "$INPUT" > "$SAMPLE_LIST"

########################################################################
# 1. Determine which requested samples are actually present
########################################################################

PRESENT_DROP_SAMPLES=()

for SAMPLE in "${DROP_SAMPLES[@]}"; do
    if grep -Fxq "$SAMPLE" "$SAMPLE_LIST"; then
        PRESENT_DROP_SAMPLES+=("$SAMPLE")
    else
        echo "Warning: sample not found and will be ignored: $SAMPLE" >&2
    fi
done

STEP_IN="$INPUT"

if (( ${#PRESENT_DROP_SAMPLES[@]} > 0 )); then
    DROP_CSV="$(IFS=,; echo "${PRESENT_DROP_SAMPLES[*]}")"

    echo "Removing samples:"
    printf '  %s\n' "${PRESENT_DROP_SAMPLES[@]}"

    bcftools view \
        --samples "^${DROP_CSV}" \
        -Oz \
        -o "$STEP_EXCL" \
        "$INPUT"

    STEP_IN="$STEP_EXCL"
else
    echo "None of the requested samples were found; no samples removed."
fi

########################################################################
# 2. Mask heterozygous genotypes
########################################################################

echo "Masking heterozygous genotypes..."

bcftools +setGT "$STEP_IN" \
    -Oz \
    -o "$STEP_NOHET" \
    -- \
    -t q \
    -n . \
    -i 'GT="het"'

########################################################################
# 3. Remove sites containing any missing genotype
########################################################################

echo "Removing sites with missing genotypes..."

bcftools view \
    -e 'COUNT(GT="mis")>0' \
    -Oz \
    -o "$STEP_NOMISSING" \
    "$STEP_NOHET"

########################################################################
# 4. Retain variable biallelic SNPs
########################################################################

echo "Retaining variable biallelic SNPs..."

bcftools view \
    -m2 \
    -M2 \
    -v snps \
    -i 'COUNT(GT="alt")>0 && COUNT(GT="ref")>0' \
    -Oz \
    -o "$STEP_SNPS" \
    "$STEP_NOMISSING"

########################################################################
# 5. Remove INFO and retain only FORMAT/GT
########################################################################

echo "Removing INFO fields and non-GT FORMAT fields..."

bcftools annotate \
    -x 'INFO,^FORMAT/GT' \
    -Oz \
    -o "$OUTPUT" \
    "$STEP_SNPS"

########################################################################
# 6. Index output
########################################################################

bcftools index -f -t "$OUTPUT"

########################################################################
# 7. Validation
########################################################################

echo "Validating output..."

bcftools view \
    -H \
    -i 'COUNT(GT="alt")==0 || COUNT(GT="ref")==0' \
    -Ov \
    -o "$CHECK_FILE" \
    "$OUTPUT"

if [[ -s "$CHECK_FILE" ]]; then
    BAD="$(wc -l < "$CHECK_FILE" | tr -d ' ')"

    echo "Validation FAILED: $BAD monomorphic site(s) found." >&2
    echo "First offending record:" >&2
    head -n 1 "$CHECK_FILE" >&2
    exit 3
fi

N_INPUT="$(bcftools view -H "$INPUT" | wc -l | tr -d ' ')"
N_OUTPUT="$(bcftools view -H "$OUTPUT" | wc -l | tr -d ' ')"
N_SAMPLES="$(bcftools query -l "$OUTPUT" | wc -l | tr -d ' ')"

echo
echo "Validation passed."
echo "Input records:  $N_INPUT"
echo "Output SNPs:    $N_OUTPUT"
echo "Final samples:  $N_SAMPLES"
echo "Output:         $OUTPUT"
echo "Index:          ${OUTPUT}.tbi"
