#!/bin/bash

SPECIES="Quercus"
VCF_DIR="./"
OUT_DIR="./"
mkdir -p "$OUT_DIR"

set -e

for file in ${VCF_DIR}/*.vcf.gz; do
    base=$(basename "$file" .vcf.gz)
    echo "Processing $base..."

    # SNP filtering
    bcftools view -v snps "$file" | \
    bcftools filter -e 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0' -Oz -o "${OUT_DIR}/${base}.snps.filtered.vcf.gz"
	bcftools index "${OUT_DIR}/${base}.snps.filtered.vcf.gz"

    # INDEL filtering
    bcftools view -v indels "$file" | \
    bcftools filter -e 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' -Oz -o "${OUT_DIR}/${base}.indels.filtered.vcf.gz"
	bcftools index "${OUT_DIR}/${base}.indels.filtered.vcf.gz"

    # Merge filtered SNPs and INDELs
    bcftools concat -a -Oz -o "${OUT_DIR}/${base}.filtered.vcf.gz" \
        "${OUT_DIR}/${base}.snps.filtered.vcf.gz" \
        "${OUT_DIR}/${base}.indels.filtered.vcf.gz"

    # Index final output
    bcftools index "${OUT_DIR}/${base}.filtered.vcf.gz"

    # Optionally clean up intermediate files
    rm "${OUT_DIR}/${base}.snps.filtered.vcf.gz" "${OUT_DIR}/${base}.snps.filtered.vcf.gz.csi"
    rm "${OUT_DIR}/${base}.indels.filtered.vcf.gz" "${OUT_DIR}/${base}.indels.filtered.vcf.gz.csi"

    echo "Finished filtering $base"
done

echo "All files processed."
