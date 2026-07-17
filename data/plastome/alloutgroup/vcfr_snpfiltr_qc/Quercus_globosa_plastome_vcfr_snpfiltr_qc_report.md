# Quercus globosa Plastome VCF QC And Filtering Conditions

Generated: 2026-07-14 14:41:08 CDT

## Input

- VCF: Quercus_globosa_plastome.allsite.allsite.exclude8_2397024510.filtered.variant.vcf.gz
- Samples: 209
- Variant records: 442
- SNPfiltR installed: FALSE

SNPfiltR is not installed in this R library, so this run uses vcfR-based equivalent QC summaries and records the SNPfiltR status for reproducibility.

## Filtering Conditions

| step | condition                                                         | expression                                                                                                                      |
| ---- | ----------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------- |
| 1    | Exclude eight samples                                             | Quercus_414, Quercus_415, Quercus_419, Quercus_422, Quercus_424, Quercus_427, Quercus_433, Quercus_435                          |
| 2    | Keep reference sites or biallelic SNP/indel sites                 | N_ALT=0 \|\| (N_ALT=1 && (TYPE="snp" \|\| TYPE="indel"))                                                                        |
| 3    | Remove sites with >=10% missing genotypes before genotype masking | F_MISSING < 0.10                                                                                                                |
| 4    | Mask low-quality variant-site genotypes to missing                | N_ALT>0 & (FMT/DP<10 \| FMT/GQ<30 \| FMT/DP>10000)                                                                              |
| 5    | Mask low-quality reference-site genotypes to missing              | N_ALT=0 & (FMT/DP<10 \| FMT/RGQ<30 \| FMT/DP>10000)                                                                             |
| 6    | Remove sites with >=10% missing genotypes after genotype masking  | F_MISSING < 0.10                                                                                                                |
| 7    | Apply GATK-style site hard filters separately to SNPs and indels  | SNPs: QD>=2, FS<=60, MQ>=40, MQRankSum>=-12.5, ReadPosRankSum>=-8, SOR<=3; indels: QD>=2, FS<=200, ReadPosRankSum>=-20, SOR<=10 |
| 8    | Keep variant SNP/indel sites for this final variant VCF           | N_ALT>0 && (TYPE="snp" \|\| TYPE="indel")                                                                                       |

## QC Summary

| metric                            | value                                                                                |
| --------------------------------- | ------------------------------------------------------------------------------------ |
| VCF                               | Quercus_globosa_plastome.allsite.allsite.exclude8_2397024510.filtered.variant.vcf.gz |
| Generated                         | 2026-07-14 14:41:08 CDT                                                              |
| SNPfiltR installed                | FALSE                                                                                |
| Samples                           | 209                                                                                  |
| Variant records                   | 442                                                                                  |
| SNP records                       | 304                                                                                  |
| Indel records                     | 138                                                                                  |
| PASS records                      | 442                                                                                  |
| Mean site missing rate            | 0.02372                                                                              |
| Max site missing rate             | 0.09569                                                                              |
| Mean sample missing rate          | 0.02372                                                                              |
| Max sample missing rate           | 0.1697                                                                               |
| Mean site heterozygote rate       | 0.1394                                                                               |
| Max site heterozygote rate        | 0.9856                                                                               |
| Mean DP across called genotypes   | 693.6                                                                                |
| Median DP across called genotypes | 606                                                                                  |
| Mean GQ across called genotypes   | 96.42                                                                                |
| Ts/Tv for SNP records             | 1.286                                                                                |

## Output Files

- qc_summary.csv
- filtering_conditions.csv
- site_qc_vcfr.csv
- sample_qc_vcfr.csv
- bcftools_filtering_commands_from_vcf_header.txt
- variant_class_counts.png
- site_missingness_histogram.png
- sample_missingness_heterozygosity.png
- site_quality_depth.png

## Notes

- The filtering conditions are read from the final VCF header and the existing alloutgroup filtering history.
- This report describes the already-filtered variant VCF; it does not rewrite or further filter the VCF.
