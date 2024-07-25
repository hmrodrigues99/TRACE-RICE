#! /bin/bash
# -------------------------------------------------------
#
# Script name: filter_vcf.sh
#
# Purpose of the script: Apply hard-filters to a VCF file using GATK SelectVariants
#                        It uses optimized filters for SNPs and InDels separately
#
# Author: Beatriz Vieira
#
# Date Created: 2024-??-??
#
# Copyright (c) Beatriz Vieira, 2024
# Email: marybtv@itqb.unl.pt
#
# Whole-genome polymorphisms and phylogeny of rice varieties circulating in the European market
# Maria Beatriz, Hugo Rodrigues, Pedro Barros, Margarida Oliveira
#
# -------------------------------------------------------
#
# Inputs:  1) Cohort VCF file (obtained using GATK GenotypeGVCFs)
#
# Outputs: 1) Filtered cohort VCF file
#
# -------------------------------------------------------
#
# Required Tools
#
# gatk (v4.2.6.1)
#
# ------------------------------------------------------

# Get a SNP-only subset

gatk SelectVariants \
    -V cohort.vcf.gz \
    -select-type SNP \
    -O snps.vcf.gz

# Get a InDel-only subset

gatk SelectVariants \
    -V cohort.vcf.gz \
    -select-type INDEL \
    -O indels.vcf.gz

# Hard-filter SNPs

gatk VariantFiltration \
    -V snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O snps_filtered.vcf.gz

# Hard-filter InDels

gatk VariantFiltration \
    -V indels.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O indels_filtered.vcf.gz

# Merge SNP and InDels back together

gatk MergeVcfs \
    -I snps_filtered.vcf.gz \
    -I indels_filtered.vcf.gz \
    -O filtered_cohort.vcf.gz \

# Get a separate file with variants that passed all filters

gatk SelectVariants \
    -V filtered_cohort.vcf.gz \
    -O PASS_cohort.vcf.gz \
    --exclude-filtered
