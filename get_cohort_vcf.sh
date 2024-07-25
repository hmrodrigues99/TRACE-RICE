#! /bin/bash
# -------------------------------------------------------
#
# Script name: get_cohort_vcf.sh
#
# Purpose of the script: Apply GATK CombineGVCFs and GenotypeGVCFs on a list of .gvcf files
#                        from GATK Haplotype Caller
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
# Inputs:  1) align_{variety name} directories in the same directory of script run
#          2) Genome file (.fa)
#
# Outputs: 1) Cohort variant calling format file (.vcf)
#
# Run: ./get_cohort_vcf.sh
#
# -------------------------------------------------------
#
# Required Tools
#
# gatk (v4.2.6.1)
#
# -------------------------------------------------------

# Path to reference genome
ref="./Oryza_sativa.IRGSP-1.0.dna.toplevel.fa"

# -------------------------------------------------------

# Join all sample .gvcf files into a single all_gvcf.vcf.gz file

gatk   CombineGVCFs \
       -R $ref \
       -V align_Bomba/Bomba_g.vcf.gz \
       -V align_Puntal/Puntal_g.vcf.gz \
       -O all_g.vcf.gz

# Get cohort VCF with all variants for all samples

gatk   GenotypeGVCFs \
       -R $ref \
       -V all_g.vcf.gz \
       -O cohort.vcf.gz
