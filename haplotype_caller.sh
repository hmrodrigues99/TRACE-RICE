#! /bin/bash
# -------------------------------------------------------
#
# Script name: haplotype_caller.sh
#
# Purpose of the script: Apply GATK HaplotypeCaller on a list of recalibrated files (.bam)
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
# Outputs: 1) Genotype variant calling format files (.gvcf)
#
# Run: ./haplotype_caller.sh sample_list.txt
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

while read f
do

variety="$(basename ${f%_*})"
echo "Starting Haplotype Caller for $variety variety"

inrbam="align_${variety}/${variety}_recalibrated.bam"
outgvcf="align_${variety}/${variety}_g.vcf.gz"

gatk HaplotypeCaller \
  --native-pair-hmm-threads 20 \
  -R $ref \
  -I $inrbam \
  -O $outgvcf \
  -ERC GVCF

echo "Haplotype Caller for ${variety} variety is concluded"

done<$1
