#! /bin/bash
# -------------------------------------------------------
#
# Script name: map_bwamem.sh
#
# Purpose of the script: Map paired-end fastq files to a reference genome.
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
# Inputs: 1) Fastq files in the same directory of script run
#        2) Genome file (.fa)
#
# Outputs: 1) Mapped varieties (.sam)
#
# -------------------------------------------------------
#
# Required Tools
#
#bwa-mem (v0.7.17)
#
# -------------------------------------------------------


# Path to reference genome
ref="./Oryza_sativa.IRGSP-1.0.dna.toplevel.fa"

# -------------------------------------------------------

# Map paired-end fastq files to reference genome
for i in $(ls -1 *_1.fastq.gz) ;

do

# Create a specific directory for each sample/variety
variety=${i%_1.fastq.gz}
mkdir align_${variety}
# And a file listing all samples/directories created
echo "${variety}" >> sample_list.txt

echo "aligning ${variety}"

bwa mem -t 10 $ref ${variety}_1.fastq.gz ${variety}_2.fastq.gz > align_${variety}/${variety}.sam

done     
