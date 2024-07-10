#! /bin/bash
```
---------------------------------------------------------

Script name: bqsr.sh

Purpose of the script: Apply GATK BaseRecalibrator on a list of files (.bam)

Author: Beatriz Vieira

Date Created: 2024-??-??

Copyright (c) Beatriz Vieira, 2024
Email: marybtv@itqb.unl.pt

Whole-genome polymorphisms and phylogeny of rice varieties circulating in the European market
Maria Beatriz, Hugo Rodrigues, Pedro Barros, Margarida Oliveira

---------------------------------------------------------

Inputs: 1) align_{variety name} directories in the same directory of script run
        2) Genome file (.fa)

Outputs: 1) Recalibrated files (.bam)

---------------------------------------------------------

Required Tools

gatk (v4.2.6.1)

---------------------------------------------------------
```

# Path to reference genome
ref="./Oryza_sativa.IRGSP-1.0.dna.toplevel.fa"
# Path to genome variants file
varlist="../genomeVariants/oryza_sativa.vcf"

# -------------------------------------------------------

while read f
do

variety="$(basename ${f%_*})"
echo $variety

bqfile="align_${variety}/${variety}_RGdedup_recal_data.table"
  gatk BaseRecalibrator \
  -I ${f} \
  -R $ref \
  --known-sites $varlist \
  -O $bqfile \

output="align_${variety}/${variety}_recalibrated.bam"
  gatk ApplyBQSR \
  -I ${f} \
  -R $ref \
  --bqsr-recal-file $bqfile \
  -O ${output} \

done <$1
