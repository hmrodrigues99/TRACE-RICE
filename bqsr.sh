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

gatk

---------------------------------------------------------
```

# Path to reference genome
ref="./Oryza_sativa.IRGSP-1.0.dna.toplevel.fa"
# Path to genome variants file
varlist="../genomeVariants/oryza_sativa.vcf"
# Path to GATK
gatk=""/home/mariabtv/bin/gatk-4.2.6.1/gatk-package-4.2.6.1-spark.jar""

# -------------------------------------------------------

while read f
do

variety="$(basename ${f%_*})"
echo $variety
outfile=align_${variety}/${variety}_dedup_recal_data.table
bqfile=align_${variety}/${variety}_dedup_recal_data_$i.table
output=align_${variety}/${variety}_dedup_recal_$i.bam

 gatk BaseRecalibratorSpark \
-R $ref \
-I align_${variety}/${variety}_dedup.bam \
-O $outfile \
-- --spark-runner LOCAL --spark-master local[40] \
--conf $gatk=$workPath
 gatk --java-options ▒-Xmx8G▒ ApplyBQSRSpark \
-R $ref \
-I align_${variety}/${variety}_dedup.bam \
-bqsr $bqfile \
--static-quantized-quals 10 --static-quantized-quals 20 \
--static-quantized-quals 30 -O $output \
-- --spark-runner LOCAL --spark-master local[40] \
--conf $gatk=$workPath

outfile="align_${variety}/${variety}_RGdedup_recal_data.table"
  gatk BaseRecalibrator \
  -I align_${variety}/${variety}_dedup_recal_$i.bam \
  -R $ref \
  --known-sites $varlist \
  -O $outfile \


bqfile="align_${variety}/${variety}_RGdedup_recal_data.table"
output="align_${variety}/${variety}_recalibrated.bam"

  gatk ApplyBQSR \
  -R $ref \
  -I align_${variety}/${variety}_RGdedup_recal_data.table \
  --bqsr-recal-file $bqfile \
  -O ${output} \

done <$1
