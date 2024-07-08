#! /bin/bash
```
---------------------------------------------------------

Script name: Variant_Enrichment.R

Purpose of the script: Map paired-end fastq files to reference genome.

Author: Beatriz Vieira

Date Created: 2024-??-??

Copyright (c) Beatriz Vieira, 2024
Email: marybtv@itqb.unl.pt

Whole-genome polymorphisms and phylogeny of rice varieties circulating in the European market
Maria Beatriz, Hugo Rodrigues, Pedro Barros, Margarida Oliveira

---------------------------------------------------------

Inputs: 1) Fastq files in the same directory of script run
        2) Genome file (.fa) in the same directory of script run

Outputs: 1) Mapped varieties (.sam)

---------------------------------------------------------

Required Tools

bwa-mem (v0.7.17)
samtools 

---------------------------------------------------------
```

# Map paired-end fastq files to reference genome
for i in $(ls -1 *_1.fastq.gz) ;

do

variety=${i%_1.fastq.gz}

echo "aligning ${variety}"

bwa mem -t 10 ../../genome/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa  ${variety}_1.fastq.gz ${variety}_2.fastq.gz > ${variety}.sam

done

# Convert SAM to BAM files
while read f

do

variety="$(basename ${f%.*})"

inp=${variety}.sam
outp=${variety}.bam

echo "converting ${inp} to ${outp}"

cmd="samtools view -h -S -b -o  align_${variety}/${outp} align_${variety}/${inp}"
$cmd

done <$1


