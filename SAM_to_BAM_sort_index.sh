#! /bin/bash
```
---------------------------------------------------------

Script name: SAM_to_BAM_sort_index.sh

Purpose of the script: Convert SAM to BAM files, sort and index them.

Author: Beatriz Vieira

Date Created: 2024-??-??

Copyright (c) Beatriz Vieira, 2024
Email: marybtv@itqb.unl.pt

Whole-genome polymorphisms and phylogeny of rice varieties circulating in the European market
Maria Beatriz, Hugo Rodrigues, Pedro Barros, Margarida Oliveira

---------------------------------------------------------

Inputs: 1) SAM files in the same directory of script run

Outputs: 1) Mapped, sorted and indexed varieties (.bam)

---------------------------------------------------------

Required Tool(s)

samtools (v1.7)

---------------------------------------------------------
```

# For each .sam file within n align_{variety_name} directories
while read f

do

variety="$(basename ${f%.*})"

inp=${variety}.sam
outp=${variety}.bam

# Convert SAM to BAM files
echo "converting ${inp} to ${outp}"

cmd="samtools view -h -S -b -o  align_${variety}/${outp} align_${variety}/${inp}"
echo ${cmd}
$cmd

# Sort BAM files
inbam=${variety}.bam
outsorted=${variety}.sorted

echo "sorting ${inbam} into ${outsorted}"

cmdsort="samtools sort -o align_${variety}/${outsorted} align_${variety}/${inbam}"
echo ${cmdsort}
$cmdsort

# Index BAM files
insorted=${variety}.sorted
indexed=${insorted}

echo "indexing ${insorted} into ${indexed}"

cmdindex="samtools index align_${variety}/${insorted}"
echo ${cmdindex}
$cmdindex

done <$1
