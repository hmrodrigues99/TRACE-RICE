#! /bin/bash
```
---------------------------------------------------------

Script name: SAM_to_BAM_sort_index.sh

Purpose of the script: Provices options to convert SAM to BAM format, sorting and indexing.

Author: Beatriz Vieira

Date Created: 2024-??-??

Copyright (c) Beatriz Vieira, 2024
Email: marybtv@itqb.unl.pt

Whole-genome polymorphisms and phylogeny of rice varieties circulating in the European market
Maria Beatriz, Hugo Rodrigues, Pedro Barros, Margarida Oliveira

---------------------------------------------------------

Inputs: 1) align_{variety name} directories in the same directory of script run

Outputs: 1) Mapped, sorted and indexed varieties (.bam)

---------------------------------------------------------

Required Tool(s)

samtools (v1.7)

---------------------------------------------------------
```

# Change to FALSE to disable option(s)
SAM_to_BAM=TRUE
BAM_sort=TRUE
BAM_index=TRUE

# -------------------------------------------------------

while read f

do

# Convert SAM to BAM files
if [ $SAM_to_BAM == TRUE ]
then

inp=${variety}.sam
outp=${variety}.bam

echo "converting ${inp} to ${outp}"

cmd="samtools view -h -S -b -o  align_${variety}/${outp} align_${variety}/${inp}"
echo ${cmd}
$cmd
fi

# Sort BAM files
if [ $BAM_sort == TRUE ]
then

inbam=${variety}.bam
outsorted=${variety}.sorted

echo "sorting ${inbam} into ${outsorted}"

cmdsort="samtools sort -o align_${variety}/${outsorted} align_${variety}/${inbam}"
echo ${cmdsort}
$cmdsort
fi

# Index BAM files
if [ $BAM_index == TRUE ]
then

insorted=${variety}.sorted

echo "indexing ${insorted}"

cmdindex="samtools index align_${variety}/${insorted}"
echo ${cmdindex}
$cmdindex
fi

done <$1
