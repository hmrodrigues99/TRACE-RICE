#! /bin/bash
# -------------------------------------------------------
#
# Script name: SAM_to_BAM_sort_dups_index.sh
#
# Purpose of the script: Performs if set to TRUE:
#                        1) Convert SAM to BAM format
#                        2) Computes BAM file statistics (samtools stats & flagstats)
#                        3) Add mate score tags
#                        4) Sorting
#                        5) Mark duplicates
#                        6) Indexing
#
# Author: Beatriz Vieira
#
# Date Created: 2024
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
#
# Outputs: 1) Sorted, marked duplicates and indexed varieties (.bam)
#
# Run: ./sam_to_bam.sh sample_list.txt
#
# -------------------------------------------------------
#
# Required Tool
#
# samtools (v1.7)
#
# -------------------------------------------------------

# Change to FALSE to disable option(s)
SAM_to_BAM=TRUE
BAM_stats=TRUE
BAM_mark_duplicates=TRUE
BAM_sort=TRUE
BAM_index=TRUE

# For every sample (n align_variety directories with one variety.sam file each)
# -------------------------------------------------------

while read f

do

variety="$(basename ${f%_*})"
echo $variety

# Convert SAM to BAM files
# -------------------------------------------------------
if [ $SAM_to_BAM == TRUE ]
then

inp=${variety}.sam
outp=${variety}.bam

echo "converting ${inp} to ${outp}"
cnv="samtools view -h -S -b -o  align_${variety}/${outp} align_${variety}/${inp}"
echo ${cnv}
$cnv
fi

# Compute BAM statistics
# -------------------------------------------------------

if [ $BAM_stats == TRUE ]
then

stats="samtools stats -@ 20 align_${variety}/${variety}.bam > align_${variety}/${variety}_stats.txt"
echo ${stats}
${stats}
echo "stats results for ${variety} in align_${variety}/${variety}_stats.txt"

flagstats="samtools flagstat align_${variety}/${variety}.bam > align_${variety}/${variety}_flagstats.txt"
echo ${flagstats}
${flagstats}
echo "flagstats results for ${variety} in align_${variety}/${variety}_flagstats.txt"
fi

# Add mate score tags to mark duplicates later
# -------------------------------------------------------
if [ $BAM_mark_duplicates == TRUE ]
then

fix="samtools fixmate -m align_${variety}/${variety}.bam align_${variety}/${variety}_fix.bam"
echo ${fix}
$fix
echo "fixmate results for ${variety} in align_${variety}/${variety}_fix.bam"
fi

# Sort BAM files
# -------------------------------------------------------
if [ $BAM_sort == TRUE ]
then

inbam=${variety}_fix.bam
outsorted=${variety}_sorted.bam

echo "sorting ${inbam} into ${outsorted}"

cmdsort="samtools sort -o align_${variety}/${outsorted} align_${variety}/${inbam}"
echo ${cmdsort}
$cmdsort
fi

# Mark duplicates in BAM files
# -------------------------------------------------------
if [ $BAM_mark_duplicates == TRUE ]
then

insorted=${variety}_sorted.bam
outdedup=${variety}_dedup.bam

cmdmark="samtools markdup -s align_${variety}/${insorted} align_${variety}/${outdedup}"
echo $cmdmark
$cmdmark
echo  "markdup for $variety worked"

fi

# Index BAM files
# -------------------------------------------------------
if [ $BAM_index == TRUE ]
then

indedup=${variety}_dedup.bam

echo "indexing ${indedup}"

cmdindex="samtools index align_${variety}/${indedup}"
echo ${cmdindex}
$cmdindex
echo "$variety successfully indexed in align_${variety}/${variety}_dedup.bam"

fi

done <$1
