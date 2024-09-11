## ---------------------------------------------------------
##
## Script name: heatmaps.R
##
## Purpose of the script: Get a SNP distribution heatmap for all chromosomes in all samples
##
## Author: Hugo Rodrigues
##
## Date Created: 2024-02-21
##
## Copyright (c) Hugo Rodrigues, 2024
## Email: hm.rodrigues@itqb.unl.pt
##
## Whole-genome polymorphisms and phylogeny of rice varieties circulating in the European market
## Maria Beatriz, Hugo Rodrigues, Pedro Barros, Margarida Oliveira
##
## ---------------------------------------------------------
##
## Inputs: 1) SNP density table (tab-delimited)
##
## Output: 1) SNP density heatmaps
##
## ---------------------------------------------------------

## Set working directory

setwd("C:/Users/TRACE_Rice/Documents/Projects/RScripts")     # Change working directory here

## ---------------------------------------------------------

# Load Packages

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(data.table)
library(ldt)
library(tidyverse)

## ---------------------------------------------------------

# Read SNP density Table

data <- read_delim("./density.out.by_sample", delim="\t")

# Chromosomes 1 to 12

chrs <- data[data$Chromosome <= 12, ]

# Add 1 and log2 all matrix values

log2p1 <- function(x) log1p(x)/log(2)

# Make chromosome matrix heatmap

make_matrix <- function(chr_name, chr_number) {
  chr_name <- chrs[chrs$Chromosome == chr_number, ]
  chr_name <- subset(chr_name, select = -Chromosome)
  chr_name <- as.matrix(chr_name)
  chr_name <- t(chr_name)
  chr_log <- log2p1(chr_name)
  return(chr_log)
}

chr1_log <- make_matrix(chr1, 1)
chr2_log <- make_matrix(chr2, 2)
chr3_log <- make_matrix(chr3, 3)
chr4_log <- make_matrix(chr4, 4)
chr5_log <- make_matrix(chr5, 5)
chr6_log <- make_matrix(chr6, 6)
chr7_log <- make_matrix(chr7, 7)
chr8_log <- make_matrix(chr8, 8)
chr9_log <- make_matrix(chr9, 9)
chr10_log <- make_matrix(chr10, 10)
chr11_log <- make_matrix(chr11, 11)
chr12_log <- make_matrix(chr12, 12)

# Coloring de-gra-de

col_fun2 = colorRamp2(c(0, 10, 15), c("lightblue", "blue", "darkblue"))

# Generating heatmap plots

plot_heatmap <- function(heatmap, title) {
  Heatmap(heatmap, name = "log(SNPs)", col = col_fun2, cluster_columns = FALSE,
          row_title = "Varieties", show_column_dend = FALSE, row_names_side = "left",
          column_title = str(title), row_km = 4, row_dend_width = unit(2, "cm"),
          row_title_gp = gpar(fontsize = 24), column_title_gp = gpar(fontsize = 24))
}

# After each plot, save image with using Export (Plots tab) > Save as Image > SVG format

plot_heatmap(chr1_log, "Chr1")
plot_heatmap(chr2_log, "Chr2")
plot_heatmap(chr3_log, "Chr3")
plot_heatmap(chr4_log, "Chr4")
plot_heatmap(chr5_log, "Chr5")
plot_heatmap(chr6_log, "Chr6")
plot_heatmap(chr7_log, "Chr7")
plot_heatmap(chr8_log, "Chr8")
plot_heatmap(chr9_log, "Chr9")
plot_heatmap(chr10_log, "Chr10")
plot_heatmap(chr11_log, "Chr11")
plot_heatmap(chr12_log, "Chr12")
