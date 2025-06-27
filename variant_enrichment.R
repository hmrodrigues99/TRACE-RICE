## ---------------------------------------------------------
##
## Script name: variant_enrichment.R
##
## Purpose of the script: Performs a gene enrichment analysis on a set of
##                        genes previously annotated to contain HIGH impact SNPs.
##
## Author: Hugo Rodrigues
##
## Date Created: 2024-02-08
##
## Copyright (c) Hugo Rodrigues, 2024
## Email: hm.rodrigues@itqb.unl.pt
##
## Whole-genome polymorphisms and phylogeny of rice varieties circulating in the European market
## Maria Beatriz, Hugo Rodrigues, Pedro Barros, Margarida Oliveira
##
## ---------------------------------------------------------
##
## Inputs: 1) SNP table (tab-delimited) previously annotated using SnpEff
##
## Output: 1) Plot highlighting enriched GO:terms for the provided genes
##
## ---------------------------------------------------------

## Set working directory

setwd("C:/Users/TRACE_Rice/Projects/RScripts")     # Change working directory here

## ---------------------------------------------------------

# Load Packages

library(readr)
library(gprofiler2)
library(ggplot2)
library(svglite)

## ---------------------------------------------------------

# Read HIGH-impact gene list

genes <- read_lines(file = "genes_high_impact.txt")
genes <- unique(genes)

# Get gprofiler2 enrichment results
gostres <- gost(query = genes, organism = "osativa", significant = TRUE, user_threshold = 0.05)
data <- gostres$result
# Get top 15 enriched terms
data <- head(data, 15)
# Removal of redundant terms
data <- data[-c(2,5,6,7,9,12,13,15), ]

# Make plot
p <- ggplot(data, aes(x = -log2(p_value), y = reorder(term_name, p_value))) +
  geom_point(aes(size = intersection_size, color = source), stroke = 1) +
  scale_size_continuous(range = c(1,10)) +
  labs(x = "pvalue", y = "GO Term", size = "Count", color = "Category", title = "High Impact Genes Enrichment")

# Add theme
p <- p + theme(text = element_text(size = 15), axis.text.x = element_text(size = 10))

# Save figure
ggsave(filename = "enrichment_analysis.svg", width = 10, height = 8, dpi = 600)
