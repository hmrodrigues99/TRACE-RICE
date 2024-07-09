## ---------------------------------------------------------
##
## Script name: Variant_Enrichment.R
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
## Output: 1) Plot highliting enriched GO:terms for the provided genes
##
## ---------------------------------------------------------

## Set working directory

setwd("C:/Users/TRACE_Rice/Documents/Projects/RScripts")     # Change working directory here

## ---------------------------------------------------------

# Load Packages

library(data.table)
library(gprofiler2)
library(stringr)
library(tidyverse)

## ---------------------------------------------------------

# Read High Impact Variant Table

data <- read_delim("./rice_high_impact_variants.tab", delim="\t", col_names = c("Chr", "Pos", "ID", "Ref", "Alt", "Qual", "Pass", "SnpEff", "Genotype", 
                                                                                "Albatros", "Arborio", "Arelate", "Ariete", "Basmati_TypeIII", "Bomba",
                                                                                "CL28", "Caravela", "Carnaroli", "Elettra", "Gageron", "Giza177", "Giza181",
                                                                                "JSendra", "Lusitano", "Macarico", "Manobi", "Puntal", "Ronaldo",
                                                                                "Super_Basmati", "Teti", "Ulisse"))

# Get gene list
annotations <- data['SnpEff']
split_annotations <- as.data.frame(str_split_fixed(annotations$SnpEff, '\\|', 6))
genes <- split_annotations[,5]
genes <- unique(genes)

# Gene Enrichment Analysis (gprofiler2)
gostres <- gost(query = genes, organism = "osativa", significant = TRUE)

# Check top 3 GO:terms
head(gostres$result, 3)

# Uncomment for interactive use
# gostplot(gostres, capped = TRUE, interactive = TRUE)

# Generate Manhattan plot of the functional enrichment results
p <- gostplot(gostres, capped = TRUE, interactive = FALSE)
publish_gostplot(p, filename = "./HIGH_genes_top_enrichment_terms.png", width = 10, height = 6, 
                 highlight_terms = c("GO:0006952", "GO:0044419", "GO:0043207"))
