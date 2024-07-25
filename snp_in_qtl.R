## ---------------------------------------------------------
##
## Script name: snp_in_qtl.R
##
## Purpose of the script: Attributes QTL labels to SNPs, if the SNP position
##                        occurs within one or more of the QTL regions provided.
##
## Author: Hugo Rodrigues
##
## Date Created: 2024-01-05
##
## Copyright (c) Hugo Rodrigues, 2024
## Email: hm.rodrigues@itqb.unl.pt
##
## Whole-genome polymorphisms and phylogeny of rice varieties circulating in the European market
## Maria Beatriz, Hugo Rodrigues, Pedro Barros, Margarida Oliveira
##
## ---------------------------------------------------------
##
## Inputs: 1) SNP table (tab-delimited)
##         2) QTL table (tab-delimited)
##
## Output: 1) SNP table with QTL information added
##
## ---------------------------------------------------------

## Set working directory

setwd("C:/Users/TRACE_Rice/Documents/Projects/RScripts")     # Change working directory here

## ---------------------------------------------------------

# Load Packages

library(data.table)
library(tidyverse)

## ---------------------------------------------------------

# Read Tables

sample_list <- c("Chr", "Pos", "ID", "Ref", "Alt", "Qual", "Pass", "SnpEff", "Genotype",
                 "Albatros", "Arborio", "Arelate", "Ariete", "Basmati_TypeIII",
                 "Bomba", "CL28", "Caravela", "Carnaroli", "Elettra", "Gageron",
                 "Giza177", "Giza181", "JSendra", "Lusitano", "Macarico",
                 "Manobi", "Puntal", "Ronaldo", "Super_Basmati", "Teti", "Ulisse")

SNP_table <- read_delim("./rice_high_impact_variants.tab", delim="\t", col_names = sample_list)

EQ_QTLs <- read_delim("./eatingqualityQTL.txt", col_names = c("Chr", "sPos", "ePos", "EQTL"))     # Eating quality QTLs
EQ_QTLs <- EQ_QTLs[,1:4]

seed_QTLs <- read_delim("./seedQTL.txt", col_names = c("Chr", "sPos", "ePos", "SQTL"))     # Seed related QTLs
seed_QTLs <- seed_QTLs[,1:4]

# Prepare output table header

output <- data.frame(Chr=integer(), Pos=integer(), Eating_Quality_QTL=character(),
                     Seed_QTL=character(), Annotation=character())

# Define a list for the multiple chromosomes, SNP positions and annotations

chrs <- SNP_table[,1]
snps <- SNP_table[,2]
anns <- SNP_table[,8]

# Loop for each individual SNP in data (Pos column)

for (i in 1:nrow(snps)) {
  # For each SNP in SNP_table, get its chromosome (chr), position (pos) and annotation (ann)
  chr <- chrs[i,]
  chr <- chr[[1]]
  snp <- snps[i,]
  snp <- snp[[1]]
  ann <- anns[i,]
  ann <- ann[[1]]
  # For the current SNP, retrieve its associated Eating Quality QTLs
  booleans_EQ_QTLs <- data.table::between(snp, EQ_QTLs$sPos, EQ_QTLs$ePos)
  filtered_EQ_QTLs <- EQ_QTLs[booleans_EQ_QTLs,]
  filtered_EQ_QTLs <- toString(filtered_EQ_QTLs$EQTL)
  # For the current SNP, retrieve its associated Seed related QTLs
  booleans_seed_QTLs <- data.table::between(snp, seed_QTLs$sPos, seed_QTLs$ePos)
  filtered_seed_QTLs <- seed_QTLs[booleans_seed_QTLs,]
  filtered_seed_QTLs <- toString(filtered_seed_QTLs$SQTL)
  # Add SNP row to output table
  output[nrow(output) + 1,] = c(chr, snp, filtered_EQ_QTLs, filtered_seed_QTLs, ann)
}

# Process output Annotation column

split <- within(output, Annotation<-data.frame(do.call('rbind', strsplit(as.character(Annotation), '|', fixed=TRUE))))

core <- split[, c("Chr", "Pos", "Eating_Quality_QTL", "Seed_QTL")]
Modification <- split$Annotation$X2
Gene_Name <- split$Annotation$X4
Gene_ID <- split$Annotation$X5
Transcript_ID <- split$Annotation$X7

final_output <- cbind(Gene_ID, Gene_Name, Transcript_ID, core, Modification)
final_output <- final_output[, c("Gene_ID", "Gene_Name", "Transcript_ID", "Chr", "Pos",
                                 "Modification", "Eating_Quality_QTL", "Seed_QTL")]

# Write final output table

write.csv(final_output, file = "SNP_QTL_TABLE.csv", row.names = FALSE)
