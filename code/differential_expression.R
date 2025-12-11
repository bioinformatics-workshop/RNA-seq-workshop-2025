#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

## Load helper functions
source("code/rnaseq_analysis_Fct.R")

## load libraries
library(org.Mm.eg.db)
library(tidyverse)
library(DESeq2)
library(glue)
library(EnhancedVolcano)
library(plotly)
library(clusterProfiler)
library(pheatmap)

# FCT <- args[1]
FCT <- "factor"

## Create output directory
# if( !dir.exists("analysis/deseq2") ){
#   dir.create("analysis/deseq2", recursive = TRUE)
# }

## Load annotation file
mmu <- org.Mm.eg.db
annot <- clusterProfiler::bitr(keys(mmu), fromType = "ENTREZID", toType = c("SYMBOL", "GENENAME"), OrgDb = mmu)

## Load sample info and metadata
sampleinfo <- read_csv("metadata/metadata.csv", col_names = TRUE) |>
  column_to_rownames("samplename")
sample_order <- rownames(sampleinfo)

## Load count data
counts <- fcountMerge()
counts <- counts[,sample_order]

print(glue("Processing factor: {FCT}"))

print("generating comparisons")
# get a list of all possible combinations from sample groupings
cmp <- gtools::combinations(n = length(unique(sampleinfo[,FCT])),
                            r = 2,
                            v = sampleinfo[,FCT])
cmp <- as.list(data.frame(t(cmp)))

print("create output directories")
# create output directories
genelistsDir <<- glue("analysis/deseq2/{FCT}/genelists")
plotsDir <<- glue("analysis/deseq2/{FCT}/plots")
GOenrichDir <- glue("analysis/deseq2/{FCT}/GOenrichment") 
outDir <- list(genes=genelistsDir, plots=plotsDir, GO=GOenrichDir)

lapply(outDir, check_dir)

print(glue("Output directories: {genelistsDir}"))
print(glue("Output directories: {plotsDir}"))
print(glue("Output directories: {GOenrichDir}"))

## Generate DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                colData = as.matrix(sampleinfo),
                                design = formula(glue("~ {FCT}")))
dds <- DESeq(dds) # run differential analysis for all comparisons

# extracting comparison results
deseqDF <- lapply(cmp, get_deseq_df, dds, FCT)
deseqDF <- bind_cols(deseqDF)
 
deg_list <- filterDEGs(deseqDF)
saveRDS(deg_list, file = glue("{genelistsDir}/deseq2_results.RDS"))

print("exporting data")
print(deg_list$Summary)

list_export(deg_list, annot, genelistsDir) # DEG lists
SavePlots(plots = deg_list[["Plot"]], name = glue("{plotsDir}/DEG_summary")) # summary plot
pcaPlot(dds, FCT) # PCA plots
upsetPlot(deg_list, plotsDir) # upset plots (if there are more than one comparisons)
lapply(cmp, eVolcanoPlot, dds, FCT) # enhanced volcano plots
lapply(cmp, heatmapDEG, dds, FCT) # heatmaps
# Enrichment analysis
lapply(cmp, GOEnrichmentAnalysis, degList = deg_list)










