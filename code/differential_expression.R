#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# Load helper functions
source("code/utility_scripts.R")

# load libraries
library(tidyverse)
library(DESeq2)
library(glue)
library(EnhancedVolcano)
library(plotly)


FCT = args[1]

# Load annotation file
annot <- read_csv("genome/ath_gene_annot.csv")

# Load sample info and metadata
sampleinfo <- read_csv("metadata/metadata.csv", col_names = TRUE)

sampleinfo <- sampleinfo %>%
    column_to_rownames("samplename")

sample_order <- rownames(sampleinfo)

# Load count matrix file from featurecounts and clean up the header
counts <- fcount_merge()

counts <- counts[,sample_order]

# Obtain deseq2 output per factor comparisons
deseq_per_factor <- function(fct){

    print(glue("Processing factor: {fct}"))

    print("generating comparisons")
    # get a list of all possible combinations from sample groupings
    cmp <- gtools::combinations(n = length(unique(sampleinfo[,fct])),
            r = 2,
            v = sampleinfo[,fct])
    cmp <- as.list(data.frame(t(cmp)))

    print("create output directories")
    # create output directories
    genelistsDir <<- glue("analysis/deseq2/{fct}/genelists")
    plotsDir <<- glue("analysis/deseq2/{fct}/plots")
    outDir <- list(genes=genelistsDir, plots=plotsDir)

    lapply(outDir, check_dir)

    print(glue("Output directories: {genelistsDir}"))
    print(glue("Output directories: {plotsDir}"))

    print("generating dds object")
    # Construct DESeqDataSet
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                colData = as.matrix(sampleinfo),
                                design = formula(glue("~ {fct}")))

    dds

    print("running deseq2")
    # Run DESeq2
    dds <- DESeq(dds)

    print("generating PCA")
    generate_pca(dds, fct)

    print("extracting comparisons results")
    # Extract results
    deseqDF <- lapply(cmp, get_deseq_df, dds, fct)
    deseqDF <- bind_cols(deseqDF)

    class(deseqDF)

    print("generating filtered gene lists")
    ## Export filtered output
    deg_list <- filterDEGs(degDF = deseqDF)
    saveRDS(deg_list, file = glue("{genelistsDir}/deseq2_results.RDS"))

    ## Exporting data

    print("exporting list")
    print(deg_list$Summary)
    # export DEG lists
    list_export(deg_list, annot, genelistsDir)

    print("exporting tables and plots")
    # export summary table
    write.table(deg_list$Summary,
        file = glue("{genelistsDir}/deg_summary.csv"),
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE)

    # export summary plot
    SavePlots(plots = deg_list[["Plot"]], 
        name = glue("{plotsDir}/DEG_summary"))

    # export upset plots (if there are more than one comparisons)
    upset_plots(deg_list, plotsDir)

    # export enhanced volcano plots
    lapply(cmp, e_volcano_plot, dds, fct)

}

deseq_per_factor(FCT)