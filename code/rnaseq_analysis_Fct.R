#######################################
## RNA-seq analysis ##
#######################################
## last update 18-Sept-2025


## source required functions
# source("code/rnaseq_analysis_Fct.R")

######################################
## Save ggplots in multiple formats ##
######################################
SavePlots <- function(plots, name, width=13, height=8){   
  ggsave(plots, filename=paste0(name,".pdf"), width = width, height = height) 
  ggsave(plots, filename=paste0(name,".png"), width = width, height = height) 
}

## Usage
# SavePlots(plt, name = "expression_matrix", width = 8, height = 4)
# SavePlots(plt, name = "analysis/output/expression_matrix")

##################################
## Check and create directories ##
##################################
check_dir <- function(dir){
  
  if( !dir.exists(dir) ){
    dir.create(dir, recursive = TRUE)
  }
}

## Usage
# check_dir("analysis/fastqc")

###############################
## Extract featurecount data ##
###############################
extractFeaturecountData <- function(file){
  data <- readr::read_delim(file,
                     id = "sample",
                     comment = "#",
                     delim = "\t",
                     col_names = TRUE) %>%
    dplyr::select(Geneid, Counts = last_col(), sample) %>%
    mutate(sample = str_extract(sample,"(?<=featurecounts/).*(?=.fcnts)"))
  
  return(data)
}

## Usage
# extractFeaturecountData("analysis/featurecounts/EM1.fcnts.txt")

######################################
## Merge all featurecount files ##
######################################
fcountMerge <- function(dir = "analysis/featurecounts"){
  list_of_files <- list.files(dir,
                              pattern = "\\.fcnts.txt$",
                              full.names = TRUE)
  
  alldata <- map_dfr(list_of_files, extractFeaturecountData)
  
  alldata <- alldata %>%
    pivot_wider(
      names_from = sample,
      values_from = Counts,
      names_sort = TRUE) %>%
    column_to_rownames("Geneid")
  
  return(alldata)
}

## Usage
# fcountdata <- fcountMerge("analysis/featurecounts")
# fcountdata

###############################################
## Obtain DESeq2 outputs for each comparison ##
###############################################
get_deseq_df <- function(cmp, dds, fct) {
  deseq_df <- data.frame(row.names = rownames(counts))
  res <- DESeq2::results(dds, contrast = c(fct, cmp))
  ## Set NAs to reasonable values to avoid errors in downstream filtering steps
  res[is.na(res[, "padj"]), "padj"] <- 1
  res[is.na(res[, "log2FoldChange"]), "log2FoldChange"] <- 0
  deg <- as.data.frame(res)
  colnames(deg)[colnames(deg) %in% c("log2FoldChange", "padj")] <- c("log2FC", "FDR")
  colnames(deg) <- paste(paste(cmp, collapse = "-"), colnames(deg), sep = "_")
  deseq_df <- cbind(deseq_df, deg[rownames(deseq_df), ])
  
  return(deseq_df)
}

## Usage
# deseqDF <- lapply(cmp, get_deseq_df, dds = dds, fct = FCT)


################################################
## Extract DESeq2 outputs for each comparison ##
################################################
extract_deseq_data <- function(cmp, deg_df, gene_list, filter = TRUE){
  data <- deg_df[, grep(cmp, colnames(deg_df)), drop = FALSE]
  if (filter == TRUE) {
    data <- data[rownames(data) %in% unlist(gene_list[cmp]), ]
    colnames(data) <- gsub(glue("{cmp}_"), "", colnames(data))
  } else {
    colnames(data) <- gsub(glue("{cmp}_"), "", colnames(data))
  }
  
  return(data)
}

## Usage
# deg_data <- lapply(colnames(pf), extract_deseq_data)
# names(deg_data) <- colnames(pf)

###############################################
## Filter DESeq2 outputs for each comparison ##
###############################################
filterDEGs <- function(degDF, FDR = 0.05, Fold = 2){

  pval <- degDF[, grep("_FDR$", colnames(degDF)), drop = FALSE]
  log2FC <- degDF[, grep("_log2FC$", colnames(degDF)), drop = FALSE]

  ## DEGs that are up or down regulated
  pf <- pval <= FDR & (log2FC >= log2(Fold) | log2FC <= -log2(Fold))
  colnames(pf) <- gsub("_FDR", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistUPorDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[,x,drop=FALSE],,drop = FALSE]), simplify = FALSE)

  ## DESeq2 results for DEGs only
  deg_data <- lapply(colnames(pf), function(current_cmp) {
    extract_deseq_data(cmp = current_cmp, deg_df = degDF, gene_list = DEGlistUPorDOWN, filter = TRUE)
  })
  names(deg_data) <- colnames(pf)

  ## DESeq2 results for all genes
  all_data <- lapply(colnames(pf), function(current_cmp) {
    extract_deseq_data(cmp = current_cmp, deg_df = degDF, gene_list = NULL, filter = FALSE)
  })
  names(all_data) <- colnames(pf)

  ## DEGs that are up regulated
  pf <- pval <= FDR & log2FC >= log2(Fold)
  colnames(pf) <- gsub("_FDR", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistUP <- sapply(colnames(pf), function(x) rownames(pf[pf[,x,drop=FALSE],,drop=FALSE]), simplify=FALSE)

  ## DEGs that are down regulated
  pf <- pval <= FDR & log2FC <= -log2(Fold)
  colnames(pf) <- gsub("_FDR", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[,x,drop=FALSE],,drop=FALSE]), simplify=FALSE)

  ## Summary dataframe
  df <- data.frame(Comparisons=names(DEGlistUPorDOWN),
                   Counts_Up_or_Down=sapply(DEGlistUPorDOWN, length), Counts_Up=sapply(DEGlistUP, length),
                   Counts_Down=sapply(DEGlistDOWN, length))

  ## summary plots
  df_plot <- data.frame(Comparisons=rep(as.character(df$Comparisons), 2), Counts=c(df$Counts_Up, df$Counts_Down), Type=rep(c("Up", "Down"), each=length(df[,1])))

  plt <- summaryPlot(df_plot, Fold = Fold, FDR = FDR)

  print(plt)

  resultlist <- list(AllData = all_data,
                     DEGData = deg_data,
                     UporDown = DEGlistUPorDOWN,
                     Up = DEGlistUP,
                     Down = DEGlistDOWN,
                     Summary = df,
                     Plot = plt)

  return(resultlist)

}

## Usage
# deseqDF <- lapply(cmp, get_deseq_df, dds, FCT)
# deseqDF <- bind_cols(deseqDF)
# deg_list <- filterDEGs(deseqDF)

############################
## GO Enrichment Analysis ##
############################
runEnrichGO <- function(ont, comp_name = compName, DF = df, orgDB = mmu) {

  print(glue("Running GO enrichment: {ont}"))
  res <- enrichGO(gene = DF$ENTREZID,
                  OrgDb = orgDB,
                  ont = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  
  print(glue("writing {ont} results to table"))
  
  write.table(res@result,
              glue::glue("{GOenrichDir}/{comp_name}.GO_enrichment.{ont}.csv"),
              sep = ",",
              row.names = F,
              quote = F)
  
  res_enrich <- dplyr::filter(res, qvalue < 0.05 & p.adjust < 0.05)
  
  # check if there are enriched results and plot
  if( dim(res_enrich)[1] > 0) {
    print("non-empty matrix: proceed to plotting")
    
    # plots of enriched results
    d_plt <- dotplot(res)
    b_plt <- barplot(res)
    c_plt <- cnetplot(res)
    
    # Exporting plots
    pdf(glue::glue("{GOenrichDir}/{comp_name}.GO_plots.{ont}.pdf"))
    my_plts <- list(d_plt, b_plt, c_plt)
    invisible(lapply(my_plts, print))
    dev.off()
  }
}

GOEnrichmentAnalysis <- function(cmp, degList, orgDB = mmu){
  
  require(clusterProfiler)
  require(org.Mm.eg.db)
  
  compName <- glue("{cmp[1]}-{cmp[2]}")
  # df <- as_tibble(degList[["DEGData"]][[compName]], rownames = "geneid")
  df <- degList[["DEGData"]][[compName]]
  df$geneid <- rownames(df)
  ids <- clusterProfiler::bitr(geneID = df$geneid, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgDB)
  df <- left_join(df, ids, by = c("geneid" = "SYMBOL"))
  df <- na.omit(df)
    
  print(glue::glue("Processing {compName}"))
  
  GO_ont <- c("BP", "MF", "CC")
  lapply(GO_ont, runEnrichGO, comp_name = compName, DF = df, orgDB = mmu)
}

## Usage
# GOEnrichmentAnalysis(cmp = cmp, degList = deg_list, orgDB = mmu)




##############################
## GSEA Enrichment Analysis ##
##############################



#############################
## ReactomePA (Hg/Mm only) ##
#############################

# reactomeAnalysis <- function(cmp,org = "mouse"){
#   
#   case_when( 
#     org == "mouse" ~ ORG <- "org.Mm.eg.db",
#     org == "human" ~ ORG <- "org.Hs.eg.db",
#     org == "tair" ~ ORG <- "org.At.tair.db",
#     TRUE ~ "Unknown Genome"
#   )
#   
#   df <- read_csv(file.path(dataDir, sampleName, "diff_expression.csv"))
#   sig_prots <- df |> dplyr::filter(abs(log2_foldchange) >= 0.58 & p_value < 0.05)
#   ids <- bitr(sig_prots$protein, fromType = "SYMBOL", toType = c("ENTREZID", "SYMBO
# L"), OrgDb = ORG)
#   sig_prots <- left_join(sig_prots, ids, by = c("protein" = "UNIPROT"))
#   sig_prots <- dplyr::select(sig_prots, c(ENTREZID, UNIPROT = protein, SYMBOL, log2_
#                                           foldchange)) |> distinct()
#   
#   x <- enrichPathway(sig_prots$ENTREZID, pvalueCutoff = 0.05, readable = TRUE, organism = "mouse")
#   
#   outDir <- file.path(dataDir,sampleName,"ReactomePA")
#   
#   if(!dir.exists(outDir)){
#     dir.create(outDir)}
#   
#   write_csv(x@result, file.path(outDir,"ReactomePA.csv"), col_names = TRUE, quote = "needed")
# 
#   print("Plotting top 20 pathways")
#   
#   if ( length(x@result$Description) < 11 ) {
#     len <- length(x@result$Description)
#   } else {
#     len <- 10
#   }
#   
#   for(i in 1:len){
#     
#     cat("Processing",i," of ", len, "\n")    
#     cln_name <- gsub(" ", "_", x@result$Description[i])
#     cln_name <- gsub("/", "-", cln_name)
#     plt <- viewPathway(x@result$Description[i], organism = "mouse", readable = TRUE)
#     ggsave(file.path(outDir, glue("{cln_name}.pdf")), plot = plt, width = 8, height = 11)
#   }
# }



###################################################################
######################### VISUALIZATION ###########################
###################################################################

##################################################
## Summary Plot of Up/Down DEGs per comparisons ##
##################################################
summaryPlot <- function(df, Fold = 2, FDR = 0.05){
  df |>
    dplyr::mutate(Counts = if_else(Type == "Down", -1*Counts, Counts)) |>
    
  ggplot(aes(x = Comparisons, y = Counts, fill = Type)) +
  geom_col() +
  geom_hline(yintercept = 0, color = "black") +
  labs(y = "Number of DEGs",
       x = "",
       title = glue("DEG Counts (Fold = {Fold} & FDR = {FDR})")) +
  coord_flip() +
  facet_grid(~Type, scales = "free_x") +
  scale_y_continuous(expand = c(0,0),
                     labels = function(x) abs(x)) +
  theme(legend.position = "none",
        panel.spacing.x = unit(0, "pt"),
        strip.background = element_rect(color = "black", fill = "grey"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))
}

## Usage
# plt <- SummaryPlot(df_plot)
# plt

############################
## Enhanced Volcano Plots ##
############################
eVolcanoPlot <- function(cmp, dds, fct){
  res <- DESeq2::results(dds, contrast = c(fct, cmp))
  cond1 <- cmp[[1]]
  cond2 <- cmp[[2]]
  
  title <- glue("Factor: {fct} \n Comparison: {cond1} vs {cond2}") 
  subtitle <- paste("LFC > 1, FDR < 0.05")
  
  plt <- EnhancedVolcano(res,
                         lab = rownames(res),
                         x = 'log2FoldChange',
                         y = 'pvalue',
                         title = title,
                         subtitle = subtitle,
                         pCutoff = 0.05,
                         pCutoffCol = 'padj', # p-value
                         FCcutoff = 0,
                         pointSize = 4.0,
                         labSize = 6.0,
                         drawConnectors = TRUE,
                         widthConnectors = 0.75)
  
  SavePlots(plt, name = glue("{plotsDir}/Enhancedvolcano_{cond1}-vs-{cond2}"))
}

## Usage
# lapply(cmp, eVolcanoPlot, dds, FCT)


########################
## PCA Plots (Plotly) ##
########################
# modified 11/26/2025 - Added parameters to plot batch effect-removed data
pcaPlot <- function(dds, fct, outfile_prefix="PCA", batch_effect=FALSE, batch_condition="batch") {
  
  vst <- vst(dds)
  if( batch_effect == TRUE){
  # rld <- rlog(dds)
    vst_mat <- assay(vst)
    vst_corrected <- limma::removeBatchEffect(vst_mat, batch = vst[[batch_condition]])
    vst_cor_obj <- vst
    assay(vst_cor_obj) <- vst_corrected
    pltdata <- plotPCA(vst_cor_obj, intgroup = fct, returnData = TRUE)
    percentVar <- round(100 * attr(pltdata, "percentVar"))
    plt <- ggplot(pltdata, aes(PC1, PC2, color=group, label=name)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_fixed()
    pltly <- plotly::ggplotly(plt)
    htmlwidgets::saveWidget(pltly, glue("{plotsDir}/{outfile_prefix}_batchcorrected.html"))
    SavePlots(plt, glue("{plotsDir}/{outfile_prefix}_batchcorrected"))
  }  
    
  pltdata <- plotPCA(vst, intgroup = fct, returnData = TRUE)
  percentVar <- round(100 * attr(pltdata, "percentVar"))
  plt <- ggplot(pltdata, aes(PC1, PC2, color=group, label=name)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed()
  pltly <- plotly::ggplotly(plt)
  htmlwidgets::saveWidget(pltly, glue("{plotsDir}/{outfile_prefix}.html"))
  SavePlots(plt, glue("{plotsDir}/{outfile_prefix}"))
}

## Usage
# pcaPlot(dds, FCT, outfile_prefix="PCA", batch_effect=TRUE, batch_condition="batch")
# pcaPlot(dds, FCT, outfile_prefix="PCA", batch_effect=FALSE)

#################
## UpSet Plots ##
#################
upsetPlot <- function(deg_list, plotsDir) {
  
  plot_list <- c("UporDown", "Up", "Down")
  
  for(degs in plot_list) {
    
    if(length(names(deg_list[[degs]])) > 1){
      my_plot <- UpSetR::upset(UpSetR::fromList(deg_list[[degs]]),
                               order.by = "freq",
                               nsets = length(names(deg_list[[degs]])),
                               mainbar.y.label = "Number of DEGs Overlapped",
                               sets.x.label = "Number of DEGs",
                               group.by = "sets",
                               mb.ratio = c(0.7,0.3)
      ) 
      
      pdf(file = glue("{plotsDir}/upset_{degs}.pdf"))
      print(my_plot)
      dev.off()
      
      png(file = glue("{plotsDir}/upset_{degs}.png"))
      print(my_plot)
      dev.off()
    }
  }  
}

## Usage
# upsetPlot(deg_list, plotsDir)

#################
## Heatmaps ##
#################
heatmapDEG <- function(cmp, dds, fct){
  
  cond1 <- cmp[[1]]
  cond2 <- cmp[[2]]
  
  # extract VST matrix from dds
  vsd <- vst(dds, blind = FALSE)
  expr_mat <- assay(vsd)
  
  # obtain up to top 25 DEGs 
  res <- DESeq2::results(dds, contrast = c(fct, cmp))
  res <- res[complete.cases(res$padj), ]
  deg <- subset(res, padj < 0.05 & abs(log2FoldChange) >= 1 )
  deg <- deg[order(deg$log2FoldChange, decreasing = TRUE),]
  up_genes <- rownames(head(deg, 25))
  dn_genes <- rownames(tail(deg, 25))
  sig_genes <- unique(c(up_genes, dn_genes))
  
  # Subset matrix with DEGs and samplename
  if(length(sig_genes) > 0){
    sig_genes <- sig_genes[sig_genes %in% rownames(expr_mat)]
    samplesToKeep <- rownames(sampleinfo[sampleinfo$factor %in% c(cond1, cond2),])
    subset_mat <- expr_mat[sig_genes, samplesToKeep, drop = FALSE]
    # subset_mat <- as.matrix(subset_mat)
    
    if(nrow(subset_mat) == 0){
      message("No valid genes for heatmap.")
      return(NULL)
    }
    
    subset_mat <- subset_mat[rowSums(is.finite(subset_mat)) == ncol(subset_mat), , drop = FALSE]
    # plotting heatmap
    cluster_rows = FALSE
    if( nrow(subset_mat) >= 2){
      cluster_rows = TRUE
    }
    
    pheatmap(subset_mat,
             cluster_rows = cluster_rows,
             show_rownames = TRUE,
             cluster_cols = FALSE,
             scale = "row",
             # cutree_rows = 2,
             # cutree_cols = 2,
             # annotation_col = col_annotation,
             # annotation_row = row_annotation,
             main = glue("{cond1} vs {cond2}"),
             filename = glue("{plotsDir}/Heatmap_top25DEGs_{cond1}-vs-{cond2}.pdf")
             )
  }
}

## Usage
# lapply(cmp, heatmapDEG, dds, FCT)

#################################################################
######################### DATA EXPORT ###########################
#################################################################


#######################
## Export gene lists ##
#######################
list_export <- function(deg_list, annotations, genelistsDir){
  
  print(glue("Printing to outDir: {genelistsDir}"))
  deg_set <- c("Up", "Down")
  
  for(deg in deg_set){
    for(cmp in names(deg_list[[deg]])){
      genelist <- deg_list[[deg]][[cmp]] %>% as_tibble()
      colnames(genelist) <- "Gene ID"
      
      if(!is.null(deg_list[[deg]][[cmp]])){
        genelist <- genelist %>%
          left_join(annotations, by = c("Gene ID" = "SYMBOL"))
      }
      write_csv(genelist, 
                file = glue("{genelistsDir}/{cmp}.{deg}DEG.annotated.csv"),
                quote = "needed",
                col_names = TRUE
      )
    }
  }
  
  for(cmp in names(deg_list[["DEGData"]])){
    
    degdata <- deg_list[["DEGData"]][[cmp]] %>%
      as_tibble(rownames = "gene_id")
    
    write_csv(degdata,
              file = glue("{genelistsDir}/{cmp}.deseq2_output.csv"),
              quote = "needed",
              col_names = TRUE)
    
  }
  
}

## Usage
# mmu <- org.Mm.eg.db
# annot <- clusterProfiler::bitr(keys(mmu), fromType = "ENTREZID", toType = c("SYMBOL", "GENENAME", "GENETYPE"), OrgDb = mmu)
# deseqDF <- lapply(cmp, get_deseq_df, dds, FCT)
# deseqDF <- bind_cols(deseqDF)
# deg_list <- filterDEGs(degDF = deseqDF)
# list_export(deg_list, annot, genelistsDir)

