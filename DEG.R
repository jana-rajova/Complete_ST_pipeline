if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install("DESeq2")
BiocManager::install("ArrayExpress")
install.packages("dplyr")
install.packages("pheatmap")
install.packages("ggrepel")
install.packages("RColorBrewer")
BiocManager::install("limma")
BiocManager::install("edgeR")
install.packages("enrichR")
install.packages("gridExtra")
install.packages("stringr")

library(dplyr)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(limma)
library(edgeR)
library(enrichR)
library(gridExtra)
library(stringr)

dataset <- "CN56"
methods <- c('seurat', 'DESC', 'scanorama')
path <- "/media/MNM-NetStorage/OrganizedSeqFiles/ST/Complete_ST_pipeline/results/cluster_gene_expression/"
i <- 'DESC'
for (i in methods){
  countTable <-read.table(paste0(path, dataset, '/', dataset, '_', i, '_cluster_expression.tsv'), sep='\t', header=TRUE)
  print(dataset)
  print(i)
  countTable$cluster <- as.factor(countTable$cluster)
  by_cluster <- countTable %>% group_by(cluster)
  TH.top.cluster <- by_cluster %>% summarise(TH=sum(TH))
  TH.top.cluster <- TH.top.cluster[order(TH.top.cluster$TH, decreasing = TRUE),]
  print(TH.top.cluster)
  TH.top.cluster <- as.character(TH.top.cluster$cluster)
  x=''
  for (j in TH.top.cluster[-1]){
    x <- c(x, paste0('cluster', TH.top.cluster[1],'-cluster',j))
  }
  x <- x[-1]
  
  rownames(countTable) <- paste0(countTable$well, '-', countTable$cluster)
  sampleTable <- data.frame('cluster' = countTable$cluster, 'well' = countTable$well)
  rownames(sampleTable) <- rownames(countTable)
  countTable <- t(countTable)
  countTable <- countTable[-c(1,2),]
  mode(countTable) <- 'numeric'
  
  dds <-DESeqDataSetFromMatrix(as.matrix(countTable),design =~0+cluster,colData=sampleTable)
  print(dds)
  #Normalize
  normCounts <-vst(dds, blind =TRUE)
  assay(normCounts)[1:5, 1:5]
  #Distribution
  hist(assay(normCounts), main=paste0(dataset, '-', i))
  
  #Step 1: Define design matrix
  designMatrix <-model.matrix(~0+cluster, data =sampleTable)
  head(designMatrix)
  contrastMatrix <-makeContrasts(contrasts=x ,levels =designMatrix)
  head(contrastMatrix)

  dge <-DGEList(countTable)
  dge <-calcNormFactors(dge)
  dge <-estimateDisp(dge, designMatrix, robust =TRUE)
  fit <-glmQLFit(dge, designMatrix, robust =TRUE)
  
  #Step 4: Perform hypothesis testing
  res <-glmQLFTest(fit, contrast =contrastMatrix)
  res <-topTags(res, n =nrow(countTable))
  sigRes <- res$table[1:13]
  sigRes$logFC.ave <- rowMeans(sigRes)
  sigRes$FDR <- res$table$FDR
  sigRes.sub <-subset(sigRes, FDR <0.5&abs(logFC.ave)> 1)
  write.table(sigRes.sub, paste0(path, dataset, '/', 'DEG_', i, '.tsv'), sep="\t")
  head(sigRes.sub, 10)
  nrow(sigRes.sub)
  
  #Visualize results
  volcanoPlot <-ggplot(sigRes,
                       aes(x =logFC.ave,
                           y =-log10(FDR),
                           color =ifelse(FDR <0.05&abs(logFC.ave) >1, "darkred", "grey"))) +
    geom_point() +
    xlab(expression("Fold Change, Log"[2]*"")) +
    ylab(expression("Adjusted P value, Log"[10]*"")) +
    geom_vline(xintercept =c(-1, 1), linetype ="dotted", size =1) +
    geom_hline(yintercept =-log10(0.05), linetype ="dotted", size =1) +
    theme_minimal() +theme(legend.position ="none") +
    scale_colour_manual(values =c("darkred", "grey", "steelblue")) +
    geom_text_repel(aes(x =logFC.ave, y =-log10(FDR), 
                        label = rownames(sigRes.sub)[1:25],
                        size =2, 
                        color ="steelblue"),
                    data =sigRes.sub[1:25, ])
  jpeg(paste0(path, dataset, '/volcano_plot_', i,'.png'))
  print(volcanoPlot)
  dev.off()
  print(volcanoPlot)

    
  #Enrichment analysis
  head(listEnrichrDbs())
  enrichmentRes <-enrichr(
    genes =rownames(sigRes.sub),
    databases =c("Allen_Brain_Atlas_10x_scRNA_2021", "KEGG_2019_Human"))
  # check entries from the two databases 
  enrichmentRes$Allen_Brain_Atlas_10x_scRNA_2021[1:5, 1:4]
  enrichmentRes$KEGG_2019_Human[1:5, 1:4]
  write.table(enrichmentRes$KEGG_2019_Human, paste0(path, dataset, '/', 'KEGG_2019_human_GO_', i, '.tsv'), sep="\t")
  write.table(enrichmentRes$Allen_Brain_Atlas_10x_scRNA_2021, paste0(path, dataset, '/', 'Allen_Brain_atlas_GO_', i, '.tsv'), sep="\t")
  
  #Visualize significant terms
  enrichPlots <-list()
  for(res in names(enrichmentRes)){
    plotDF <-enrichmentRes[[res]][1:10, ]
    plotDF$Term =str_trunc(plotDF$Term, 45)
    plotDF$Term <-factor(plotDF$Term, levels =rev(plotDF$Term))
    enrichPlot <-ggplot(plotDF, aes(x =Term, y =-log10(Adjusted.P.value))) +
      geom_bar(stat ="identity", width =0.05) +geom_point(size =3) +
      theme_minimal() +
      theme(text =element_text(size =10),
            plot.title =element_text(hjust =(10/nchar(res)) *2),
            plot.margin =margin(t =5, r =50, b =5, unit ="pt"),
            axis.text.y =element_text(size =8)) +
      coord_flip() +geom_hline(yintercept =-log10(0.05), linetype ="dashed", color ="gray") +
      labs(title =res, y =expression("Adjusted P value, Log"[10]*""))
    enrichPlots[[res]] <-enrichPlot }
  plot(arrangeGrob(grobs =enrichPlots, layout_matrix =rbind(c(1, 2))))
}
