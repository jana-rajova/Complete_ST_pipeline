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

datasets <- c("CN56", "SN-c", "TX", "SN-TX-c")
methods <- c('seurat', 'DESC', 'scanorama')
path <- "/media/MNM-NetStorage/OrganizedSeqFiles/ST/Complete_ST_pipeline/results/cluster_gene_expression/"

gene <- 'TH'
for (dataset in datasets){
  for (i in methods){
    countTable <-read.table(paste0(path, dataset, '/', dataset, '_', i, '_cluster_expression.tsv'), sep='\t', header=TRUE)
    print(dataset)
    print(i)
    countTable$cluster <- as.factor(countTable$cluster)
    by_cluster <- countTable %>% group_by(cluster)
    gene.top.cluster <- by_cluster %>% summarize(TH=mean(TH/feature_counts))
    hist.TH <- as.data.frame(gene.top.cluster)
    hist.TH <- hist.TH[order(hist.TH$TH, decreasing = TRUE),]
    
    gene.top.cluster <- gene.top.cluster[order(gene.top.cluster$TH, decreasing = TRUE),]
    jpeg(paste0(path, dataset, '/', dataset, '_top_TH_cluster_plot_', i,'.png'))
    barplot(height=hist.TH$TH, names=hist.TH$cluster, 
            col=heat.colors(length(hist.TH$TH)),
            xlab="cluster", 
            ylab="TH count"
    )
    dev.off()
    
    write.table(gene.top.cluster, paste0(path, dataset, '/', 'TH_top_clusters_', i, '.tsv'), sep="\t", quote = FALSE)
    MBP.top.cluster <- by_cluster %>% summarise(MBP=mean(MBP/feature_counts))
    hist.MBP <- as.data.frame(MBP.top.cluster)
    hist.MBP <- hist.MBP[order(hist.MBP$MBP, decreasing = TRUE),]
    jpeg(paste0(path, dataset, '/', dataset, '_top_MBP_cluster_plot_', i,'.png'))
    barplot(height=hist.MBP$MBP, names=hist.MBP$cluster, 
                       col=heat.colors(length(hist.MBP$MBP)),
                       xlab="cluster", 
                       ylab="MBP count"
    )
    
    dev.off()
    write.table(gene.top.cluster, paste0(path, dataset, '/', 'MBP_top_clusters_', i, '.tsv'), sep="\t", quote = FALSE)
    print(gene.top.cluster)
    gene.top.cluster <- as.character(gene.top.cluster$cluster)
    x=''
    for (j in gene.top.cluster[-1]){
      x <- c(x, paste0('cluster', gene.top.cluster[1],'-cluster',j))
    }
    x <- x[-1]
    x
    rownames(countTable) <- paste0(countTable$well, '-', countTable$cluster)
    sampleTable <- data.frame('cluster' = countTable$cluster, 'well' = countTable$well)
    rownames(sampleTable) <- rownames(countTable)
    countTable <- t(countTable)
    countTable <- countTable[-c(1:4),]
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
    sigRes <- res$table[1:length(x)]
    sigRes$logFC.ave <- rowMeans(sigRes)
    sigRes$FDR <- res$table$FDR
    sigRes.sub <-subset(sigRes, FDR <0.01&abs(logFC.ave)> 1)
    sigRes.sub$log.FDR <- sigRes.sub$logFC.ave / -log10(sigRes.sub$FDR)
    write.table(sigRes.sub, paste0(path, dataset, '/', 'DEG_', i, '.tsv'), sep="\t")
    write.table(sigRes.sub[1:100,], paste0(path, dataset, '/', 'DEG_top100_', i, '.tsv'), sep="\t")
    write.table(rownames(sigRes.sub)[1:101], paste0(path, dataset, '/', 'DEG_top100_genes_', i, '.tsv'), sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
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
                      max.overlaps = 100,
                      data =sigRes.sub[1:25, ])
    jpeg(paste0(path, dataset, '/', dataset, '_volcano_plot_', i,'.png'))
    print(volcanoPlot)
    dev.off()
    print(volcanoPlot)
  
      
    #Enrichment analysis
    head(listEnrichrDbs())
    enrichmentRes <-enrichr(
      genes =rownames(sigRes.sub),
      databases =c("Allen_Brain_Atlas_10x_scRNA_2021", "KEGG_2019_Human"))
    # check entries from the two databases 
    print(enrichmentRes$Allen_Brain_Atlas_10x_scRNA_2021[1:5, 1:4])
    print(enrichmentRes$KEGG_2019_Human[1:5, 1:4])
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
    jpeg(paste0(path, dataset, '/', dataset, '_GO_terms_', i,'.png'))
    plot(arrangeGrob(grobs =enrichPlots, layout_matrix =rbind(c(1, 2))))
    dev.off()
    plot(arrangeGrob(grobs =enrichPlots, layout_matrix =rbind(c(1, 2))))
  }

}


