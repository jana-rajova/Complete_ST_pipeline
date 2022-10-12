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



datasets <- c('CN56_ST3')
methods <- c('seurat_DEG_analysis')
path <- '/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/results/cluster_gene_expression/'
marker.file <- '/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/data/PanglaoDB_markers_27_Mar_2020.tsv'

# test setup variables
dataset <- 'TX_1'
method <- 'seurat_DEG_analysis'
gene <- 'TH'
for (dataset in datasets){
  for (method in methods){
    output.path <- paste0(path, '/', method, '/', dataset, '/')
    countTable <-read.table(paste0(output.path, 'merged_cluster_expression.tsv'), sep='\t', header=TRUE)
    output.path <- paste0(output.path, '/results/')
    dir.create(output.path, showWarnings = FALSE)
    print(dataset)
    print(method)
    countTable$cluster <- as.factor(countTable$cluster)
    by_cluster <- countTable %>% group_by(cluster)
    
    #select top TH cluster - this used to be to determine dopaminergic neuron rich cluster
    # we changed to to the cluster with thi highest human content
    gene.top.cluster <- by_cluster %>% summarize(TH=mean(TH/feature_counts))
    hist.TH <- as.data.frame(gene.top.cluster)
    hist.TH <- hist.TH[order(hist.TH$TH, decreasing = TRUE),]
    
    gene.top.cluster <- gene.top.cluster[order(gene.top.cluster$TH, decreasing = TRUE),]
    gene.top.cluster[gene] <- scale(gene.top.cluster[gene])
    
    jpeg(paste0(output.path, dataset, '_top_TH_cluster_plot.png'))
    barplot(height=hist.TH$TH, names=hist.TH$cluster, 
            col=heat.colors(length(hist.TH$TH)),
            xlab="cluster", 
            ylab="TH count",
            las=2, cex.names=.8
    )
    dev.off()

    write.table(gene.top.cluster, paste0(output.path, 'TH_top_clusters_', method, '.tsv'), sep="\t", quote = FALSE)
    
    
    #select striatal cluster(s)
    str.genes <-c('PENK', 'ADORA2A', 'PPP1R1B')
    STR_genes_countTable <- countTable[, c('cluster', 'feature_counts', str.genes)]
    for (gene in str.genes){ 
      STR_genes_countTable[gene]  <- STR_genes_countTable[gene]/STR_genes_countTable$feature_counts
      STR_genes_countTable[gene] <- scale(STR_genes_countTable[gene])
    }
    STR_genes_countTable$str.score <- rowMeans(STR_genes_countTable[, str.genes])

    
    STR_genes_countTable_grouped <- STR_genes_countTable %>% group_by(cluster) %>% summarise_at(c("str.score", str.genes), mean, na.rm = TRUE)
    
    write.table(STR_genes_countTable_grouped, paste0(output.path, 'str_top_score_clusters_', method, '.tsv'), sep="\t", quote = FALSE, row.names = FALSE)
    
    jpeg(paste0(output.path, dataset, '_top_STR_clusters_plot.png'))
    barplot(height=STR_genes_countTable_grouped$str.score, 
            names=STR_genes_countTable_grouped$cluster, 
            col=heat.colors(length(STR_genes_countTable_grouped$str.score)),
            xlab="cluster", 
            ylab="striatal score (PENK/ADORA2A/DARPP-32 based)",
            las=2, cex.names=.8
    )
    dev.off()
    
    #MBP cluster is there for control
    MBP.top.cluster <- by_cluster %>% summarise(MBP=mean(MBP/feature_counts))
    hist.MBP <- as.data.frame(MBP.top.cluster)
    hist.MBP <- hist.MBP[order(hist.MBP$MBP, decreasing = TRUE),]
    jpeg(paste0(output.path, dataset, '_top_MBP_cluster_plot_', i,'.png'))
    barplot(height=hist.MBP$MBP, names=hist.MBP$cluster, 
                       col=heat.colors(length(hist.MBP$MBP)),
                       xlab="cluster", 
                       ylab="MBP count",
            las=2, cex.names=.8
    )
    
    dev.off()
    write.table(gene.top.cluster, paste0(output.path, 'MBP_top_clusters_', i, '.tsv'), sep="\t", quote = FALSE)
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
    class(countTable) <- 'numeric'
    countTable <- round(countTable*10, digits=0)
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
    sigRes.sub.up <-rownames(subset(sigRes, FDR <0.01&logFC.ave> 1))
    sigRes.sub.down <-rownames(subset(sigRes, FDR <0.01&logFC.ave< 1))
    sigRes.sub$log.FDR <- sigRes.sub$logFC.ave / -log10(sigRes.sub$FDR)
    write.table(sigRes.sub, paste0(output.path, 'DEG_', i, '.tsv'), sep="\t")
    write.table(sigRes.sub[1:100,], paste0(output.path, 'DEG_top100_', i, '.tsv'), sep="\t")
    write.table(rownames(sigRes.sub)[1:101], paste0(output.path, 'DEG_top100_genes_', i, '.tsv'), sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
    head(sigRes.sub, 10)
    nrow(sigRes.sub)
    
    # create a heatmap

    top.100.genes <-rownames(sigRes.sub[1:100,]) 
    gene.expression.cluster.sum <- by_cluster[, c(top.100.genes, 'cluster')]
    gene.expression.cluster.sum <- as.data.frame(gene.expression.cluster.sum %>% summarise_all(funs(sum)))
    row.names(gene.expression.cluster.sum) <- as.character(gene.expression.cluster.sum$cluster)
    gene.expression.cluster.sum$cluster <- NULL
    gene.expression.cluster.sum <- gene.expression.cluster.sum[, top.100.genes]
    gene.expression.cluster.sum <- scale(gene.expression.cluster.sum)
    #gene.expression.cluster.sum <- sapply(gene.expression.cluster.sum, as.numeric)
    gene.expression.cluster.sum.t <- as.data.frame(t(gene.expression.cluster.sum))
   
   # colnames(gene.expression.cluster.sum.t) <- rownames(gene.expression.cluster.sum)
   #rownames(gene.expression.cluster.sum.t) <- colnames(gene.expression.cluster.sum)
    cluster <- as.character(gene.top.cluster[1])
    gene.expression.cluster.sum.t <- gene.expression.cluster.sum.t[order(gene.expression.cluster.sum.t[as.numeric(cluster)], decreasing = TRUE), ]

    pdf(file=paste0(output.path, "pheatmap_top100genes.pdf"), width = 16, height = 30 ,onefile = T)
    pheatmap(gene.expression.cluster.sum.t, cellheight = 10)
    dev.off()
    
    marker.genes <- read.table(marker.file, header = TRUE, sep = '\t')
    marker.genes <- marker.genes %>% select(official.gene.symbol, cell.type)
    marker.genes <- as.data.frame(marker.genes)
    sigRes.sub.up <- as.data.frame(sigRes.sub.up)
    colnames(sigRes.sub.up) <- "official.gene.symbol"
    joined <- dplyr::left_join(sigRes.sub.up, marker.genes, by = "official.gene.symbol")
    
    cell.types <- joined$cell.type
    cell.types.counts <- as.data.frame(table(cell.types))
    cell.types.counts <- cell.types.counts[order(cell.types.counts$Freq, decreasing = TRUE),]
    write.table(cell.types.counts, paste0(output.path, dataset, '_cell_types_counts_', i, '.tsv'), sep="\t")
    
    ggplot(cell.types.counts[1:10,], aes(x=reorder(cell.types, Freq), y=Freq)) + 
      geom_bar(stat="identity", aes(fill = cell.types)) +
      coord_flip() +
      theme(
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.position = "none"
      ) +  
      scale_y_continuous(position = "right")
    ggsave(path = paste0(output.path), filename = paste0(dataset, '_cell_types_counts_top10_', i,'.png'), device='png', dpi=800,  bg = "transparent",  width = 10, height = 10)
    
    
    ggplot(cell.types.counts[1:25,], aes(x=reorder(cell.types, Freq), y=Freq)) + 
      geom_bar(stat="identity", aes(fill = cell.types))+
      coord_flip() +
      theme(
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.position = "none"
      ) + 
      scale_y_continuous(position = "right")
    ggsave(path = paste0(output.path), filename = paste0(dataset, '_cell_types_counts_top25_', i,'.png'), device='png', dpi=800,  bg = "transparent",  width = 10, height = 10)

    
    #Visualize results
    ggplot(sigRes,
                         aes(x =logFC.ave,
                             y =-log10(FDR),
                             color =ifelse(FDR <0.05&abs(logFC.ave) >1, "darkred", "grey"))) +
      geom_point() +
      xlab(expression("Fold Change, Log"[2]*"")) +
      ylab(expression("Adjusted P value, Log"[10]*"")) +
      geom_vline(xintercept =c(-1, 1), linetype ="dotted", size =1) +
      geom_hline(yintercept =-log10(0.05), linetype ="dotted", size =1) +
      theme_minimal() +
      theme(legend.position = "none", aspect.ratio=1) + 
      scale_colour_manual(values =c("darkred", "grey", "steelblue")) +
      geom_text_repel(aes(x =logFC.ave, y =-log10(FDR),
                          label = rownames(sigRes.sub)[1:25],
                          size =2, 
                          color ="steelblue"),
                      max.overlaps = 100,
                      data =sigRes.sub[1:25, ])
    #jpeg(paste0(path, dataset, '/', dataset, '_volcano_plot_', i,'.png'))
    ggsave(path = paste0(output.path), filename = paste0(dataset, '_volcano_plot_', i,'.png'), device='png',  width = 6, height = 6, dpi=800)
  
      
    #Enrichment analysis
    head(listEnrichrDbs())
    enrichmentRes <-enrichr(
      genes =rownames(sigRes.sub),
      databases =c("Allen_Brain_Atlas_10x_scRNA_2021", "KEGG_2019_Human"))
    # check entries from the two databases 
    print(enrichmentRes$Allen_Brain_Atlas_10x_scRNA_2021[1:5, 1:4])
    print(enrichmentRes$KEGG_2019_Human[1:5, 1:4])
    write.table(enrichmentRes$KEGG_2019_Human, paste0(output.path, 'KEGG_2019_human_GO_', i, '.tsv'), sep="\t")
    write.table(enrichmentRes$Allen_Brain_Atlas_10x_scRNA_2021, paste0(output.path, 'Allen_Brain_atlas_GO_', i, '.tsv'), sep="\t")
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
    ggsave(path = paste0(output.path), filename = paste0(dataset, '_GO_terms_', i,'.png'), device='png', dpi=800, width = 10, height = 8)
  }

}

