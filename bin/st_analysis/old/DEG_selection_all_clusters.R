if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install("DESeq2")
install.packages("dplyr")
install.packages("pheatmap")
install.packages("ggrepel")
install.packages("RColorBrewer")
BiocManager::install("limma")
BiocManager::install("edgeR")
install.packages("enrichR")
install.packages("gridExtra")
install.packages("stringr")
install.packages("statmod")

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
library(statmod)

datasets <- c("TX")
methods <- c('seurat')
path <- "../results/cluster_gene_expression/"
marker.file <- "../data/PanglaoDB_markers_27_Mar_2020.tsv"

# gene <- 'TH'
for (dataset in datasets){
  for (i in methods){
    #load datasets and get clusters
    countTable <-read.table(paste0(path, dataset, '/', dataset, '_', i, '_cluster_expression.tsv'), sep='\t', header=TRUE)
    print(dataset)
    print(i)
    countTable$cluster <- as.factor(countTable$cluster)
    clusters <- unique(countTable$cluster)
    
    # create IDs (well - cluster) and a connected sample table
    rownames(countTable) <- paste0(countTable$well, '-', countTable$cluster)
    sampleTable <- data.frame('cluster' = countTable$cluster, 'well' = countTable$well)
    rownames(sampleTable) <- rownames(countTable)
    countTable <- t(countTable)
    countTable <- countTable[-c(1:4),]
    mode(countTable) <- 'numeric'
    
    # create a basis for contrast matrix 
    for (cluster in clusters){
      print(cluster)
      x=''
      for (j in clusters[clusters != cluster]){
        x <- c(x, paste0('cluster', cluster,'-cluster',j))
      }
      x <- x[-1]
      x
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
      sigRes.sub <-subset(sigRes, FDR <0.1&abs(logFC.ave)> 1)
      sigRes.sub.up <-rownames(subset(sigRes, FDR <0.1&logFC.ave> 1))
      sigRes.sub.down <-rownames(subset(sigRes, FDR <0.1&logFC.ave< 1))
      print(length(sigRes.sub.up))
      sigRes.sub$log.FDR <- sigRes.sub$logFC.ave / -log10(sigRes.sub$FDR)
      write.table(sigRes.sub, paste0(path, dataset, '/', 'DEG_', i, '_cluster-', cluster, '.tsv'), sep="\t")
      write.table(sigRes.sub[1:100,], paste0(path, dataset, '/', 'DEG_top100_', i, '_cluster-', cluster, '.tsv'), sep="\t")
      write.table(rownames(sigRes.sub)[1:101], paste0(path, dataset, '/', 'DEG_top100_genes_', i, '_cluster-', cluster, '.tsv'), sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
      head(sigRes.sub, 10)
      nrow(sigRes.sub)
      
      marker.genes <- read.table(marker.file, header = TRUE, sep = '\t')
      marker.genes <- marker.genes %>% select(official.gene.symbol, cell.type)
      marker.genes <- as.data.frame(marker.genes)

      
      for (reg.part in list(sigRes.sub.up, sigRes.sub.down)){
        if (reg.part == sigRes.sub.up){pat <- 'up'
        } else if (reg.part == sigRes.sub.down) {
          pat <- 'down'}
        print('WORKING ON MARKERS AND ENRICHR')
        print(cluster)
        print(pat)
        print(reg.part)
        
        reg.part <- as.data.frame(reg.part)
        colnames(reg.part) <- "official.gene.symbol"
        joined <- dplyr::left_join(reg.part, marker.genes, by = "official.gene.symbol")
        joined <- joined[complete.cases(joined),]
        
        if (nrow(joined) >=1){
          print(paste0("Cell types marked by upregulated genes: ", as.character(nrow(joined))))
          cell.types <- joined$cell.type
          cell.types.counts <- as.data.frame(table(cell.types))
          cell.types.counts <- cell.types.counts[order(cell.types.counts$Freq, decreasing = TRUE),]
          print(head(joined))
          write.table(cell.types.counts, paste0(path, dataset, '/', dataset, '_cell_types_counts_', i, '_cluster-', cluster, '.tsv'), sep="\t")
          
          if (nrow(reg.part) >= 10){dim <- 10
          } else {
          dim <- as.numeric(nrow(reg.part))}
            ggplot(cell.types.counts[1:dim,], aes(x=reorder(cell.types, Freq), y=Freq)) + 
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
          ggsave(path = paste0(path, dataset, '/'), filename = paste0(dataset, '_cell_types_counts_top10_', i, '_cluster-', cluster,'_', pat, '.png'), device='png', dpi=800,  bg = "transparent")
          
          if (length(reg.part) >= 25){dim <- 25
          } else {
          dim <- length(reg.part)}
          ggplot(cell.types.counts[1:dim,], aes(x=reorder(cell.types, Freq), y=Freq)) + 
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
          ggsave(path = paste0(path, dataset, '/'), filename = paste0(dataset, '_cell_types_counts_top25_', i, '_cluster-', cluster,'_', pat,'.png'), device='png', dpi=800,  bg = "transparent")
          
          #Enrichment analysis
          print("ENRICHR ANALYSIS")
          print(cluster)
          print(pat)
          head(listEnrichrDbs())
          enrichmentRes <-enrichr(
            genes = unlist(reg.part, use.names = FALSE),
            databases =c("Allen_Brain_Atlas_10x_scRNA_2021", "KEGG_2019_Human"))
          # check entries from the two databases 
          print(enrichmentRes$Allen_Brain_Atlas_10x_scRNA_2021[1:5, 1:4])
          print(enrichmentRes$KEGG_2019_Human[1:5, 1:4])
          write.table(enrichmentRes$KEGG_2019_Human, paste0(path, dataset, '/', 'KEGG_2019_human_GO_', i, '_cluster-', cluster,'_', pat, '.tsv'), sep="\t")
          write.table(enrichmentRes$Allen_Brain_Atlas_10x_scRNA_2021, paste0(path, dataset, '/', 'Allen_Brain_atlas_GO_', i, '_cluster-', cluster,'_', pat, '.tsv'), sep="\t")
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
          ggsave(path = paste0(path, dataset, '/'), filename = paste0(dataset, '_GO_terms_', i, '_cluster-', cluster,'_', pat,'.png'), device='png', dpi=800)
      
          }
      }
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
      ggsave(path = paste0(path, dataset, '/'), filename = paste0(dataset, '_volcano_plot_', i, '_cluster-', cluster,'.png'), device='png',  width = 6, height = 6, dpi=800)
    
        

  }
}
}


