if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
# remotes::install_github("mojaveazure/seurat-disk")

library(Seurat)
# library(SeuratDisk)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize= 54428800000)

#dataset is the name of the folder and will be added as the output directory and as a signifier to the output graphs and files
resolution = 1
dataset <- "TX-excluded"
st.file <- "/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/data/ST_files/TX_excluded.txt"
stdata.folder <- "/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/data/st_data_pre_processed/stdata_tsv/"
output.dir <- paste0("/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/results/Batch_corrections/seurat/", dataset, "/")

timestamp <- strsplit(as.character(Sys.time()), ' ')[[1]][1]


# Read in the sample names to be processed
samples.table <- read.table(st.file, sep="\n", header=FALSE)
samples <- as.character(samples.table$V1)
print(samples)

# Create a list with Seurat objects
ST.list <- list()
ST.feature.list <- list()
for (i in samples){
  print(i)
  table_st <- read.table(file=paste0(stdata.folder, i,"_stdata.tsv"), sep="\t", header=TRUE)
  rownames(table_st) <- table_st$feature
  table_st$feature <- NULL
  table_st <- as.data.frame(t(table_st))
  ST.feature.list[[i]] <- colnames(table_st)
  colnames(table_st) <- paste0(i, '-', colnames(table_st))
  
  ST.list[[i]] <-  CreateSeuratObject(counts=table_st, project=as.character(i))
  ST.list[[i]]@meta.data$well <- as.character(i)
  ST.list[[i]]@meta.data$feature <- ST.feature.list[[i]]
  print(ncol(table_st))
  
}
# Normalize the samples with SCT transform

ST.list <- ST.list[samples]
for (i in 1:length(ST.list)) {
  ST.list[[i]] <- SCTransform(ST.list[[i]], verbose = FALSE)
}

#Select features for downstream integration
ST.features <- SelectIntegrationFeatures(object.list = ST.list, nfeatures = 3000)
ST.list <- PrepSCTIntegration(object.list = ST.list, anchor.features = ST.features, 
                              verbose = FALSE)

# Find integration anchors and integrate samples (takes a long time)
ST.anchors <- FindIntegrationAnchors(object.list = ST.list, normalization.method = "SCT", 
                                     anchor.features = ST.features, verbose = FALSE)
ST.integrated <- IntegrateData(anchorset = ST.anchors, normalization.method = "SCT", 
                               verbose = FALSE)
#Dimenstionally reduce and visualize
ST.integrated <- RunPCA(ST.integrated, verbose = FALSE)

ST.integrated <- RunUMAP(ST.integrated, dims = 1:11)
DimPlot(ST.integrated, group.by = "orig.ident")

plots <- DimPlot(ST.integrated, group.by = "orig.ident")


plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
                                                                     override.aes = list(size = 3)))

DimPlot(ST.integrated, group.by = "orig.ident", split.by = "orig.ident", ncol = 3)
# DefaultAssay(ST.integrated) <- "RNA"

#Find variable features and plot
ST.integrated <- FindVariableFeatures(ST.integrated, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(ST.integrated), 10)
plot1 <- VariableFeaturePlot(ST.integrated)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
ST.integrated <- RunPCA(ST.integrated, features = VariableFeatures(object = ST.integrated))
ElbowPlot(ST.integrated, ndims = 51)
pca_number <- ST.integrated@reductions[["pca"]]@stdev[1:50] - ST.integrated@reductions[["pca"]]@stdev[2:51]
pca_number <- max(which(pca_number > 0.1))
#Jackstraw - more computationally expensive
# ST.integrated <- JackStraw(ST.integrated, num.replicate = 100)
# ST.integrated <- ScoreJackStraw(ST.integrated, dims = 1:12)
# JackStrawPlot(ST.integrated, dims = 40)




ST.integrated <- FindNeighbors(ST.integrated, dims = 1:pca_number, k.param = 20, prune.SNN = 1/15, n.trees = 100)

# ST.integrated <- FindNeighbors(ST.integrated, dims = 1:10)

ST.integrated <- FindClusters(ST.integrated, algorithm = 4, resolution = 1)
DimPlot(ST.integrated, reduction = "umap")
DimPlot(ST.integrated, reduction = "umap", split.by = 'seurat_clusters', ncol=3)
#ST.integrated <- RunUMAP(ST.integrated, dims = 1:12)


### If you want to explore different Leiden resolutions, use the following line.
# It will create folders with the datasets and resolutions as a name and put the results in the folder
# for (resolution in c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)){ 

# For our purposes, resolution 1 gives good results without excessive fragmentatioN

  
print(output.dir)
dir.create(paste0(output.dir, 'plt/'), showWarnings = FALSE, recursive = TRUE)
#ST.integrated <- FindClusters(ST.integrated, resolution = resolution)
ST.integrated <- RunUMAP(ST.integrated, dims = 1:pca_number)

jpeg(file=paste0(output.dir, "plt/UMAP_well_seurat_", dataset, ".jpg"))
DimPlot(ST.integrated, group.by = "orig.ident")
dev.off()


DimPlot(ST.integrated, reduction = "umap")

jpeg(file=paste0(output.dir, "plt/UMAP_clusters_seurat_", dataset, ".jpg"))
DimPlot(ST.integrated, reduction = "umap")
dev.off()

jpeg(file=paste0(output.dir, "plt/UMAP_genes_seurat_", dataset, ".jpg"))
# Visualize the expression on genes of interest
FeaturePlot(ST.integrated, features = c("TH", "PENK", "PDYN", "SLC6A3", "MBP", "CCK"))
# Export the cluster data into a csv file and save the Seurat object
dev.off()

jpeg(file=paste0(output.dir, "plt/UMAP_by-well_seurat_", dataset, ".jpg"))
# Visualize the expression on genes of interest
DimPlot(ST.integrated, group.by = "well")
# Export the cluster data into a csv file and save the Seurat object
dev.off()
  
UMAP_embedding <- as.data.frame(Embeddings(ST.integrated@reductions$umap))
df <- data.frame("feature"=ST.integrated@meta.data$feature,
                 "sample"=ST.integrated@meta.data$well,
                 "cluster"=as.numeric(ST.integrated@meta.data[["seurat_clusters"]]),
                 "umap1" = UMAP_embedding$UMAP_1,
                 "umap2" = UMAP_embedding$UMAP_2)

file.name <- paste0(output.dir, dataset, "_seurat_clusters_combined.tsv")
write.table(df, file = file.name, row.names=FALSE, sep="\t")  




