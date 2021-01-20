# Step 4: Normalize your data using the "scran" approach.
# A. Start by converting your expression matrix to a "scran" object.

sce = newSCESet(countData=data.frame(exp.values))

# B. Perform a "rough" clustering to normalize your data, do not use a very small min.size.

norm.clusters = quickCluster(sce,  min.size=50)

# C. Check number of clusters and number of features in smallest cluster.

param.A <- length(unique(norm.clusters))

param.B <- floor(length(norm.clusters[norm.clusters == length(unique(norm.clusters))])/2)

# D. Compute factors for normalization, sizes should have same length as number
# of clusters and each number should be unique. The maximum value should not be
# more than 50% of the number of features in the smallest cluster. In this case
# we have 5 clusters and minimum number of features is 85. That means we could
# use sizes = c(20, 25, 30, 35, 40) since that's 5 numbers and the maximum = 40
# < 85 * 0.5.

setSize <- as.integer(seq(as.integer(param.B/param.A),param.B,length.out= param.A))
sce = computeSumFactors(sce, clusters=norm.clusters, positive = TRUE, sizes = setSize, assay="counts")

# E. Get normalized values and converts the normalized "scran" object to a matrix.

log2.norm.exp.values = normalize(sce)

log2.norm.exp.values.matrix = as.matrix(log2.norm.exp.values)

#removes erroneous normalization error producing NaNs
log2.norm.exp.values.matrix <- log2.norm.exp.values.matrix[ , colSums(is.na(log2.norm.exp.values.matrix)) == 0]


# F. Save the normalized expression matrixes for later use.

save(sce, file = paste(sample, '_sce', sep=''))

save(log2.norm.exp.values, file = paste(sample, '_log2_norm_exp_values', sep=''))

save(log2.norm.exp.values.matrix, file = paste(sample, '_log2_norm_exp_values_matrix', sep=''))

# Step 5: Find variable and correlated genes.
# A. Find varied genes based on Variance.

variance = apply(log2.norm.exp.values.matrix, 1, var)
mean.exp = rowMeans(log2.norm.exp.values.matrix)
var.m.exp = cbind(variance, mean.exp)
var.genes = var.m.exp[order(var.m.exp[,1], decreasing = T),]

# B. Plot the 25 most variable genes and their Variance and mean expression. Dashed red lines show the average value across all features.

num.genes = 25

par(las=2)
par(mfrow=c(1, 2))
par(pin = c(2, 2))

barplot(var.genes[1:num.genes,1], xlab='Variance', names.arg=rownames(var.genes)[1:num.genes], cex.names = 0.5, main = 'Genes with high Variance', horiz =T)
abline(v=mean(var.genes[,1]), col='red', lwd=2, lty=2)

barplot(var.genes[1:num.genes,2], xlab='Mean expression', names.arg=rownames(var.genes)[1:num.genes], cex.names = 0.5, main = 'Genes with high Variance', horiz =T)
abline(v=mean(var.genes[,2]), col='red', lwd=2, lty=2)

# Fing gene name
symbol = toTable(org.Hs.egSYMBOL)
annotation = toTable(org.Hs.egENSEMBL)
gene.list <- merge(symbol, annotation, by = "gene_id")
EN.id = "ENSG00000162989"
gene.list[which(gene.list$ensembl_id == EN.id), ]

# C. Select a gene of interest and find other genes that have similar expression pattern (based on correlation).

gene = 'MBP'
cor.genes = t(cor(log2.norm.exp.values.matrix[gene,], t(log2.norm.exp.values.matrix)))
mean.exp = rowMeans(log2.norm.exp.values.matrix)
cor.m.exp = cbind(cor.genes, mean.exp)
cor.genes = cor.m.exp[order(cor.m.exp[,1], decreasing = T),]
cor.genes = cor.genes[2:nrow(cor.genes),]

# B. Plot the top 25 genes with similar expression to your gene of interest. Dashed red lines show the average value across all features.

num.genes = 25

par(las=2)
par(mfrow=c(1, 2))
par(pin = c(2, 2))

barplot(cor.genes[1:num.genes,1], xlab='Correlation', names.arg=rownames(cor.genes)[1:num.genes], cex.names = 0.5, xlim = c(0:1), main = paste('Genes with high correlation to', gene), horiz =T)
abline(v=mean(cor.genes[,1]), col='red', lwd=2, lty=2)

barplot(cor.genes[1:num.genes,2], xlab='Mean expression', names.arg=rownames(cor.genes)[1:num.genes], cex.names = 0.5, main = paste('Genes with high correlation to', gene), horiz =T)
abline(v=mean(log2.norm.exp.values.matrix), col='red', lwd=2, lty=2)

# C. Run PCA with your data in order to separate regions with different expression pattern.

out.pca = prcomp(t(log2.norm.exp.values.matrix))

# D. As an alternative to PCA, run t-SNE with your data, test different values
# on perplexity and/or initial dimensions (50 is default). Use set.seed if you
# want reproducible results.

perplexity.1 = 50
perplexity.2 = 20

dims.1 = 50
dims.2 = 50

set.seed(1)

out.tsne.1 = Rtsne(as.matrix(t(log2.norm.exp.values.matrix)), dims = 2, initial_dims = dims.1, theta = 0.0, check_duplicates = FALSE, pca = TRUE, perplexity = perplexity.1, max_iter = 1000, verbose = F)

out.tsne.2 = Rtsne(as.matrix(t(log2.norm.exp.values.matrix)), dims = 2, initial_dims = dims.2, theta = 0.0, check_duplicates = FALSE, pca = TRUE, perplexity = perplexity.2, max_iter = 1000, verbose = F)

set.seed(NULL, normal.kind = "default")

# E. Visualize expression of a selected gene in PCA and t-SNE. Start by
# selecting a gene and make a color scale (for the expression heatmap).

gene = 'MBP'

exp.col = c('blue', 'yellow', 'red')

# F. Find features that express the gene.

gene.exp = as.matrix(log2.norm.exp.values.matrix[match(gene, rownames(log2.norm.exp.values.matrix)), ])


z = gene.exp[,1]

min = min(gene.exp[,1])
max = max(gene.exp[,1])

breaks = seq(from = (min - 0.00000001), to = max, length = 1001)
rb.pal = colorRampPalette(exp.col)
dat.col = rb.pal(1000)[as.numeric(cut(z, breaks = breaks))]

# G. Plot PCA with the expression.

par(mfrow=c(1, 2))
par(pin = c(2, 2))

plot(out.pca$x[,1], out.pca$x[,2], pch = 19, cex =0.75, xlab='Component 1', ylab='Component 2', cex.lab = 1, col = dat.col, main = paste('PCA with expression of', gene), xpd = TRUE)
legend_image = as.raster(matrix(rb.pal(100), ncol=1))

rasterImage(legend_image, 1.1 * max(out.pca$x[,1]), 0.5 * (min(out.pca$x[,2]) + max(out.pca$x[,2])), 1.2 * max(out.pca$x[,1]), max(out.pca$x[,2]), xpd = TRUE)

text(x=(1.40 * max(out.pca$x[,1])), y = seq(0.5 * (min(out.pca$x[,2]) + max(out.pca$x[,2])), max(out.pca$x[,2]), l=2), labels = seq(round(max(gene.exp[,1]), digits = 2),  round(min(gene.exp[,1]), digits = 2),l=2), xpd = TRUE)

plot(out.pca$x[,1], out.pca$x[,3], pch = 19, cex =0.75, xlab='Component 1', ylab='Component 3', cex.lab = 1, col = dat.col, main = paste('PCA with expression of', gene), xpd = TRUE)
legend_image = as.raster(matrix(rb.pal(100), ncol=1))

rasterImage(legend_image, 1.1 * max(out.pca$x[,1]), 0.5 * (min(out.pca$x[,3]) + max(out.pca$x[,3])), 1.2 * max(out.pca$x[,1]), max(out.pca$x[,3]), xpd = TRUE)

text(x=(1.40 * max(out.pca$x[,1])), y = seq(0.5 * (min(out.pca$x[,3]) + max(out.pca$x[,3])), max(out.pca$x[,3]), l=2), labels = seq(round(max(gene.exp[,1]), digits = 2),  round(min(gene.exp[,1]), digits = 2),l=2), xpd = TRUE)

# H. Plot t-SNE with the expression.
par(mfrow=c(1, 2))
par(pin = c(2, 2))

plot(out.tsne.1$Y[,1], out.tsne.1$Y[,2], pch = 19, cex =0.75, xlab='t-SNE 1', ylab='t-SNE 2', cex.lab = 1, col = dat.col, main = paste('t-SNE with expression of', gene, '\n', 'perplexity =', perplexity.1, '\n', 'initial dims =', dims.1), xpd = TRUE)
legend_image = as.raster(matrix(rb.pal(100), ncol=1))

rasterImage(legend_image, 1.1 * max(out.tsne.1$Y[,1]), 0.5 * (min(out.tsne.1$Y[,2]) + max(out.tsne.1$Y[,2])), 1.2 * max(out.tsne.1$Y[,1]), max(out.tsne.1$Y[,2]), xpd = TRUE)

text(x=(1.40 * max(out.tsne.1$Y[,1])), y = seq(0.5 * (min(out.tsne.1$Y[,2]) + max(out.tsne.1$Y[,2])), max(out.tsne.1$Y[,2]), l=2), labels = seq(round(max(gene.exp[,1]), digits = 2),  round(min(gene.exp[,1]), digits = 2),l=2), xpd = TRUE)

plot(out.tsne.2$Y[,1], out.tsne.2$Y[,2], pch = 19, cex =0.75, xlab='t-SNE 1', ylab='t-SNE 2', cex.lab = 1, col = dat.col, main = paste('t-SNE with expression of', gene, '\n', 'perplexity =', perplexity.2, '\n', 'initial dims =', dims.2), xpd = TRUE)
legend_image = as.raster(matrix(rb.pal(100), ncol=1))

rasterImage(legend_image, 1.1 * max(out.tsne.2$Y[,1]), 0.5 * (min(out.tsne.2$Y[,2]) + max(out.tsne.2$Y[,2])), 1.2 * max(out.tsne.2$Y[,1]), max(out.tsne.2$Y[,2]), xpd = TRUE)

text(x=(1.40 * max(out.tsne.2$Y[,1])), y = seq(0.5 * (min(out.tsne.2$Y[,2]) + max(out.tsne.2$Y[,2])), max(out.tsne.2$Y[,2]), l=2), labels = seq(round(max(gene.exp[,1]), digits = 2),  round(min(gene.exp[,1]), digits = 2),l=2), xpd = TRUE)

# I. Now, move on to show spatial expression of a selected gene in the tissue image, start by loading the HE-image.

he.image = readJPEG('CN56_C2_HE.jpg')


# J. Select a gene you want to display and a color scale (for the expression heatmap).

gene = 'MBP'

exp.col = c('blue', 'yellow', 'red')

# K. Find features that express the gene.

gene.exp = as.matrix(log2.norm.exp.values.matrix[match(gene, rownames(log2.norm.exp.values.matrix)), ])
coords = rownames(gene.exp)
coords = unlist(strsplit(coords, "_"))
coords = as.numeric(gsub('X', '', coords))
coords = t(matrix(coords,nrow = 2))
gene.exp = cbind(gene.exp, coords)

x = gene.exp[,2]
y = gene.exp[,3]
z = gene.exp[,1]

min = min(gene.exp[,1])
max = max(gene.exp[,1])

breaks = seq(from = (min - 0.00000001), to = max, length = 1001)
rb.pal = colorRampPalette(exp.col)
dat.col = rb.pal(1000)[as.numeric(cut(z, breaks = breaks))]

# L. Generate a pdf-with the expression pattern plotted on top of the tissue image.

out.file = paste(sample, '_', gene, '.pdf', sep = '')
pdf(out.file, onefile = TRUE, useDingbats = FALSE, height = 35 * 0.5, width = 33 * 0.5)

plot(c(1,35), c(33,1), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = 'white')
rasterImage(he.image, 1, 1, 35, 33)
par(new = T)
plot(x, y, cex = 3, ylim = c(35, 1), xlim = c(1, 33), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = alpha(dat.col, 0.5), pch = 16)

legend_image = as.raster(matrix(rb.pal(100), ncol=1))

rasterImage(legend_image, -1, 0, 0, 5, xpd = TRUE)

text(x=1, y = seq(0,5,l=2), labels = seq(round(max(gene.exp[,1]), digits = 2),  round(min(gene.exp[,1]), digits = 2), l=2), cex = 2, xpd = TRUE)

dev.off()


# Initial step: Load your previous data.
# A. Load the packages.

library(igraph)
library(SQUAREM)
library(CountClust)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(ReactomePA)
library(Rtsne)
library(scran)
library(jpeg)
library(raster)
library(akima)
library(vioplot)
library(devtools)
library(dynamicTreeCut)
library(ggplot2)
library(gplots)
library(edgeR)
library(cellTree)

# B. Set your working directory.

setwd("~/")

# C. Load the normalized expression matrixes.

sample = 'CN56_C2'

load(file = 'CN56_C2_exp_values')
load(file = 'CN56_C2_sce')
load(file = 'CN56_C2_log2_norm_exp_values') 
load(file = 'CN56_C2_log2_norm_exp_values_matrix')

# Step 6: Convert tissue image to binary dots and generate an grid for
# interpolation. A. As an alternative to showing expression and other
# information (clusters etc.) as features, you can register the HE-image (grey
# scale version, see separate instructions) as binary dots and interpolate
# across it. Start by loading the greyscale-image.

bw.image = readJPEG('CN56_C2_mask.jpg')

# B. Covert to binary format, play around with limit until you are satisfied
# with the binary dots you get.

limit = 0.65

img = which(bw.image < limit , arr.ind = TRUE)

test.y = attributes(bw.image)$dim[1] / 34
test.x = attributes(bw.image)$dim[2] / 32

V1 = (img[,2] / test.x) + 1
V2 = (img[,1] / test.y) + 1
img = cbind(V1, V2)

coords = colnames(exp.values)
coords = unlist(strsplit(coords, "_"))
coords = as.numeric(gsub('X', '', coords))
coords = t(matrix(coords,nrow = 2))

max.x = max(coords[,1]) + 0.5
min.x = min(coords[,1]) - 0.5

max.y = max(coords[,2]) + 0.5
min.y = min(coords[,2]) - 0.5

par(pin = c(max(img[,1]) * 0.06, max(img[,2]) * 0.06))

plot(img[,1], -img[,2], cex = 0.05, pch = 19, xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE)

# C. Set the percentage of dots you want to keep, nrow(img.2) should not be more
# than 50,000 (if the plot of tissue morphology looks OK?).

limit = 0.0135

set.seed(1)

img.2 = img[sample(nrow(img),size=limit*nrow(img),replace=FALSE),]

nrow(img.2)

set.seed(NULL, normal.kind = "default")

par(pin = c(max(img[,1]) * 0.06, max(img[,2]) * 0.06))

plot(img.2[,1], -img.2[,2], cex = 0.1, pch = 19, xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE)

# D. Generate a "grid cell"-size for x to interpolate over.

coords = colnames(exp.values)
coords = unlist(strsplit(coords, "_"))
coords = as.numeric(gsub('X', '', coords))
coords = t(matrix(coords,nrow = 2))

x = round((max(coords[,1]) - min(coords[,1]) + 1), digits = 0) * 4

# E. Assign each binary dot to a grid-cell; this step takes time (10-15 min).

res = (max.x - min.x) / x

r = raster(xmn = min.x, ymn = min.y, xmx = max.x, ymx = max.y, res = res)

r[] = 0

tab = table(cellFromXY(r, img.2))

r[as.numeric(names(tab))] = tab

pixel.centers = coordinates(r)
set1 = img.2[,1:2]
set2 = pixel.centers[,1:2]

distp1p2 = function(p1,p2) {
  dst = sqrt((p1[1]-p2[1])^2+(p1[2]-p2[2])^2)
  return(dst)
}

dist2b = function(y) which.min(apply(set2, 1, function(x) min(distp1p2(x,y))))
# dist2b = function(y) which.min(mclapply(set2, 1, function(x) min(distp1p2(x,y))))

new.set2 = apply(set1, 1, dist2b)


section.input = cbind(img.2[,1], img.2[,2], new.set2)

y = r@nrows

colnames(section.input) = c(x, y, 'grid.cell')

# F. Save the new object for later use.

save(section.input, file = paste(sample, "_section_input", sep=""))

# load(file = paste(sample, "_section_input", sep=""))
# Fing gene name
symbol = toTable(org.Hs.egSYMBOL)
annotation = toTable(org.Hs.egENSEMBL)
gene.list <- merge(symbol, annotation, by = "gene_id")
EN.id = "ENSG00000135046"
gene.list[which(gene.list$ensembl_id == EN.id), ]


# G. Plot the interpolated gene expression across the binary dots. First select
# a gene you want to display and a color scale (for the expression heatmap).

plotInterpol <- function(gene) {
  
  if (!is.na(match(gene, toupper(rownames(log2.norm.exp.values.matrix))))) {
    exp.col = c('blue', 'yellow', 'red')
    
    # H. Now, generate the interpolation.
    
    x = as.numeric(colnames(section.input)[1])
    y = as.numeric(colnames(section.input)[2])
    
    gene.exp = as.matrix(log2.norm.exp.values.matrix[match(gene, toupper(rownames(log2.norm.exp.values.matrix))), ])
    coords = rownames(gene.exp)
    coords = unlist(strsplit(coords, "_"))
    coords = as.numeric(gsub('X', '', coords))
    coords = t(matrix(coords,nrow = 2))
    gene.exp = cbind(gene.exp, coords)
    
    x.1 = gene.exp[,2]
    y.1 = gene.exp[,3]
    z.1 = gene.exp[,1]
    
    min = min(gene.exp[,1])
    max = max(gene.exp[,1])
    
    s1 =  interp(x.1, y.1, z.1, nx = x, ny = y)
    mat.1 = s1$z
    mat.1 = as.data.frame(mat.1)
    colnames(mat.1) = c(1:y)
    mat.1 = mat.1[,rev(colnames(mat.1)),]
    mat.1 = as.matrix(mat.1)
    mat.1 = as.numeric(mat.1)
    set = mat.1[section.input[,3]]
    section.value = cbind(section.input[,1:2], set)
    
    min = min(gene.exp[,1])
    max = max(gene.exp[,1])
    
    breaks = seq(from = (min - 0.00000001), to = max, length = 1001)
    rb.pal = colorRampPalette(exp.col)
    dat.col = rb.pal(1000)[as.numeric(cut(section.value[,3], breaks = breaks))]
    dat.col[is.na(dat.col)] = 'grey'
    
    # I. Finally, plot the interpolated expression. As an optional step, this image can be combined with the grey-scale image for better visualization (see separate instructions).
    
    out.file = paste(sample, '_', gene,'_interpolated.pdf', sep = '')
    pdf(out.file, onefile = TRUE, useDingbats = FALSE, height = y * 0.5, width = x * 0.5)
    
    plot(section.value[,1], -section.value[,2], cex = 1.5, pch = 19, col = alpha(dat.col, 1), xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE)
    
    dev.off()
  } 
}

geneList <- read.table("genelist.txt", header = FALSE, skip = 0, sep="\t",
                       stringsAsFactors = FALSE, fill=TRUE)


invisible(mclapply(geneList$V1,function(x) plotInterpol(x), mc.cores = detectCores()))

# Step 7: Hierarchical clustering. A. Here we are going to carry out
# hierarchical clustering and visualization based on distance derived from
# correlation between features. Here the data is already normalized so you can
# set a very small value on the minimum numbers of features in a cluster
# (min.size).

clusters = quickCluster(log2.norm.exp.values,  min.size=5)


##Added part of the code, which adds dataset exrtacted from scanorama (dimensionally reduced datasets after batch correction)
##The file is called "dataset-dimred" + the well designation (CN65C2 in this case) + "expdata.csv"
#it still has to be made into R data.frame, columns must be renamed (0-99) and first column (-1) has to be made rownames and then dropped from the data.frame content

length(unique(clusters))

# B. Alternatively, you can carry out hierarchical clustering based on Euclidean
# distances between features, using Ward's criterion to minimize the total
# variance within each cluster. Features with similar expression patterns are
# grouped together. This method is less robust to noise but more sensitive to
# subtle changes in the expression. Then, define clusters by applying a dynamic
# tree cut (Langfelder, Zhang, and Horvath 2008) to the dendrogram. This
# approach uses the shape of the branches in the dendrogram to refine the
# cluster definitions.

# my.dist = dist(t(log2.norm.exp.values.matrix))
# 
# my.tree = hclust(my.dist, method="ward.D2")
# 
# clusters = unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), verbose=0))
# 
# length(unique(clusters))

# C. Generate colors for the clusters, here we selected these colors for 12
# theoretical clusters (even thou we only have 8-9 here).

cluster.colors.matrix = c('red', 'blue', 'magenta', 'orange', 'maroon', 'Cadetblue', 'yellow', 'cyan', 'burlywood', 'black', 'purple', 'dark grey')

cluster.colors = clusters
for (i in length(unique(clusters)):1){
  cluster.colors = gsub(i, cluster.colors.matrix[i], cluster.colors)
  print(cluster.colors[i])
}

# D. Run PCA with your data.
out.pca = prcomp(t(log2.norm.exp.values.matrix))

# E. Visualize PCA with features colored according to hierarchical clusters.
# Here we visualize component 1 vs. component 2 and component 1 vs. component 3.

par(mfrow=c(1, 2))
par(pin = c(2, 2))

plot(out.pca$x[,1], out.pca$x[,2], pch = 19, cex =0.75, xlab='Component 1', ylab='Component 2', cex.lab = 1, col = cluster.colors, main = 'PCA with clusters', xpd = TRUE)

legend(1.1 * max(out.pca$x[,1]), 1.1 * max(out.pca$x[,2]), c(1:length(unique(clusters))), col = cluster.colors.matrix[1:length(unique(clusters))], pch = 19, xpd = TRUE, bty = 'n', cex = 0.75)

plot(out.pca$x[,1], out.pca$x[,3], pch = 19, cex =0.75, xlab='Component 1', ylab='Component 3', cex.lab = 1, col = cluster.colors, main = 'PCA with clusters', xpd = TRUE)

legend(1.1 * max(out.pca$x[,1]), 1.1 * max(out.pca$x[,3]), c(1:length(unique(clusters))), col = cluster.colors.matrix[1:length(unique(clusters))], pch = 19, xpd = TRUE, bty = 'n', cex = 0.75)

# F. Run t-SNE with your data, test different values on perplexity and/or initial dimensions (50 is default). Use set.seed if you want reproducible results.

perplexity.1 = 5
perplexity.2 = 20

dims.1 = 100
dims.2 = 100

set.seed(1)

out.tsne.1 = Rtsne(as.matrix(t(log2.norm.exp.values.matrix)), dims = 2, initial_dims = dims.1, theta = 0.0, check_duplicates = FALSE, pca = TRUE, perplexity = perplexity.1, max_iter = 1000, verbose = F)

out.tsne.2 = Rtsne(as.matrix(t(log2.norm.exp.values.matrix)), dims = 2, initial_dims = dims.2, theta = 0.0, check_duplicates = FALSE, pca = TRUE, perplexity = perplexity.2, max_iter = 1000, verbose = F)

set.seed(NULL, normal.kind = "default")

# G. Visualize t-SNE with features colored according to hierarchical clusters.

par(mfrow=c(1, 2))
par(pin = c(2, 2))

plot(out.tsne.1$Y[,1], out.tsne.1$Y[,2], pch = 19, cex =0.75, xlab='t-SNE 1', ylab='t-SNE 2', cex.lab = 1, col = cluster.colors, main = paste('t-SNE clusters\n', 'perplexity =', perplexity.1, '\n', 'initial dims =', dims.1), xpd = TRUE)

legend(1.1 * max(out.tsne.1$Y[,1]), 1.1 * max(out.tsne.1$Y[,2]), c(1:length(unique(clusters))), col = cluster.colors.matrix[1:length(unique(clusters))], pch = 19, xpd = TRUE, bty = 'n', cex = 0.75)

plot(out.tsne.2$Y[,1], out.tsne.2$Y[,2], pch = 19, cex =0.75, xlab='t-SNE 1', ylab='t-SNE 2', cex.lab = 1, col = cluster.colors, main = paste('t-SNE clusters\n', 'perplexity =', perplexity.2, '\n', 'initial dims =', dims.2), xpd = TRUE)

legend(1.1 * max(out.tsne.2$Y[,1]), 1.1 * max(out.tsne.2$Y[,2]), c(1:length(unique(clusters))), col = cluster.colors.matrix[1:length(unique(clusters))], pch = 19, xpd = TRUE, bty = 'n', cex = 0.75)

# H. Now, plot the spatial locations of each cluster. Start by loading your HE-image.

he.image = readJPEG('CN56_C2_HE.jpg')

# I. With this code we assign a cluster to each coordinate.

clust = as.matrix(as.numeric(clusters))
clust = as.matrix(clust)
rownames(clust) = colnames(exp.values)

coords = rownames(clust)
coords = unlist(strsplit(coords, "_"))
coords = as.numeric(gsub('X', '', coords))
coords = t(matrix(coords,nrow = 2))
clust = cbind(clust, coords)

x = clust[,2]
y = clust[,3]

# J. Do the spatial cluster plotting.

out.file = paste(sample, '_clusters_original.pdf', sep = '')
pdf(out.file, onefile = TRUE, useDingbats = FALSE, height = 35 * 0.3, width = 33 * 0.3)

plot(c(1,35), c(33,1), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = 'white')
rasterImage(he.image, 1, 1, 35, 33)
par(new = T)
plot(x, y, cex = 3, ylim = c(35, 1), xlim = c(1, 33), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = alpha(cluster.colors, 0.5), pch = 16)

legend(-1, 0, c(1:length(unique(clusters))), col = cluster.colors.matrix[1:length(unique(clusters))] , pch = 19, xpd = TRUE)

dev.off()

# K. As for the gene expression, you can also interpolate the clusters (color)
# and show the binary dots. Start by generating a color palette for transparent
# dots (between clusters).

addalpha = function(colors, alpha=1.0) {
  r = col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] = alpha*255
  r = r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

colorRampPaletteAlpha = function(colors, n=32, interpolate='linear') {
  cr = colorRampPalette(colors, interpolate=interpolate)(n)
  a = col2rgb(colors, alpha=T)[4,]
  if (interpolate=='linear') {
    l = approx(a, n=n)
  } else {
    l = spline(a, n=n)
  }
  l$y[l$y > 255] = 255 # Clamp if spline is > 255
  cr = addalpha(cr, l$y/255.0)
  return(cr)
}

# L. Now, generate the interpolation.
#original script without scanorama data
clusters.table = matrix(nrow = ncol(exp.values) , ncol = length(unique(clusters)))
rownames(clusters.table) = colnames(exp.values)

for (i in 1:length(unique(clusters))){
  clusters.temp = as.character(clusters)
  clusters.temp[clusters.temp != i]  = '0'
  clusters.temp[clusters.temp == i]  = '1'
  clusters.table[,i] = clusters.temp
}


x = as.numeric(colnames(section.input)[1])
y = as.numeric(colnames(section.input)[2])

cols.list = cluster.colors.matrix
mycol = c()
cluster.int = matrix(nrow = 0, ncol = 3)

for (i in 1:length(unique(clusters))){
  
  coords = colnames(exp.values)
  coords = unlist(strsplit(coords, "_"))
  coords = as.numeric(gsub('X', '', coords))
  coords = t(matrix(coords,nrow = 2))
  cluster.exp = cbind(coords, clusters.table[,i])
  
  x.1 = as.numeric(cluster.exp[,1])
  y.1 = as.numeric(cluster.exp[,2])
  z.1 = as.numeric(cluster.exp[,3])
  
  s1 =  interp(x.1, y.1, z.1, nx = x, ny = y)
  mat.1 = s1$z
  
  mat.1 = as.data.frame(mat.1)
  colnames(mat.1) = c(1:y)
  mat.1 = mat.1[,rev(colnames(mat.1)),]
  mat.1 = as.matrix(mat.1)
  
  col.C1 = as.numeric(mat.1)
  set = col.C1[section.input[,3]]
  section.C1.xyz.value = cbind(section.input[,1:3], set)
  
  cols = cols.list[i]
  
  breaks = seq(from = 0, to = 1, length = 10)
  rb.pal = colorRampPaletteAlpha(c(addalpha(cols, 0), addalpha(cols, 1)), 10)
  
  dat.col = rb.pal[as.numeric(cut(section.C1.xyz.value[,4],breaks = breaks))]
  
  mycol = cbind(mycol, dat.col)
  
  cluster.int.1 = cbind(section.C1.xyz.value[,1], section.C1.xyz.value[,2], section.C1.xyz.value[,4])
  
  cluster.int = rbind(cluster.int, cluster.int.1)
  
}

# M. Finally, plot the interpolated clusters. As an optional step, this image
# can be combined with the grey-scale image for better visualization (see
# separate instructions).

out.file = paste(sample, '_clusters_interpolated_original.pdf', sep = '')
pdf(out.file, onefile = TRUE, useDingbats = FALSE, height = y * 0.075, width = x * 0.075)

plot(cluster.int[,1], -cluster.int[,2], cex = 1.5, pch = 19, col = mycol, xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE)

dev.off()

##In this section is the scanorama integrated data in order to get the old interpolations and the new interpolations from the same run 
#and thus avoind extra added variability in the result
#writeout into abother direcrory with other 
setwd("./dim-exp")
for (union in c("True", "False")){
  for (dim in c(50, 100, 250, 500, 750)){
    scanorama_data = data.frame(read.csv(paste("dataset-dimred_CN56_C2_expdata_", union, dim, ".csv", sep = '')))
    colnames(scanorama_data) <- -1:(dim-1)
    rownames(scanorama_data) <-scanorama_data$"-1"
    scanorama_data = as.matrix(scanorama_data[,!(names(scanorama_data) %in% "-1")])
    
    
    clusters = quickCluster(t(scanorama_data),  min.size=5)
    
    
    length(unique(clusters))
    
    # B. Alternatively, you can carry out hierarchical clustering based on Euclidean
    # distances between features, using Ward's criterion to minimize the total
    # variance within each cluster. Features with similar expression patterns are
    # grouped together. This method is less robust to noise but more sensitive to
    # subtle changes in the expression. Then, define clusters by applying a dynamic
    # tree cut (Langfelder, Zhang, and Horvath 2008) to the dendrogram. This
    # approach uses the shape of the branches in the dendrogram to refine the
    # cluster definitions.
    
    # my.dist = dist(t(log2.norm.exp.values.matrix))
    # 
    # my.tree = hclust(my.dist, method="ward.D2")
    # 
    # clusters = unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), verbose=0))
    # 
    # length(unique(clusters))
    
    # C. Generate colors for the clusters, the first row is if you have color you want to use, 
    # the second one is a wes color palette (but that doesn't work well for too many clusters)
    # I decided to go with the standard rainbow as that gives you the best visual distinction
    # package pals seems like the best one, but is not available for R 3.3.2
    
    #  cluster.colors.matrix = c('red','violet','blueviolet', 'lightsalmon', 'mistyrose', 'olivedrab', 'maroon', 'springgreen', 'rosybrown','seagreen', 'blue', 'magenta', 'orange','sienna' ,'green', 'darkorchid', 'darkslateblue','goldenrod','coral', 'darkkhaki','darkred', 'darkseagreen', 'darkmagenta' , 'maroon', 'Cadetblue', 'yellow', 'cyan', 'burlywood', 'black', 'purple', 'dark grey', 'firebrick')
    # cluster.colors.matrix = wes_palette("Zissou", length(unique(clusters)), type = 'continuous')
    cluster.colors.matrix = rainbow(length(unique(clusters)))
    print(cluster.colors.matrix)
    cluster.colors = clusters
    for (i in length(unique(clusters)):1){
      print(i)
      cluster.colors = gsub(paste0("^", i), cluster.colors.matrix[i], cluster.colors)
    }
    print(cluster.colors)
    # D. Run PCA with your data.
    #out.pca = prcomp(t(log2.norm.exp.values.matrix))
    out.pca = prcomp(scanorama_data) 
    # E. Visualize PCA with features colored according to hierarchical clusters.
    # Here we visualize component 1 vs. component 2 and component 1 vs. component 3.
    
    par(mfrow=c(1, 2))
    par(pin = c(2, 2))
    
    plot(out.pca$x[,1], out.pca$x[,2], pch = 19, cex =0.75, xlab='Component 1', ylab='Component 2', cex.lab = 1, col = cluster.colors, main = 'PCA with clusters', xpd = TRUE)
    
    legend(1.1 * max(out.pca$x[,1]), 1.1 * max(out.pca$x[,2]), c(1:length(unique(clusters))), col = cluster.colors.matrix[1:length(unique(clusters))], pch = 19, xpd = TRUE, bty = 'n', cex = 0.75)
    
    plot(out.pca$x[,1], out.pca$x[,3], pch = 19, cex =0.75, xlab='Component 1', ylab='Component 3', cex.lab = 1, col = cluster.colors, main = 'PCA with clusters', xpd = TRUE)
    
    legend(1.1 * max(out.pca$x[,1]), 1.1 * max(out.pca$x[,3]), c(1:length(unique(clusters))), col = cluster.colors.matrix[1:length(unique(clusters))], pch = 19, xpd = TRUE, bty = 'n', cex = 0.75)
    
    # F. Run t-SNE with your data, test different values on perplexity and/or initial dimensions (50 is default). Use set.seed if you want reproducible results.
    
    perplexity.1 = 5
    perplexity.2 = 20
    
    dims.1 = 100
    dims.2 = 100
    
    set.seed(1)
    
    #out.tsne.1 = Rtsne(as.matrix(t(log2.norm.exp.values.matrix)), dims = 2, initial_dims = dims.1, theta = 0.0, check_duplicates = FALSE, pca = TRUE, perplexity = perplexity.1, max_iter = 1000, verbose = F)
    
    #out.tsne.2 = Rtsne(as.matrix(t(log2.norm.exp.values.matrix)), dims = 2, initial_dims = dims.2, theta = 0.0, check_duplicates = FALSE, pca = TRUE, perplexity = perplexity.2, max_iter = 1000, verbose = F)
    
    out.tsne.1 = Rtsne(scanorama_data, dims = 2, initial_dims = dims.1, theta = 0.0, check_duplicates = FALSE, pca = TRUE, perplexity = perplexity.1, max_iter = 1000, verbose = F)
    
    out.tsne.2 = Rtsne(scanorama_data, dims = 2, initial_dims = dims.2, theta = 0.0, check_duplicates = FALSE, pca = TRUE, perplexity = perplexity.2, max_iter = 1000, verbose = F)
    
    
    set.seed(NULL, normal.kind = "default")
    
    # G. Visualize t-SNE with features colored according to hierarchical clusters.
    
    par(mfrow=c(1, 2))
    par(pin = c(2, 2))
    
    plot(out.tsne.1$Y[,1], out.tsne.1$Y[,2], pch = 19, cex =0.75, xlab='t-SNE 1', ylab='t-SNE 2', cex.lab = 1, col = cluster.colors, main = paste('t-SNE clusters\n', 'perplexity =', perplexity.1, '\n', 'initial dims =', dims.1), xpd = TRUE)
    
    legend(1.1 * max(out.tsne.1$Y[,1]), 1.1 * max(out.tsne.1$Y[,2]), c(1:length(unique(clusters))), col = cluster.colors.matrix[1:length(unique(clusters))], pch = 19, xpd = TRUE, bty = 'n', cex = 0.75)
    
    plot(out.tsne.2$Y[,1], out.tsne.2$Y[,2], pch = 19, cex =0.75, xlab='t-SNE 1', ylab='t-SNE 2', cex.lab = 1, col = cluster.colors, main = paste('t-SNE clusters\n', 'perplexity =', perplexity.2, '\n', 'initial dims =', dims.2), xpd = TRUE)
    
    legend(1.1 * max(out.tsne.2$Y[,1]), 1.1 * max(out.tsne.2$Y[,2]), c(1:length(unique(clusters))), col = cluster.colors.matrix[1:length(unique(clusters))], pch = 19, xpd = TRUE, bty = 'n', cex = 0.75)
    
    # H. Now, plot the spatial locations of each cluster. Start by loading your HE-image.
    
    he.image = readJPEG('../CN56_C2_HE.jpg')
    
    # I. With this code we assign a cluster to each coordinate.
    
    clust = as.matrix(as.numeric(clusters))
    clust = as.matrix(clust)
    #originalscript
    #rownames(clust) = colnames(exp.values)
    rownames(clust) = colnames(t(scanorama_data)) 
    
    coords = rownames(clust)
    coords = unlist(strsplit(coords, "_"))
    coords = as.numeric(gsub('X', '', coords))
    coords = t(matrix(coords,nrow = 2))
    clust = cbind(clust, coords)
    
    x = clust[,2]
    y = clust[,3]
    
    # J. Do the spatial cluster plotting.
    
    out.file = paste(sample, '_clusters_scanorama', union, dim, '.pdf', sep = '')
    pdf(out.file, onefile = TRUE, useDingbats = FALSE, height = 35 * 0.3, width = 33 * 0.3)
    
    plot(c(1,35), c(33,1), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = 'white')
    rasterImage(he.image, 1, 1, 35, 33)
    par(new = T)
    plot(x, y, cex = 3, ylim = c(35, 1), xlim = c(1, 33), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = alpha(cluster.colors, 0.5), pch = 16)
    
    legend(-1, 0, c(1:length(unique(clusters))), col = cluster.colors.matrix[1:length(unique(clusters))] , pch = 19, xpd = TRUE)
    
    dev.off()
    
    # K. As for the gene expression, you can also interpolate the clusters (color)
    # and show the binary dots. Start by generating a color palette for transparent
    # dots (between clusters).
    
    addalpha = function(colors, alpha=1.0) {
      r = col2rgb(colors, alpha=T)
      # Apply alpha
      r[4,] = alpha*255
      r = r/255.0
      return(rgb(r[1,], r[2,], r[3,], r[4,]))
    }
    
    colorRampPaletteAlpha = function(colors, n=32, interpolate='linear') {
      cr = colorRampPalette(colors, interpolate=interpolate)(n)
      a = col2rgb(colors, alpha=T)[4,]
      if (interpolate=='linear') {
        l = approx(a, n=n)
      } else {
        l = spline(a, n=n)
      }
      l$y[l$y > 255] = 255 # Clamp if spline is > 255
      cr = addalpha(cr, l$y/255.0)
      return(cr)
    }
    
    # L. Now, generate the interpolation.
    #original script without scanorama data
    #clusters.table = matrix(nrow = ncol(exp.values) , ncol = length(unique(clusters)))
    #rownames(clusters.table) = colnames(exp.values)
    
    
    clusters.table = matrix(nrow = ncol(t(scanorama_data)) , ncol = length(unique(clusters)))
    rownames(clusters.table) = colnames(t(scanorama_data))
    
    for (i in 1:length(unique(clusters))){
      clusters.temp = as.character(clusters)
      clusters.temp[clusters.temp != i]  = '0'
      clusters.temp[clusters.temp == i]  = '1'
      clusters.table[,i] = clusters.temp
    }
    
    
    x = as.numeric(colnames(section.input)[1])
    y = as.numeric(colnames(section.input)[2])
    
    cols.list = cluster.colors.matrix
    mycol = c()
    cluster.int = matrix(nrow = 0, ncol = 3)
    
    for (i in 1:length(unique(clusters))){
      
      coords = colnames(exp.values)
      coords = unlist(strsplit(coords, "_"))
      coords = as.numeric(gsub('X', '', coords))
      coords = t(matrix(coords,nrow = 2))
      cluster.exp = cbind(coords, clusters.table[,i])
      
      x.1 = as.numeric(cluster.exp[,1])
      y.1 = as.numeric(cluster.exp[,2])
      z.1 = as.numeric(cluster.exp[,3])
      
      s1 =  interp(x.1, y.1, z.1, nx = x, ny = y)
      mat.1 = s1$z
      
      mat.1 = as.data.frame(mat.1)
      colnames(mat.1) = c(1:y)
      mat.1 = mat.1[,rev(colnames(mat.1)),]
      mat.1 = as.matrix(mat.1)
      
      col.C1 = as.numeric(mat.1)
      set = col.C1[section.input[,3]]
      section.C1.xyz.value = cbind(section.input[,1:3], set)
      
      cols = cols.list[i]
      
      breaks = seq(from = 0, to = 1, length = 10)
      rb.pal = colorRampPaletteAlpha(c(addalpha(cols, 0), addalpha(cols, 1)), 10)
      
      dat.col = rb.pal[as.numeric(cut(section.C1.xyz.value[,4],breaks = breaks))]
      
      mycol = cbind(mycol, dat.col)
      
      cluster.int.1 = cbind(section.C1.xyz.value[,1], section.C1.xyz.value[,2], section.C1.xyz.value[,4])
      
      cluster.int = rbind(cluster.int, cluster.int.1)
      
    }
    
    # M. Finally, plot the interpolated clusters. As an optional step, this image
    # can be combined with the grey-scale image for better visualization (see
    # separate instructions).
    
    out.file = paste(sample, '_clusters_interpolated_scanorama_', union, dim, '.pdf', sep = '')
    pdf(out.file, onefile = TRUE, useDingbats = FALSE, height = y * 0.075, width = x * 0.075)
    
    plot(cluster.int[,1], -cluster.int[,2], cex = 1.5, pch = 19, col = mycol, xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE)
    
    dev.off()
  }}
# Step 8: Differential expression analysis (DEA) A. For DEA we use the edgeR
# package (Robinson, McCarthy, and Smyth 2010), which uses negative binomial
# distributions to model the transcript, counts for each cluster. Start by
# setting up a design matrix specifying which features belong to each cluster.

cluster = factor(clusters)
design = model.matrix(~0 + cluster)
colnames(design)

# B. Convert your normalized data to a edgeR-object, estimate dispersion and fit
# the model. In regular bulk data, a dispersion of 0.05-0.2 is common whilst
# single cell data can have dispersion above 0.5 (due to technical noise).

y = convertTo(sce, type="edgeR")
y = estimateDisp(y, design)
fit = glmFit(y, design)
summary(y$tagwise.dispersion)

# C. Choose two clusters that you want to compare. Also set a p-value, a
# log2-fold change limit and the number of genes your want to mark in the plot
# later on (on each side, 5 = 10 genes)

clust.A = 2
clust.B = 7
log2.fold = 1
pvalue = 0.01
text.genes = 5

# D. Carry out the DEA and filter the genes so only genes that meet your criteria are kept.

contrast = numeric(ncol(design))
contrast[clust.A] = 1
contrast[clust.B] = -1
res = glmLRT(fit, contrast=contrast)
res.deg = cbind(res$table$logFC, res$table$PValue)
rownames(res.deg) = rownames(log2.norm.exp.values.matrix)
res.deg.pvalue = res.deg[res.deg[,2] <= pvalue, ]
res.deg.plus = res.deg.pvalue[res.deg.pvalue[,1] >= log2.fold, ]
res.deg.plus = res.deg.plus[order(res.deg.plus[,1], decreasing = T),]
res.deg.minus = res.deg.pvalue[res.deg.pvalue[,1] <= -log2.fold, ]
res.deg.minus = res.deg.minus[order(res.deg.minus[,1]),]
res.deg.final = rbind(res.deg.plus, res.deg.minus)
nrow(res.deg.final)

# E. Make a Volcano plot with your DEGs colored as red (genes that meet your criteria).

par(pin = c(2.5, 2.5))

plot(res.deg[,1], -log10(res.deg[,2]), pch = 19, cex = 0.5, xlim = c(-max(rbind(-res.deg.minus, res.deg.plus)),max(rbind(-res.deg.minus, res.deg.plus))), xlab = 'Log2-fold change', ylab = '-Log10(p-value)')
points(res.deg.final[,1], -log10(res.deg.final[,2]), pch = 19, cex = 0.5, col = 'red')
text(res.deg.plus[1:text.genes,1], -log10(res.deg.plus[1:text.genes,2]), labels = rownames(res.deg.plus[1:text.genes,]), cex = 0.75, pos = 4, col = 'cyan')
text(res.deg.minus[1:text.genes,1], -log10(res.deg.minus[1:text.genes,2]), labels = rownames(res.deg.minus[1:text.genes,]), cex = 0.75, pos = 4, col = 'blue')

# F. Make a heatmap showing the DEGs genes between your clusters.

ordered.exp = log2.norm.exp.values.matrix
colnames(ordered.exp) = paste(cluster, colnames(ordered.exp))
ordered.exp = ordered.exp[,sort(colnames(ordered.exp))]

ordered.colors = sort(clusters)
numbers.clusters = c(clust.A, clust.B)
ordered.colors = ordered.colors[ordered.colors %in% numbers.clusters]

for (i in 1:length(unique(clusters))){
  ordered.colors = gsub(i, cluster.colors.matrix[i], ordered.colors)
}

deg.matrix = ordered.exp[rownames(ordered.exp) %in% rownames(res.deg.final), ]
deg.matrix = deg.matrix[,1:length(ordered.colors)]

# G. Choose your color scale for the heatmap.

mycol = colorRampPalette(c('darkblue', 'cyan', 'yellow', 'red'))(256)

# H. Plot a heatmap with the marker genes sorted into each cluster.

par(pin = c(3, 2))

heatmap.2(deg.matrix, key = T , key.xlab = "Log2(norm expression)", key.ylab = "", scale = "none", density.info = "none", trace = "none", symbreaks = F, revC = T, cexRow = 0.5, cexCol = 0.5, symkey = F, symm = F, dendrogram = "none", Rowv = T, Colv = F, col = mycol, labRow = NULL, labCol = NULL, ColSideColors = ordered.colors)

# I. By running the DEA for all clusters against each other you can find marker
# genes for each cluster. To only get clear marker genes it's a good idea to
# increase the log2.fold.

log2.fold = 2
pvalue = 0.01

res.total = matrix(ncol = 2, nrow =0)
for (i in 1:length(unique(clusters))){
  clust.A = i
  for (j in 1:length(unique(clusters))){
    if (i == j){
      break
    }
    else {
      clust.B = j
      contrast = numeric(ncol(design))
      contrast[clust.A] = 1
      contrast[clust.B] = -1
      res = glmLRT(fit, contrast=contrast)
      res.deg = cbind(res$table$logFC, res$table$PValue)
      rownames(res.deg) = rownames(log2.norm.exp.values.matrix)
      res.deg.pvalue = res.deg[res.deg[,2] <= pvalue, ]
      res.deg.plus = res.deg.pvalue[res.deg.pvalue[,1] >= log2.fold, ]
      res.deg.minus = res.deg.pvalue[res.deg.pvalue[,1] <= -log2.fold, ]
      res.deg.final = rbind(res.deg.plus, res.deg.minus)
      nrow(res.deg.final)
      res.total = rbind(res.total, res.deg.final)
    }
  }
}

res.total = unique(rownames(res.total))
length(res.total)

ordered.exp = log2.norm.exp.values.matrix
colnames(ordered.exp) = paste(cluster, colnames(ordered.exp))
ordered.exp = ordered.exp[,sort(colnames(ordered.exp))]

ordered.colors = sort(clusters)

for (i in 1:length(unique(clusters))){
  ordered.colors = gsub(i, cluster.colors.matrix[i], ordered.colors)
}

deg.matrix = ordered.exp[rownames(ordered.exp) %in% res.total, ]

# J. Choose your color scale for the heatmap.

mycol = colorRampPalette(c('darkblue', 'cyan', 'yellow', 'red'))(256)

# K. Plot a heatmap of with the marker genes sorted into each cluster.

par(pin = c(3, 2))

heatmap.2(deg.matrix, key = T , key.xlab = "Log2(norm expression)", key.ylab = "", scale = "none", density.info = "none", trace = "none", symbreaks = F, revC = T, cexRow = 0.5, cexCol = 0.5, symkey = F, symm = F, dendrogram = "none", Rowv = T, Colv = F, col = mycol, labRow = NULL, labCol = NULL, ColSideColors = ordered.colors)

# L. Make a violin plot, start by selecting a gene of interest.

gene = 'MBP'

# M. Make the plot.

vio.list = ordered.exp[gene, ]

x1 = vio.list[1:length(clusters[clusters == 1])]
x2 = vio.list[(length(x1) + 1):(length(clusters[clusters == 2]) + length(x1))]
x3 = vio.list[(length(x1) + length(x2) + 1):(length(clusters[clusters == 3]) + length(x1) + length(x2))]
x4 = vio.list[(length(x1) + length(x2) + length(x3) + 1):(length(clusters[clusters == 4]) + length(x1) + length(x2) + length(x3))]
x5 = vio.list[(length(x1) + length(x2) + length(x3) + length(x4) + 1):(length(clusters[clusters == 5]) + length(x1) + length(x2) + length(x3) + length(x4))]
x6 = vio.list[(lengRplotth(x1) + length(x2) + length(x3) + length(x4) + length(x5) + 1):(length(clusters[clusters == 6]) + length(x1) + length(x2) + length(x3) + length(x4) + length(x5))]
x7 = vio.list[(length(x1) + length(x2) + length(x3) + length(x4) + length(x5) + length(x6) + 1):(length(clusters[clusters == 7]) + length(x1) + length(x2) + length(x3) + length(x4) + length(x5) + length(x6))]
x8 = vio.list[(length(x1) + length(x2) + length(x3) + length(x4) + length(x5) + length(x6) + length(x7) + 1):(length(clusters[clusters == 8]) + length(x1) + length(x2) + length(x3) + length(x4) + length(x5) + length(x6) + length(x7))]
x9 = vio.list[(length(x1) + length(x2) + length(x3) + length(x4) + length(x5) + length(x6) + length(x7) + length(x8) + 1):(length(clusters[clusters == 9]) + length(x1) + length(x2) + length(x3) + length(x4) + length(x5) + length(x6) + length(x7) + length(x8))]

plot(1,1, xlim=c(0, length(unique(cluster)) + 1), ylim=range(c(0,vio.list)), type="n", xlab='Cluster', ylab='Log2(norm expression)', axes=FALSE)

axis(side=1,at=1:length(unique(cluster)))
axis(side=2)

vioplot(x1, at=1,col=cluster.colors.matrix[1], add=TRUE)
vioplot(x2, at=2,col=cluster.colors.matrix[2], add=TRUE)
vioplot(x3, at=3,col=cluster.colors.matrix[3], add=TRUE)
vioplot(x4, at=4,col=cluster.colors.matrix[4], add=TRUE)
vioplot(x5, at=5,col=cluster.colors.matrix[5], add=TRUE)
vioplot(x6, at=6,col=cluster.colors.matrix[6], add=TRUE)
vioplot(x7, at=7,col=cluster.colors.matrix[7], add=TRUE)
vioplot(x8, at=8,col=cluster.colors.matrix[8], add=TRUE)
vioplot(x9, at=9,col=cluster.colors.matrix[9], add=TRUE)

title(main = paste('Expression of ', gene, sep = ''))

