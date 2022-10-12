library(jpeg)
library(png)
library(raster)
library(akima)
library(future.apply)
library(stringr)
library(scales)
plan(multisession)

Image.locations <- '/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Images_rev1/'
interpol.mask.location <- '/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/data/ST_interpolation_masks/'

alg <- 'seurat'
cluster.folder <- paste0('/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/results/Batch_corrections/', alg, '/')
dataset <- 'TX'
cluster.folder <- paste0(cluster.folder, dataset, '/')
plot.folder <- paste0(cluster.folder, 'plt/')
dir.create(plot.folder, recursive = TRUE, showWarnings = FALSE)
cluster.file <- grep('_clusters_combined.tsv', list.files(cluster.folder), value = TRUE)
cluster.table <- read.table(paste0(cluster.folder, cluster.file), header = TRUE)
#load cluster files
sample_ <- 'ST3_D2'

cluster.colors.matrix.general <- c('#9e0142', '#a40743', '#a90d45', '#af1346', '#b41a47', '#ba2049', '#bf264a', '#c52c4b', '#ca324d', '#d0384e', '#d53e4f', '#d8434e', '#dc484c', '#df4d4b', '#e2514a', '#e55649', '#e85b48', '#eb6046', '#ee6445', '#f16944', '#f46e44', '#f57547', '#f67b4a', '#f7824d', '#f88950', '#f98f53', '#fa9656', '#fb9c59', '#fba35c', '#fca95f', '#fdb062', '#fdb567', '#fdba6b', '#fdbf6f', '#fdc473', '#fec977', '#fece7c', '#fed380', '#fed884', '#fedd88', '#fee18d', '#fee492', '#fee898', '#feeb9d', '#feeea2', '#fff1a7', '#fff4ad', '#fff7b2', '#fffab7', '#fffdbc', '#fefebd', '#fbfdb9', '#f9fcb5', '#f6fbb1', '#f4faad', '#f1f9a9', '#eff8a5', '#ecf7a1', '#eaf69e', '#e7f59a', '#e2f499', '#dcf19a', '#d7ef9b', '#d1ec9c', '#cbea9e', '#c5e79f', '#bfe5a0', '#b9e3a1', '#b3e0a2', '#addea4', '#a6dba4', '#9fd8a4', '#98d6a4', '#91d3a4', '#8ad0a4', '#83cda5', '#7ccba5', '#75c8a5', '#6ec5a5', '#67c3a5', '#62bda7', '#5db7a9', '#57b2ac', '#52acae', '#4da6b1', '#48a0b3', '#429ab5', '#3d94b8', '#388eba', '#3389bd', '#3683bb', '#3a7db8', '#3f77b5', '#4372b2', '#486cb0', '#4c66ad', '#5160aa', '#555ba7', '#5a55a5', '#5e4fa2'
)

samples <- unique(cluster.table$sample)

if (alg != 'non_corrected'){
  #no of colors needed: # you can generate the colors from a cmap of your choice in "color picker.py"
  no.clusters <- length(unique(cluster.table$cluster))
  print(no.clusters)
  cluster.colors.matrix <- cluster.colors.matrix.general[round(seq(1, length(cluster.colors.matrix.general), length=no.clusters),0)]
} else {
  cluster.colors.matrix <- NULL
}
  

for (sample_ in samples){
  if (file.exists(paste0(Image.locations,'corrected_png/', sample_,'_HE.png'))){
    he.image <- readPNG(paste0(Image.locations,'corrected_png/', sample_,'_HE.png'))
  } else {
    he.image <- readJPEG(paste(Image.locations,sample_,'_HE.jpg', sep=""))
  }
  bw.image <- readJPEG(paste0(Image.locations, sample_, '_mask.jpg'))

  #cluster.dataframe.file <- paste0(cluster.folder, grep(sample_, cluster.files, value=TRUE))
  
  cluster.dataframe <- subset(cluster.table, sample == sample_)
  clusters <- cluster.dataframe$cluster
  
  if (is.null(cluster.colors.matrix)){
    cluster.colors.matrix <- cluster.colors.matrix.general[round(seq(2, length(cluster.colors.matrix.general), 
                                                                     length=length(unique(clusters) +1)), 0)]
  }
  
  # create image intrapolation variables
  limit = 0.65
  
  img = which(bw.image < limit , arr.ind = TRUE)
  
  test.y = attributes(bw.image)$dim[1] / 34
  test.x = attributes(bw.image)$dim[2] / 32
  
  V1 = (img[,2] / test.x) + 1
  V2 = (img[,1] / test.y) + 1
  img = cbind(V1, V2)
  
  coords = cluster.dataframe$feature
  coords = unlist(strsplit(substring(coords, 2), "_"))
  coords <- as.numeric(coords)
  coords = t(matrix(coords,nrow = 2))
  #img.2.size <- round(50000*(length(coords)/2/1000*(1.1)), 0)
  img.2.size <- limit*nrow(img)
  
  
  max.x = max(coords[,1]) + 0.5
  min.x = min(coords[,1]) - 0.5
  
  max.y = max(coords[,2]) + 0.5
  min.y = min(coords[,2]) - 0.5
  
  par(pin = c(max(img[,1]) * 0.1, max(img[,2]) * 0.1))
  
  # plot(img[,1], -img[,2], cex = 0.05, pch = 19, xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE)
  
  # C. Set the percentage of dots you want to keep, nrow(img.2) should not be more
  # than 50,000 (if the plot of tissue morphology looks OK?).
  
  limit = 0.01
  
  set.seed(1)
  
  
  #img.2 = img[sample(nrow(img),size=limit*nrow(img),replace=FALSE),]
  img.2 = img[sample(nrow(img), size=img.2.size, replace=FALSE),]
  
  print(paste0("Number of spots: ", nrow(img.2)))
  
  set.seed(NULL, normal.kind = "default")
  
  par(pin = c(max(img[,1]) * 0.1, max(img[,2]) * 0.1))
  
  # plot(img.2[,1], -img.2[,2], cex = 0.1, pch = 19, xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE)
  
  # D. Generate a “grid cell”-size for x to interpolate over.
  
  
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
  
  
  new.set2.file <- paste0(interpol.mask.location, sample_, '_newset2.tsv')
 if (file.exists(new.set2.file)){
   new.set2 <- as.integer(read.table(new.set2.file)$V1)
   if (length(new.set2) != img.2.size){
     print(paste0("Calculating 'new.set2' for ", sample_, " this might take a while..."))
     new.set2 = future_apply(set1, 1, dist2b)
     dir.create(interpol.mask.location, showWarnings = FALSE, recursive = TRUE)
     write.table(x = new.set2, file = paste0(interpol.mask.location, sample_, '_newset2.tsv'),
                 row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
   }
   else {
     print(paste0("Loading 'new.set2' for ", sample_))
   new.set2 <- as.integer(read.table(new.set2.file)$V1)
   print(length(new.set2))
   }
 } 
 if (!file.exists(new.set2.file)) {
    print(paste0("Calculating 'new.set2' for ", sample_, " this might take a while..."))
    new.set2 = future_apply(set1, 1, dist2b)
    dir.create(interpol.mask.location, showWarnings = FALSE, recursive = TRUE)
    write.table(x = new.set2, file = new.set2.file,
                row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
  }
  section.input = cbind(img.2[,1], img.2[,2], new.set2)
  
  y = r@nrows
  
  colnames(section.input) = c(x, y, 'grid.cell')
  
  #if there is cluster no0- push the while queue, R has no idea what is index 0
  if (min(cluster.table$cluster)==0){
    clusters <- clusters + 1
  }

  cluster.colors <- clusters
  for (i in unique(clusters)){
    cluster.colors = gsub(paste0('^', i, '$'), cluster.colors.matrix[i], cluster.colors)
  }
  
  # I. With this code we assign a cluster to each coordinate.
  
  clust = as.matrix(as.numeric(clusters))
  #clust = as.matrix(clust)
  rownames(clust) <- cluster.dataframe$feature
  
  coords = rownames(clust)
  coords = unlist(strsplit(coords, "_"))
  coords = as.numeric(gsub('X', '', coords))
  coords = t(matrix(coords,nrow = 2))
  clust = cbind(clust, coords)
  
  x = clust[,2]
  y = clust[,3]
  
  # J. Do the spatial cluster plotting.
  
  out.file = paste(sample_, '_clusters_original.pdf', sep = '')
  pdf(paste0(plot.folder, out.file), onefile = TRUE, useDingbats = FALSE, height = 35 * 0.3, width = 33 * 0.3)
  
  plot(c(1,35), c(33,1), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = 'white')
  rasterImage(he.image, 1, 1, 35, 33)
  par(new = T)
  plot(x, y, cex = 3, ylim = c(35, 1), xlim = c(1, 33), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = alpha(cluster.colors, 0.5), pch = 16)
  
  legend(-1, 0, c(unique(clusters)), col = cluster.colors.matrix[unique(clusters)] , pch = 19, xpd = TRUE)
  
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
  clusters.table = matrix(nrow = nrow(cluster.dataframe), ncol = max(unique(clusters)))
  rownames(clusters.table) = cluster.dataframe$feature
  
  for (i in unique(clusters)){
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
  
  for (i in unique(clusters)){
    coords = cluster.dataframe$feature
    coords = unlist(strsplit(coords, "_"))
    coords = as.numeric(gsub('X', '', coords))
    coords = t(matrix(coords,nrow = 2))
    cluster.exp = cbind(coords, clusters.table[,i])
    
    x.1 = as.numeric(cluster.exp[,1])
    y.1 = as.numeric(cluster.exp[,2])
    z.1 = as.numeric(cluster.exp[,3])
    
    s1 =  interp(x.1, y.1, z.1, nx = x, ny = y, duplicate = 'strip' )
    mat.1 = s1$z
    
    mat.1 = as.data.frame(mat.1)
    colnames(mat.1) = c(1:y)
    mat.1 = mat.1[,rev(colnames(mat.1)),]
    mat.1 = as.matrix(mat.1)
    
    col.C1 = as.numeric(mat.1)
    set = col.C1[section.input[,3]]
    section.C1.xyz.value = cbind(section.input[,1:3], set)
    
    cols = cols.list[i]
    
    breaks = seq(from = 0, to = 1, length = 100)
    rb.pal = colorRampPaletteAlpha(c(addalpha(cols, 0), addalpha(cols, 1)), 100)
    
    dat.col = rb.pal[as.numeric(cut(section.C1.xyz.value[,4],breaks = breaks))]
    
    mycol = cbind(mycol, dat.col)
    
    cluster.int.1 = cbind(section.C1.xyz.value[,1], section.C1.xyz.value[,2], section.C1.xyz.value[,4])
    
    cluster.int = rbind(cluster.int, cluster.int.1)
    
  }
  
  # M. Finally, plot the interpolated clusters. As an optional step, this image
  # can be combined with the grey-scale image for better visualization (see
  # separate instructions).
  
  out.file = paste(sample_, '_clusters_interpolated_', alg, '.pdf', sep = '')
  pdf(paste0(plot.folder, out.file), onefile = TRUE, useDingbats = FALSE, height = y * 0.075, width = x * 0.075)
  
 # plot(c(1,35), c(33,1), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = 'white')
#  rasterImage(he.image, 1, 1, 35, 33)
#  par(new = T)
  plot(cluster.int[,1], -cluster.int[,2], cex = 0.35, ylim = c(-35, -1), xlim = c(1, 33), pch = 19, col = mycol, xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE)
  
  dev.off()

  }

