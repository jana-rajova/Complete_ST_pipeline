library(jpeg)
library(png)
library(raster)
library(akima)
library(future.apply)
library(stringr)
library(scales)
plan(multisession, workers=20)

Image.locations <- '/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Images_rev1/'
interpol.mask.location <- '/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/data/ST_interpolation_masks/'

alg <- 'seurat'
cluster.folder <- paste0('/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/results/Batch_corrections/', alg, '/')
dataset <- 'SN'
cluster.folder <- paste0(cluster.folder, dataset, '/')
plot.folder <- paste0(cluster.folder, 'plt/')
dir.create(paste0(plot.folder,'/png'), recursive = TRUE, showWarnings = FALSE)
cluster.file <- grep('_clusters_combined.tsv', list.files(cluster.folder), value = TRUE)
cluster.table <- read.table(paste0(cluster.folder, cluster.file), header = TRUE)
#load cluster files
sample_ <- 'ST3_D2'

cluster.colors.matrix.general <- c('#000080', '#00008b', '#000097', '#0000a3', '#0000ae', '#0000ba', '#0000c6', '#0000d1', '#0000dd', '#0000e9', '#0000f5', '#0000ff', '#0000ff', '#0006ff', '#0011ff', '#001bff', '#0025ff', '#0030ff', '#003aff', '#0044ff', '#004fff', '#0059ff', '#0063ff', '#006dff', '#0078ff', '#0082ff', '#008cff', '#0097ff', '#00a1ff', '#00abff', '#00b6ff', '#00c0ff', '#00caff', '#00d5ff', '#00dffc', '#03e9f4', '#0bf3ec', '#14fee3', '#1cffdb', '#24ffd3', '#2cffca', '#35ffc2', '#3dffba', '#45ffb1', '#4effa9', '#56ffa1', '#5eff98', '#67ff90', '#6fff88', '#77ff80', '#80ff77', '#88ff6f', '#90ff67', '#98ff5e', '#a1ff56', '#a9ff4e', '#b1ff45', '#baff3d', '#c2ff35', '#caff2c', '#d3ff24', '#dbff1c', '#e3ff14', '#ecff0b', '#f4f903', '#fcef00', '#ffe600', '#ffdc00', '#ffd300', '#ffc900', '#ffc000', '#ffb600', '#ffad00', '#ffa300', '#ff9900', '#ff9000', '#ff8600', '#ff7d00', '#ff7300', '#ff6a00', '#ff6000', '#ff5700', '#ff4d00', '#ff4400', '#ff3a00', '#ff3100', '#ff2700', '#ff1d00', '#ff1400', '#f50a00', '#e90100', '#dd0000', '#d10000', '#c60000', '#ba0000', '#ae0000', '#a30000', '#970000', '#8b0000', '#800000'

)

samples <- unique(cluster.table$sample)

if (dataset != 'separately_clustered'){
  #no of colors needed: # you can generate the colors from a cmap of your choice in "color picker.py"
  no.clusters <- length(unique(cluster.table$cluster))
  print(no.clusters)
  cluster.colors.matrix <- cluster.colors.matrix.general[round(seq(1, length(cluster.colors.matrix.general), length=no.clusters),0)]
} else {
  cluster.colors.matrix <- NULL
}
  

for (sample_ in samples){
  print(sample_)
  if (file.exists(paste0(Image.locations,'corrected_png/', sample_,'_HE.png'))){
    he.image <- readPNG(paste0(Image.locations,'corrected_png/', sample_,'_HE.png'))
  } else {
    he.image <- readJPEG(paste(Image.locations,sample_,'_HE.jpg', sep=""))
  }
  bw.image <- readJPEG(paste0(Image.locations, sample_, '_mask.jpg'))

  #cluster.dataframe.file <- paste0(cluster.folder, grep(sample_, cluster.files, value=TRUE))
  
  cluster.dataframe <- subset(cluster.table, sample == sample_)
  clusters <- cluster.dataframe$cluster
  
  if (dataset == 'separately_clustered'){
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
  
  
  
  max.x = max(coords[,1]) + 0.5
  min.x = min(coords[,1]) - 0.5
  
  max.y = max(coords[,2]) + 0.5
  min.y = min(coords[,2]) - 0.5
  
  par(pin = c(max(img[,1]) * 0.1, max(img[,2]) * 0.1))
  
  # plot(img[,1], -img[,2], cex = 0.05, pch = 19, xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE)
  
  # C. Set the percentage of dots you want to keep, nrow(img.2) should not be more
  # than 50,000 (if the plot of tissue morphology looks OK?).
  
  # limit = 0.1
  feature.multiplier = 1
  set.seed(1)
  
  img.2.size <- feature.multiplier * length(clusters)
  img.2 = img[sample(nrow(img),size=img.2.size,replace=FALSE),]
 
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
  # if (file.exists(new.set2.file)){
  #   print(paste0('Found ', new.set2.file))
  #   new.set2 <- as.integer(read.table(new.set2.file)$V1)
  #   if (length(new.set2) != img.2.size){
  #     plan(multisession, workers = 20)
  #     print(paste(length(new.set2), img.2.size))
  #     print(paste0("Calculating 'new.set2' for ", sample_, " this might take a while..."))
  #     new.set2 = future_apply(set1, 1, dist2b)
  #     dir.create(interpol.mask.location, showWarnings = FALSE, recursive = TRUE)
  #     write.table(x = new.set2, file = paste0(interpol.mask.location, sample_, '_newset2.tsv'),
  #                 row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
  #   }
  #   else {
  #     print(paste0("Loading 'new.set2' for ", sample_))
  #     new.set2 <- as.integer(read.table(new.set2.file)$V1)
  #     print(length(new.set2))
  #   }
  # } 
  # if (!file.exists(new.set2.file)) {
    # plan(multisession, workers = 20)
    # print(paste0("Calculating 'new.set2' for ", sample_, " this might take a while..."))
    # new.set2 = future_apply(set1, 1, dist2b)
    # dir.create(interpol.mask.location, showWarnings = FALSE, recursive = TRUE)
    # write.table(x = new.set2, file = new.set2.file,
                # row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
  # }
  
  # section.input = cbind(img.2[,1], img.2[,2], new.set2)
  section.input = cbind(img.2[,1], img.2[,2])
  y = r@nrows
  
  # colnames(section.input) = c(x, y, 'grid.cell')
  colnames(section.input) = c(x, y)
  
  #if there is cluster no0- push the while queue, R has no idea what is index 0
  if (min(cluster.table$cluster)==0){
    clusters <- clusters + 1
  }

  cluster.colors <- clusters
  for (i in unique(clusters)){
    # print(typeof(i))
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
  plot(x, y, cex = 2.5, ylim = c(35, 1), xlim = c(1, 33), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = alpha(cluster.colors, 0.5), pch = 16)
  
  legend(-1, 0, c(sort(unique(clusters))), col = cluster.colors.matrix[sort(unique(clusters))] , pch = 19, xpd = TRUE)
  
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
  
  
  # x = as.numeric(colnames(section.input)[1])
  # y = as.numeric(colnames(section.input)[2])
  
  x = 500
  y = round(x * (as.numeric(colnames(section.input)[2])/as.numeric(colnames(section.input)[1])), 0)
  cols.list = cluster.colors.matrix
  mycol = c()
  cluster.int = matrix(nrow = 0, ncol = 3)
  
  coords = cluster.dataframe$feature
  coords = unlist(strsplit(coords, "_"))
  coords = as.numeric(gsub('X', '', coords))
  coords = t(matrix(coords,nrow = 2))
  cluster.exp = cbind(coords,cluster.dataframe$cluster)
  if (min(cluster.exp[, 3])==0){
    cluster.exp[, 3] <- cluster.exp[, 3] + 1
  }
  
  x.1 = as.numeric(cluster.exp[,1])
  y.1 = as.numeric(cluster.exp[,2])
  cluster.ex.df = data.frame(row.names = seq(1, x*y))
  print('Calculating interpolation')
  for (cluster in sort(unique(cluster.exp[,3]))){
    z.1 = as.numeric(as.numeric(cluster.exp[,3]) == as.numeric(cluster))
    s1 =  interp(x.1, y.1, z.1, nx = x, ny = y, linear = TRUE)
    x.coord = s1$x
    y.coord = s1$y
    s1 = s1$z
    mat.1 = as.data.frame(s1)
    colnames(mat.1) = c(1:y)
    mat.1 = mat.1[,rev(colnames(mat.1)),]
    mat.1 = as.matrix(mat.1)
    col.C1 = as.numeric(mat.1)
    # set = col.C1[section.input[,3]]
    cluster.ex.df[, as.character(cluster)] = col.C1
  }
  cluster.ex.df = apply(cluster.ex.df,1, which.max)
  cols  = c()
  count = 0
  print('Mapping clusters')
  for (val in cluster.ex.df){
    if (identical(val, integer(0))){
      cols = c(cols, NA)
      count = count + 1
    } else {
    cols= c(cols, as.numeric(names(val)))
    }
  }

  a = rep(x.coord, length(y.coord))
  b = c()
  for (y in y.coord){
    b = c(b, rep(y, length(x.coord)))
  }
  colx = cols
  for (i in na.omit(sort(unique(cols)))){
    print(i)
    colx = gsub(paste0('^', i, '$'), cluster.colors.matrix[i], colx)
  }
  df = data.frame(a,b, colx)
  df.nona = na.omit(df)
  # plot(df.nona[,1], df.nona[,2], col = df.nona[,3], cex=0.1, pch=19)
  # 
  # z.1 = as.numeric(as.numeric(cluster.exp[,3]))
  # 
  # s1 =  interp(x.1, y.1, z.1, nx = x, ny = y, linear = FALSE)
  # mat.1 = s1$z
  # 
  # mat.1 = as.data.frame(mat.1)
  # colnames(mat.1) = c(1:y)
  # mat.1 = mat.1[,rev(colnames(mat.1)),]
  # mat.1 = as.matrix(mat.1)
  # 
  # col.C1 = as.numeric(mat.1)
  # set = col.C1[section.input[,3]]
  # section.C1.xyz.value = cbind(section.input[,1:3], cols)
  # # section.C1.xyz.value = cbind(section.C1.xyz.value, as.numeric(lapply(section.C1.xyz.value[,4], round, 0)))
  # section.C1.xyz.value = cbind(section.C1.xyz.value, section.C1.xyz.value[, 4])
  # for (i in na.omit(sort(unique(section.C1.xyz.value[, 4])))){
  #   print(i)
  #   section.C1.xyz.value[, 5] = gsub(paste0('^', i, '$'), cluster.colors.matrix[i], section.C1.xyz.value[, 5])
  # }
  
  # M. Finally, plot the interpolated clusters. As an optional step, this image
  # can be combined with the grey-scale image for better visualization (see
  # separate instructions).
  
  out.file = paste(sample_, '_clusters_interpolated_', alg, '.pdf', sep = '')
  pdf(paste0(plot.folder, out.file), onefile = TRUE, useDingbats = FALSE)
  
  # plot(c(1,35), c(33,1), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = 'white')
  #  rasterImage(he.image, 1, 1, 35, 33)
  #  par(new = T)
  plot(as.numeric(df.nona[,1]), -as.numeric(df.nona[,2]), cex = 0.00005, ylim = c(-max(df.nona[,2]), -min(df.nona[,2])), xlim = c(min(df.nona[,1]), max(df.nona[,1])), pch = 19, xlab='', ylab='', col = df.nona[, 3], xaxt='n', yaxt='n', axes=FALSE)

  dev.off()
  
  out.file = paste(sample_, '_clusters_interpolated_', alg, '.png', sep = '')
  png(paste0(plot.folder, '/png/', out.file))
  
  # plot(c(1,35), c(33,1), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = 'white')
  #  rasterImage(he.image, 1, 1, 35, 33)
  #  par(new = T)
  plot(as.numeric(df.nona[,1]), -as.numeric(df.nona[,2]), cex = 0.00005, ylim = c(-max(df.nona[,2]), -min(df.nona[,2])), xlim = c(min(df.nona[,1]), max(df.nona[,1])), pch = 19, xlab='', ylab='', col = df.nona[, 3], xaxt='n', yaxt='n', axes=FALSE)
  legend(-1, 0, c(sort(unique(clusters))), col = cluster.colors.matrix[sort(unique(clusters))] , pch = 19, xpd = TRUE)
  dev.off()

}
