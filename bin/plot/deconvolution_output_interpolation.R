library(jpeg)
library(raster)
library(akima)
library(future.apply)
library(stringr)
library(scales)


Image.locations <- '/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Images_rev1/'
interpol.mask.location <- '/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/data/ST_interpolation_masks/'

deconvolution.folder <- c('/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/results/cell2loc_res/L5_CTX_M_STR_CNS_Tax4_selection_1000_astro-merge_upd/TX/cell2location_map/')

singleR_max = FALSE

features <- 'all'
deconv.file <- paste0(deconvolution.folder, 'celltype_content_', features, '_features.tsv')
deconv.file <- read.table(deconv.file, sep='\t', header=TRUE)
rownames(deconv.file) <- deconv.file$X
deconv.file$X <- NULL
output = paste0(deconvolution.folder, 'celltype_interpolation/')

if (singleR_max){
  non.cts.cols = c('region', 'sample', 'cluster', 'feature')
  singleR.cts = as.data.frame(read.table('/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/results/SingleR/hESC_TX_sc_SingleR_analysisL5_CTX_M_STR_CNS_description_selection_1000_astmerged_Tax4.tsv', sep='\t'))
  singleR.cts = sort(table(singleR.cts$Celltype_SingleR))
# the following ones were picked to show their relative proportion in the figures in relation to each other
  singleR.cts = c('Non-glutamatergic neuroblasts', 'Vascular leptomeningeal cells; mouse', 'Subventricular zone radial glia-like cells', 'Dopaminergic neurons; mouse', 'Astrocytes; mouse')
  singleR.cts = gsub('-', '.', singleR.cts)
  singleR.cts = gsub(' ', '.', singleR.cts)
  singleR.cts = gsub(';','.', singleR.cts)
  deconv.file = deconv.file[c(non.cts.cols, singleR.cts)]
  output = paste0(output, 'figure_selection_max/')
}

samples <- unique(deconv.file$sample)
samples = c('ST3_D2')
for (sample_ in samples){
  bw.image = readJPEG(paste(Image.locations,sample_,'_mask.jpg', sep=""))
  sample.deconv.df <- subset(deconv.file, sample == sample_)
  print(sample_)
  # B. Covert to binary format, play around with limit until you are satisfied
  # with the binary dots you get.
  
  limit = 0.65
  
  img = which(bw.image < limit , arr.ind = TRUE)
  
  test.y = attributes(bw.image)$dim[1] / 34
  test.x = attributes(bw.image)$dim[2] / 32
  
  V1 = (img[,2] / test.x) + 1
  V2 = (img[,1] / test.y) + 1
  img = cbind(V1, V2)
  
  coords = sample.deconv.df$feature
  coords = unlist(strsplit(substring(coords, 2), "_"))
  coords <- as.numeric(coords)
  coords = t(matrix(coords,nrow = 2))
  
  
  max.x = max(coords[,1]) + 0.5
  min.x = min(coords[,1]) - 0.5
  
  max.y = max(coords[,2]) + 0.5
  min.y = min(coords[,2]) - 0.5
  
  par(pin = c(max(img[,1]) * 0.1, max(img[,2]) * 0.1))
  
  # plot(img[,1], -img[,2], cex = 0.05, pch = 19, xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE)
  
  # C. Set the percentage of dots you want to keep, nrow(img.2) should not be more
  # than 50,000 (if the plot of tissue morphology looks OK?).
  
  
  limit = 0.2
  multiplier = 150
  
  set.seed(1)
  
  img.2.size <- as.integer(length(coords[,1]) * multiplier)
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
  if (file.exists(new.set2.file)){
    print(paste0('Found ', new.set2.file))
    new.set2 <- as.integer(read.table(new.set2.file)$V1)
    if (length(new.set2) != img.2.size){
      plan(multisession, workers = 20)
      print(paste(length(new.set2), img.2.size))
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
    plan(multisession, workers = 20)
    print(paste0("Calculating 'new.set2' for ", sample_, " this might take a while..."))
    new.set2 = future_apply(set1, 1, dist2b)
    dir.create(interpol.mask.location, showWarnings = FALSE, recursive = TRUE)
    write.table(x = new.set2, file = new.set2.file,
                row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
  }
  
  section.input = cbind(img.2[,1], img.2[,2], new.set2)
  
  y = r@nrows
  
  colnames(section.input) = c(x, y, 'grid.cell')
  
  # F. Save the new object for later use.
  
  # save(section.input, file = paste(sample_, "_section_input", sep=""))
  
  # G. Plot the interpolated gene expression across the binary dots. First select
  # a gene you want to display and a color scale (for the expression heatmap).
  
  # gene = 'MBP'
  print('Plotting distribution for:')
  cell.types <- colnames(sample.deconv.df)
  cell.types = cell.types[ - which(cell.types %in% c('sample', 'feature', 'region', 'cluster'))]
  abs_max <- max(deconv.file[, cell.types])
  
  for (cell.type in cell.types) {
    print(cell.type)
    cell.type.file <- gsub('[;()/ ,]', '_', cell.type)
    cell.type.file <- gsub('_+', '_', cell.type.file)
    exp.col = c('blue', 'yellow', 'red')
    
    # H. Now, generate the interpolation.
    
    x = as.numeric(colnames(section.input)[1])
    y = as.numeric(colnames(section.input)[2])
    
    cell.type.exp = as.matrix(sample.deconv.df[, cell.type])
    rownames(cell.type.exp) <- sample.deconv.df$feature

    
    if (!is.na((max(cell.type.exp[,1])))) {
      coords = rownames(cell.type.exp)
      coords = as.numeric(unlist(strsplit(substring(coords, 2), "_")))
      coords = t(matrix(coords,nrow = 2))
      cell.type.exp = cbind(cell.type.exp, coords)
      
      x.1 = as.numeric(cell.type.exp[,2])
      y.1 = as.numeric(cell.type.exp[,3])
      z.1 = as.numeric(cell.type.exp[,1])
      
      min = min(na.omit(cell.type.exp[,1]))
      max = max(cell.type.exp[,1])
      
      s1 =  interp(x.1, y.1, z.1, nx = x, ny = y)
      mat.1 = s1$z
      if (features == 'TX'){
        mat.1[mat.1 < 0] <- NA
      }
      mat.1 = as.data.frame(mat.1)
      colnames(mat.1) = c(1:y)
      mat.1 = mat.1[,rev(colnames(mat.1)),]
      mat.1 = as.matrix(mat.1)
      mat.1 = as.numeric(mat.1)
      set = mat.1[section.input[,3]]
      section.value = cbind(section.input[,1:2], set)
      
      breaks = seq(from = (min - 0.00000001), to = max, length = 1001)
      rb.pal = colorRampPalette(exp.col)
      dat.col = rb.pal(1000)[as.numeric(cut(section.value[,3], breaks = breaks))]
      dat.col[is.na(dat.col)] = 'grey'
      
      # I. Finally, plot the interpolated expression. As an optional step, this image can be combined with the grey-scale image for better visualization (see separate instructions).
      
      dir.create(paste0(output, features, '_features/', sample_, '/', 'PDF/relative/'), showWarnings = FALSE, recursive = TRUE)
      out.file = paste(output, features, '_features/', sample_, '/', 'PDF/relative/', sample_, '_', cell.type.file,'_relative_interpolated.pdf', sep = '')
      pdf(out.file, onefile = TRUE, useDingbats = FALSE, height = y * 0.3, width = x * 0.3)
      
      plot(section.value[,1], -section.value[,2], cex = 1.5, pch = 19, col = alpha(dat.col, 1), xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE)
      
      dev.off()
      
      

      dir.create(paste0(output, features, '_features/', sample_, '/', 'PNG_ctrl/relative/'), showWarnings = FALSE, recursive = TRUE)
      out.file = paste(output, features, '_features/', sample_, '/', 'PNG_ctrl/relative/', sample_, '_', cell.type.file,'_relative_interpolated.png', sep = '')
      png(out.file)
      
      plot(section.value[,1], -section.value[,2], cex = 0.2, pch = 19, col = alpha(dat.col, 1), xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE) + title(paste(strwrap(cell.type, width = 55), collapse = "\n"))
      
      dev.off()
      
      # do the same but replace relative min max with an universal one
      
      print('Values for absolute scaling')
      min = min(na.omit(sample.deconv.df[,cell.types]))
      print(min)
      max = abs_max
      print(max)
      
      s1 =  interp(x.1, y.1, z.1, nx = x, ny = y)
      mat.1 = s1$z
      if (features == 'TX'){
        mat.1[mat.1 < (min+abs(min/1000))] <- NA
      }
      mat.1 = as.data.frame(mat.1)
      colnames(mat.1) = c(1:y)
      mat.1 = mat.1[,rev(colnames(mat.1)),]
      mat.1 = as.matrix(mat.1)
      mat.1 = as.numeric(mat.1)
      set = mat.1[section.input[,3]]
      section.value = cbind(section.input[,1:2], set)
      
      breaks = seq(from = (min - 0.00000001), to = max, length = 1001)
      rb.pal = colorRampPalette(exp.col)
      dat.col = rb.pal(1000)[as.numeric(cut(section.value[,3], breaks = breaks))]
      dat.col[is.na(dat.col)] = 'grey'

      
      # I. Finally, plot the interpolated expression. As an optional step, this image can be combined with the grey-scale image for better visualization (see separate instructions).
      
      cbar.perc = 0.03
      
      dir.create(paste0(output, features, '_features/', sample_, '/', 'PDF/absolute/'), showWarnings = FALSE, recursive = TRUE)
      out.file = paste(output, features, '_features/', sample_, '/', 'PDF/absolute/', sample_, '_', cell.type.file,'_absolute_interpolated.pdf', sep = '')
      pdf(out.file, onefile = TRUE, useDingbats = FALSE, height = y * 0.3, width = x * (0.3 * (1 + cbar.perc + 0.1)))
      
      par(fig=c(0,0.8,0,1))
      plot(section.value[,1], -section.value[,2], cex = 1.5, pch = 19, col = alpha(dat.col, 1), xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE)
      par(fig=c(0.8,1,0,1), new = T, cex.axis=2.5)
      image(0.1, seq(0, max, len=100), matrix(1:100,nrow=1), col=rb.pal(100), axes=FALSE, xlab="", ylab="")
      axis(side=2, at=seq(0, 1, len=21), las=1)
      dev.off()
      
      dir.create(paste0(output, features, '_features/', sample_, '/', 'PNG_ctrl/absolute/'), showWarnings = FALSE, recursive = TRUE)
      out.file = paste(output, features, '_features/', sample_, '/', 'PNG_ctrl/absolute/', sample_, '_', cell.type.file,'_absolute_interpolated.png', sep = '')
      png(out.file)
      
      par(fig=c(0,0.8,0,1))
      plot(section.value[,1], -section.value[,2], cex = 0.2, pch = 19, col = alpha(dat.col, 1), xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE) + title(paste(strwrap(cell.type, width = 55), collapse = "\n"))
      par(fig=c(0.8,1,0,1), new = T)
      image(0.1, seq(0, max, len=100), matrix(1:100,nrow=1), col=rb.pal(100), axes=FALSE, xlab="", ylab="")
      axis(side=2, at=seq(0, 1, len=21), las=1)
      
      dev.off()
    }
  }

  # plot clusters
  # print('cluster')
  # exp.col = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
  # 
  # # H. Now, generate the interpolation.
  # 
  # x = as.numeric(colnames(section.input)[1])
  # y = as.numeric(colnames(section.input)[2])
  # 
  # cluster.df = read.table(paste0('/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/results/Seurat_clustered_wells_STARsolo_corrected_MultiMap_PettersREFs/TX/1/', sample_,  "_seurat_clusters_1_TX.tsv"), sep='\t', header = TRUE)
  # 
  # if (!is.na((max(clusters.exp[,1])))) {
  #   coords = cluster.df[,'feature']
  #   coords = as.numeric(unlist(strsplit(substring(coords, 2), "_")))
  #   coords = t(matrix(coords,nrow = 2))
  #   clusters.exp = cbind(cluster.df[, 'cluster'], coords)
  #   
  #   x.1 = clusters.exp[,2]
  #   y.1 = clusters.exp[,3]
  #   z.1 = clusters.exp[,1]
  #   
  #   # do the same but replace relative min max with an universal one
  #   
  #   min = min(cluster.df[, 'cluster'])
  #   max = max(cluster.df[, 'cluster'])
  #   
  #   s1 =  interp(x.1, y.1, z.1, nx = x, ny = y)
  #   mat.1 = s1$z
  #   mat.1 = as.data.frame(mat.1)
  #   colnames(mat.1) = c(1:y)
  #   mat.1 = mat.1[,rev(colnames(mat.1)),]
  #   mat.1 = as.matrix(mat.1)
  #   mat.1 = as.numeric(mat.1)
  #   set = mat.1[section.input[,3]]
  #   section.value = cbind(section.input[,1:2], set)
  # 
  #   
  #   breaks = seq(from = (min - 0.00000001), to = max, length = 1001)
  #   rb.pal = colorRampPalette(exp.col)
  #   dat.col = rb.pal(1000)[as.numeric(cut(section.value[,3], breaks = breaks))]
  #   dat.col[is.na(dat.col)] = 'grey'
  #   
  #   # I. Finally, plot the interpolated expression. As an optional step, this image can be combined with the grey-scale image for better visualization (see separate instructions).
  #   
  #   dir.create(paste0('visuals_interpolation/', sample_, '/', 'PDF/absolute/'), showWarnings = FALSE, recursive = TRUE)
  #   out.file = paste('visuals_interpolation/', sample_, '/', 'PDF/absolute/', sample_, '_', cell.type.file,'_absolute_interpolated.pdf', sep = '')
  #   pdf(out.file, onefile = TRUE, useDingbats = FALSE, height = y * 0.3, width = x * 0.3)
  #   
  #   plot(section.value[,1], -section.value[,2], cex = 1.25, pch = 19, col = alpha(dat.col, 1), xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE)
  #   
  #   dev.off()
  #   
  #   dir.create(paste0('visuals_interpolation/', sample_, '/', 'PNG_ctrl/absolute/'), showWarnings = FALSE, recursive = TRUE)
  #   out.file = paste('visuals_interpolation/', sample_, '/', 'PNG_ctrl/absolute/', sample_, '_', cell.type.file,'_absolute_interpolated.png', sep = '')
  #   png(out.file)
  #   
  #   plot(section.value[,1], -section.value[,2], cex = 0.2, pch = 19, col = alpha(dat.col, 1), xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE) + title(paste(strwrap(cell.type, width = 55), collapse = "\n"))
  #   
  #   dev.off()
    }
 