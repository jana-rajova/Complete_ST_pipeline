library(jpeg)
library(png)
library(raster)
library(akima)
library(future.apply)
library(stringr)
library(scales)
plan(multisession)

Image.locations <- '/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Images_rev1/'
interpol.mask.location <- '/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/ST_feature_analysis/data/ST_interpolation_masks/'

gene_folder <- '/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/ST_feature_analysis/data/ST_files/ST_matrix_STARsolo_correct_PetterREFs_Multimap/stdata/'


samples <- c('ST3_E1', 'ST1_C1', 'ST3_E2', 'ST1_D2', 'CN56_D2', 'CN56_E1', 'CN56_E2', 'ST3_C2', 'ST3_D1', 'ST3_D2', 'ST1_C2', 'CN57_E1', 'CN57_E2')
for (sample_ in samples) {
  #samples <-  c('ST3_E1')
  if (file.exists(paste(Image.locations,sample_,'_HE.png', sep=""))){
  he.image <- readPNG(paste(Image.locations,sample_,'_HE.png', sep=""))
  } else {
  he.image <- readJPEG(paste(Image.locations,sample_,'_HE.jpg', sep=""))
  }
  bw.image <- readJPEG(paste(Image.locations,sample_,'_mask.jpg', sep=""))
  st.file <- paste0(gene_folder, sample_, '_stdata.tsv')
  st.df <- t(read.table(st.file, sep = '\t', row.names = 1, header = TRUE, check.names = FALSE))
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
  
  coords = rownames(st.df)
  coords = unlist(strsplit(substring(coords, 2), "_"))
  coords <- as.numeric(coords)
  coords = t(matrix(coords,nrow = 2))
  img.2.size <- round(50000*(length(coords)/2/1000*(1.1)), 0)
  
  max.x = max(coords[,1]) + 0.5
  min.x = min(coords[,1]) - 0.5
  
  max.y = max(coords[,2]) + 0.5
  min.y = min(coords[,2]) - 0.5
  
  
  par(pin = c(max(img[,1]) * 0.1, max(img[,2]) * 0.1))
  
  # plot(img[,1], -img[,2], cex = 0.05, pch = 19, xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE)
  
  # C. Set the percentage of dots you want to keep, nrow(img.2) should not be more
  # than 50,000 (if the plot of tissue morphology looks OK?).
  
  limit = 0.2
  
  set.seed(1)
  
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
    print(paste0("Loading 'new.set2' for ", sample_))
    new.set2 <- as.integer(read.table(new.set2.file)$V1)
    print(length(new.set2))
  } 
  if (!file.exists(new.set2.file) | length(new.set2) != img.2.size) {
    print(paste0("Calculating 'new.set2' for ", sample_, " this might take a while..."))
    new.set2 = future_apply(set1, 1, dist2b)
    dir.create(interpol.mask.location, showWarnings = FALSE, recursive = TRUE)
    write.table(x = new.set2, file = paste0(interpol.mask.location, sample_, '_newset2.tsv'),
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
  genes <- read.table('/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/ST_feature_analysis/data/geneList.txt', sep = '\n')$V1
  print('Plotting distribution for:')
  for (gene in genes) {
    if (gene %in% colnames(st.df)){
      print(gene)
      exp.col = c('blue', 'yellow', 'red')
      
      # H. Now, generate the interpolation.
      
      x = as.numeric(colnames(section.input)[1])
      y = as.numeric(colnames(section.input)[2])
      
      gene.exp <- as.matrix(st.df[, gene])
      rownames(gene.exp) <- rownames(st.df)
      if (!is.na((max(gene.exp[,1])))) {
        coords = rownames(gene.exp)
        coords = as.numeric(unlist(strsplit(substring(coords, 2), "_")))
        coords = t(matrix(coords,nrow = 2))
        gene.exp = cbind(gene.exp, coords)
        
        x.1 = gene.exp[,2]
        y.1 = gene.exp[,3]
        z.1 = gene.exp[,1]
        
        min = min(gene.exp[,1])
        max = max(gene.exp[,1])
        print(c(min, max))
        
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
        setwd(gene_folder)
        
        # I. Finally, plot the interpolated expression. As an optional step, this image can be combined with the grey-scale image for better visualization (see separate instructions).
        
        dir.create(paste0('../gene_expression/', sample_, '/', 'PDF/'), showWarnings = FALSE, recursive = TRUE)
        out.file = paste('../gene_expression/', sample_, '/', 'PDF/', sample_, '_', gene,'_relative_interpolated.pdf', sep = '')
        pdf(out.file, onefile = TRUE, useDingbats = FALSE, height = y * 0.3, width = x * 0.3)
        
        # plot(c(1,33), c(35,1), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = 'white')
        # rasterImage(he.image, 1, 1, 33, 35)
        # par(new = T)
        plot(-section.value[,1], section.value[,2], cex = 1.5, pch = 19, col = alpha(dat.col, 1), xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE, xlim=c(-1,-33), ylim=c(35,1))
        
        dev.off()
        
        
        dir.create(paste0('../gene_expression/', sample_, '/', 'PNG_ctrl/'), showWarnings = FALSE, recursive = TRUE)
        out.file = paste('../gene_expression/', sample_, '/', 'PNG_ctrl/', sample_, '_', gene,'_relative_interpolated.png', sep = '')
        png(out.file)
        
        plot(c(1,33), c(35,1), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = 'white')
        rasterImage(he.image, 1, 1, 33, 35)
        par(new = T)
        plot(-section.value[,1], section.value[,2], cex = 0.35, pch = 19, col = alpha(dat.col, 1), xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE, xlim=c(-1,-33), ylim=c(35,1)) + title(paste(strwrap(gene, width = 55), collapse = "\n"))
        
        dev.off()
      }
    } else {
      print(paste(gene, 'not found :('))
    }
  }   
}




