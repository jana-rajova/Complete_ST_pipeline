library(jpeg)
library(raster)
library(akima)
library(future.apply)
library(stringr)
library(scales)
plan(multisession)

Image.locations <- '/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Images_rev1/'
interpol.mask.location <- '/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/stereoscope/res/ST_interpolation_masks/'

stereoscope.folders <- c('/media/Dropbox/MNM projects/Spatial transcriptomics project/Manuscript preparation/stereoscope_relevant_runs/nigral_sections/VM_STR_Ctx_reference_final/Lin-selection-VM-STR-Ctx-230821_ST1/75000/',
                         '/media/Dropbox/MNM projects/Spatial transcriptomics project/Manuscript preparation/stereoscope_relevant_runs/nigral_sections/VM_STR_Ctx_reference_final/Lin-selection-VM-STR-Ctx-230821_CN57/75000/')
#gene_folder <- '/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/test_structur/ST_feature_analysis/data/ST_files/ST_matrix_STARsolo_correct_PetterREFs_Multimap/stdata/'
cutoff <- -0.00008

for (stereoscope.run in stereoscope.folders){
  setwd(stereoscope.run)
  print(stereoscope.run)
  samples <- list.dirs(stereoscope.run, full.names = FALSE, recursive = FALSE)
  samples <- samples[grep(pattern = '_stdata', x = samples)]
  samples <- str_split(samples, '_stdata', simplify = TRUE)[, 1]
  
  for (sample_ in samples){
    if (file.exists(paste(Image.locations,sample_,'_HE.png', sep=""))){
      he.image.path <- readPNG(paste(Image.locations,sample_,'_HE.png', sep=""))
    } else {
      he.image.path <- readJPEG(paste(Image.locations,sample_,'_HE.jpg', sep=""))
    }
    bw.image = readJPEG(paste(Image.locations,sample_,'_mask.jpg', sep=""))
    stereoscope.file <- paste0(stereoscope.run, sample_, '_stdata/',
                               list.files(paste0(stereoscope.run, sample_, '_stdata/')))
    stereoscope.df <- read.table(stereoscope.file, sep = '\t', row.names = 1, header = TRUE, check.names = FALSE)
    min.max.df <- stereoscope.df
    min.max.df[min.max.df == cutoff] <- NA
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
        
        coords = rownames(stereoscope.df)
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
        
        set.seed(1)
        
        img.2 = img[sample(nrow(img),size=limit*nrow(img),replace=FALSE),]
        
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
          print(paste0("Loading 'new.set2' for ", sample_, "! <3"))
          new.set2 <- as.integer(read.table(new.set2.file)$V1)
        } else {
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
        
        #creating mask for mat1:
        mat <- NA
        for (cell.type in colnames(stereoscope.df)) {
          
          # H. generate the interpolation.
          x = as.numeric(colnames(section.input)[1])
          y = as.numeric(colnames(section.input)[2])
          
          cell.type.exp = as.matrix(stereoscope.df[, cell.type])
          rownames(cell.type.exp) <- rownames(stereoscope.df)
          if (!is.na(max(cell.type.exp[,1], na.rm=TRUE))) {
            coords = rownames(cell.type.exp)
            coords = as.numeric(unlist(strsplit(substring(coords, 2), "_")))
            coords = t(matrix(coords,nrow = 2))
            cell.type.exp = cbind(cell.type.exp, coords)
            
            x.1 = cell.type.exp[,2]
            y.1 = cell.type.exp[,3]
            z.1 = cell.type.exp[,1]
            
            s1 =  interp(x.1, y.1, z.1, nx = x, ny = y)
            mat.1 = s1$z
            mat.1 = as.data.frame(mat.1)
            colnames(mat.1) = c(1:y)
            mat.1 = mat.1[,rev(colnames(mat.1)),]
            mat.1 = as.matrix(mat.1)
            
            if (all(is.na(mat))){
              mat <- mat.1
            } else {
              mat <- mat + mat.1
            }
          }
        }
        mat <- mat/length(colnames(stereoscope.df))
        mat[mat < cutoff] <- NA
        mat[mat >= cutoff] <- 1
        
        
        print('Plotting distribution for:')
        for (cell.type in colnames(stereoscope.df)) {
          print(cell.type)
          cell.type.file <- gsub('[;()/ ,]', '_', cell.type)
          cell.type.file <- gsub('_+', '_', cell.type.file)
          exp.col = c('blue', 'yellow', 'red')
          
          # H. Now, generate the interpolation.
          
          x = as.numeric(colnames(section.input)[1])
          y = as.numeric(colnames(section.input)[2])
          
          cell.type.exp = as.matrix(stereoscope.df[, cell.type])
          rownames(cell.type.exp) <- rownames(stereoscope.df)
          if (!is.na(max(cell.type.exp[,1], na.rm=TRUE))) {
            coords = rownames(cell.type.exp)
            coords = as.numeric(unlist(strsplit(substring(coords, 2), "_")))
            coords = t(matrix(coords,nrow = 2))
            cell.type.exp = cbind(cell.type.exp, coords)
            
            x.1 = cell.type.exp[,2]
            y.1 = cell.type.exp[,3]
            z.1 = cell.type.exp[,1]
            
            s1 =  interp(x.1, y.1, z.1, nx = x, ny = y)
            mat.1 = s1$z
            mat.1 = as.data.frame(mat.1)
            colnames(mat.1) = c(1:y)
            mat.1 = mat.1[,rev(colnames(mat.1)),]
            mat.1 = as.matrix(mat.1)
            mat.1 <- mat.1 * mat
            mat.1[mat.1 < 0] <- 0
            #print(table(is.na(mat.1)))
            mat.1 = as.numeric(mat.1)
            set = mat.1[section.input[,3]]
            section.value = cbind(section.input[,1:2], set)
            
            min = min(min.max.df, na.rm=TRUE)
            max = max(min.max.df, na.rm=TRUE)
            
            breaks = seq(from = (min - 0.00000001), to = max, length = 1001)
            rb.pal = colorRampPalette(exp.col)
            dat.col = rb.pal(1000)[as.numeric(cut(section.value[,3], breaks = breaks))]
            dat.col[is.na(dat.col)] = 'grey'
            
            # I. Finally, plot the interpolated expression. As an optional step, this image can be combined with the grey-scale image for better visualization (see separate instructions).
            
            dir.create(paste0('visuals_interpolation/', sample_, '/', 'PDF/'), showWarnings = FALSE, recursive = TRUE)
            out.file = paste('visuals_interpolation/', sample_, '/', 'PDF/', sample_, '_', cell.type.file,'_absolute_interpolated.pdf', sep = '')
            pdf(out.file, onefile = TRUE, useDingbats = FALSE, height = y * 0.3, width = x * 0.3)
            plot(c(1,33), c(35,1), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = 'white')
            rasterImage(he.image, 1, 1, 33, 35)
            par(new = T)
            plot(-section.value[,1], section.value[,2], cex = 1.25, pch = 19, col = alpha(dat.col, 1), xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE, xlim=c(-1,-33), ylim=c(35,1))
            dev.off()
            print(out.file)
            dir.create(paste0('visuals_interpolation/', sample_, '/', 'PNG_ctrl/'), showWarnings = FALSE, recursive = TRUE)
            out.file = paste('visuals_interpolation/', sample_, '/', 'PNG_ctrl/', sample_, '_', cell.type.file,'_absolute_interpolated.png', sep = '')
            png(out.file)
            plot(c(1,33), c(35,1), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = 'white')
            rasterImage(he.image, 1, 1, 33, 35)
            par(new = T)
            plot(-section.value[,1], section.value[,2], cex = 1.25, pch = 19, col = alpha(dat.col, 1), xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE, xlim=c(-1,-33), ylim=c(35,1)) + title(paste(strwrap(cell.type, width = 55), collapse = "\n"))
            
            dev.off()
            if (!grepl('selection_dataframes', stereoscope.run, fixed = TRUE)){
              print('computing relative positions')
              print(sample_)
              
              min = min(set, na.rm=TRUE)
              max = max(set, na.rm=TRUE)
              
              breaks = seq(from = (min - 0.00000001), to = max, length = 1001)
              rb.pal = colorRampPalette(exp.col)
              dat.col = rb.pal(1000)[as.numeric(cut(section.value[,3], breaks = breaks))]
              dat.col[is.na(dat.col)] = 'grey'
              
              # I. Finally, plot the interpolated expression. As an optional step, this image can be combined with the grey-scale image for better visualization (see separate instructions).
              
              dir.create(paste0('visuals_interpolation/', sample_, '/PDF/relative/'), showWarnings = FALSE, recursive = TRUE)
              out.file = paste('visuals_interpolation/', sample_, '/PDF/relative/', sample_, '_', cell.type.file,'_relative_interpolated.pdf', sep = '')
              pdf(out.file, onefile = TRUE, useDingbats = FALSE, height = y * 0.3, width = x * 0.3)
              plot(c(1,33), c(35,1), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = 'white')
              rasterImage(he.image, 1, 1, 33, 35)
              par(new = T)
              plot(-section.value[,1], section.value[,2], cex = 1.25, pch = 19, col = alpha(dat.col, 1), xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE, xlim=c(-1,-33), ylim=c(35,1))
              dev.off()
              print(out.file)
              dir.create(paste0('visuals_interpolation/', sample_, '/PNG_ctrl/relative/'), showWarnings = FALSE, recursive = TRUE)
              out.file = paste('visuals_interpolation/', sample_, '/PNG_ctrl/relative/', sample_, '_', cell.type.file,'_relative_interpolated.png', sep = '')
              png(out.file)
              plot(c(1,33), c(35,1), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col = 'white')
              rasterImage(he.image, 1, 1, 33, 35)
              par(new = T)
              plot(-section.value[,1], section.value[,2], cex = 1.25, pch = 19, col = alpha(dat.col, 1), xlab='', ylab='', xaxt='n', yaxt='n', axes=FALSE, xlim=c(-1,-33), ylim=c(35,1)) + title(paste(strwrap(cell.type, width = 55), collapse = "\n"))
              
              dev.off()
            }
          }
        }
      }
}
