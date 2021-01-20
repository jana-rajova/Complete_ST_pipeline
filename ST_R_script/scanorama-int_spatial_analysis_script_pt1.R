library("jpeg")
library("rjson")
library("raster")
library("raster")
library("scales")
library("gdata")
library("tools")
library("roxygen2")
install.packages("data.table")
library("data.table")
# library("staligner")
source("https://bioconductor.org/biocLite.R")
biocLite("org.Rn.eg.db")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("org.Rn.eg.db")

# Setting variable names
# We now need to set some variables for our script, these tell R what our files are called, and where they are located.

setwd("/home/rstudio/")
source("convert.ensembl.ids.R")
sample_ <-"CN53_C1"
Image.locations <- "./"
spot.image.name <- "CN53_C1_cy.jpg"
spot.data.name <- "spot_data-selection-CN53_C1_cy.tsv"

# Loading data into R
# Now we will load data into R to use for the spot correction.


spot.image <- readJPEG(paste(Image.locations, spot.image.name, sep = ""))
spots <- read.table(spot.data.name, header = T)
stdata.file <- read.delim(paste(sample_,"_stdata.tsv",sep=""),sep="\t", header = T, row.names = 1)

# This is what our stdata file looks like currently:
head(stdata.file[1:3], n = 6L)


out.file = paste("Corrected_spots_", file_path_sans_ext(spot.image.name), ".pdf", sep = "")
pdf(out.file, onefile=TRUE,useDingbats=FALSE)

mar.default = c(6,4,4,2) + 0.1
par(mar = mar.default + c(-2, -2, -2, 1), pty = "s")

plot(c(1,35), c(33,1), ylab = '', xlab = '', xaxt='n', yaxt='n', bty = 'n', col = "white")
rasterImage(spot.image, 1, 1, 35, 33)
par(new=T)
plot(spots$new_x, spots$new_y, cex = 1, ylim=c(35, 1), xlim=c(1, 33), ylab = '', xlab = '', xaxt='n',   yaxt='n', bty = 'n', col = 'red')
dev.off()

rt.coords.x = round(spots$new_x, digits = 2)
rt.coords.y = round(spots$new_y, digits = 2)
selected.barcodes = cbind(rt.coords.x, rt.coords.y)
rownames(selected.barcodes) = paste(spots$x,spots$y, sep="x")
selected.barcodes <- selected.barcodes[!duplicated(rownames(selected.barcodes)),]
selected.barcodes <- selected.barcodes[!duplicated(selected.barcodes[,c(1,2)]),] # Removes duplicated positions, added MH&JR
selected.barcodes = as.matrix(selected.barcodes)
save(selected.barcodes, file = paste(file_path_sans_ext(spot.image.name),   "_selected_adjusted_spots", sep=""))

# We start by subsetting our stdata file with the array spots that are under our tissue.

inside.stdata.file <- merge(selected.barcodes,stdata.file,by="row.names",all.x=FALSE, all.y=FALSE)
rownames(inside.stdata.file) <- paste(inside.stdata.file$rt.coords.x,inside.stdata.file$rt.coords.y,sep="_")
inside.stdata.file[,1:3] <- NULL

# Currently, our data has genes as columns, and array spots as rows. 
# We would like to change this so that our genes are rows and the array spots are columns. 
# We do this with the transpose function t().

inside.stdata.file <- t(inside.stdata.file)

# Now our data looks like this:

head(inside.stdata.file[,1:6], n = 6L)

# Optional: Convert ensembl gene ids to gene names
# If it is easier to work with gene names rather than gene ids, it is possible to convert between them. 
# To do so, we first have to remove version numbers (if any) from the gene ids. For example, ENSG00000000003.14 will 
# become ENSG00000000003. We then convert to gene names, making sure to use the correct organism (in this case human).

row.names(inside.stdata.file) <- gsub("\\.[0-9]+","",row.names(inside.stdata.file), perl = TRUE)
inside.stdata.file <- convert.ensembl.ids(inside.stdata.file, org = "humanizedrat")

# This is what out stdata file looks like after conversion:
head(inside.stdata.file[,1:4], n = 6L)


# Writing data back to tables
# Now that we have adjusted the spots to their actual position, we can write the data back out to a table. 
# It is written in two forms, one with gene names as columns, and one with gene names as rows.
library("staligner")
write.coor.matrix(inside.stdata.file, sample_)
write.trans.matrix(inside.stdata.file, sample_)

# QC checking
# As a final step, we check our data to see the number of unique transcripts per feature and the number of unique genes per feature.

hist(colSums(inside.stdata.file), breaks = 50, col = "grey", xlab = "Number of unique transcripts", ylab = "Number of features", main = "The number of unique transcripts per feature")


hist(colSums(inside.stdata.file>0), breaks = 100, col = "forest green", xlab = "Number of unique genes", ylab = "Number of features", main = "The number of unique genes per feature")

# We also output a small text file that details the top 10 most expressed genes, and their expression levels, the mean number of transcripts per feature, the mean number of genes per feature and the number of features covered.

write("Expression levels for the top ten genes under the tissue","QC_data.txt")
write("", "QC_data.txt", append=TRUE)
write(names(head(sort(rowSums(inside.stdata.file), decreasing = T), n=10)), "QC_data.txt",ncolumns=10,append=TRUE)
write(head(sort(rowSums(inside.stdata.file), decreasing = T), n=10),"QC_data.txt",ncolumns=10,append=TRUE)
head(sort(rowSums(inside.stdata.file), decreasing = T), n=10)


write("", "QC_data.txt", append=TRUE)
write(paste("Mean number of transcripts per feature under tissue", round(mean(colSums(inside.stdata.file)), digits = 2), sep = " : "),"QC_data.txt",append=TRUE)
paste("Mean number of transcripts per feature under tissue", round(mean(colSums(inside.stdata.file)), digits = 2), sep = " : ")


write("", "QC_data.txt", append=TRUE)
write(paste("Mean number of genes per feature under tissue", round(mean(colSums(inside.stdata.file>0)), digits = 2), sep = " : "), "QC_data.txt", append=TRUE)
paste("Mean number of genes per feature under tissue", round(mean(colSums(inside.stdata.file>0)), digits = 2), sep = " : ")


write("", "QC_data.txt", append=TRUE)
write(paste("Number of features covered", dim(inside.stdata.file)[2], sep = " : "), "QC_data.txt", append = TRUE)
paste("Number of features covered", dim(inside.stdata.file)[2], sep = " : ")



# Step 1: Install and load packages.
# A. If you do not already have the required packages installed now is a good
# time to install them. The first ones are installed through "Bioconductor", the
# second ones through "Cran" and the last through "github".

# source("https://bioconductor.org/biocLite.R")
# biocLite("scran")
# biocLite("org.Mm.eg.db")
# biocLite("org.Hs.eg.db")
# biocLite("org.Rn.eg.db")
# biocLite("DOSE")
# biocLite("ReactomePA")
# biocLite("edgeR")
# biocLite("cellTree")
# 
# install.packages('igraph')
# install.packages('raster')
# install.packages('jpeg')
# install.packages('akima')
# install.packages('Rtsne')
# install.packages('SQUAREM')
# install.packages('vioplot')
# install.packages('devtools')

library(devtools)
# install_github("kkdey/Countclust")

# B. When the packages have been installed they have to be loaded. You have to
# load the packages every time you open R.

library(igraph)
library(SQUAREM)
library(CountClust)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(org.Rn.eg.db)
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

# Step 2: Load your data and carry out quality control. A. First set your
# working directory, this is there you have all input files and where all output
# files will be saved.

setwd("/home/rstudio/")

# B. Next, name your sample and load the expression matrix.

sample = "CN53_C1"

exp.values = read.table("CN53_C1_stdata_aligned_counts_IDs.txt", header = TRUE, row.names = 1)

# C. Replace "." with "_" in the gene names because the dot causes a bug later on.

rownames(exp.values) = gsub('\\.', '\\_', rownames(exp.values))

# D. Check how many genes and features you have in your expression matrix.

dim(exp.values)

# E. Plot histograms showing the number genes and unique transcripts with mean
# (red) and SD (black) marked in each plot. These values will help you to
# determine if a sample is good or shitty! For mouse brain, 4000-6000 genes per
# feature are usually normal. For clinical samples, the number is usually lower
# but be observant if you have a very low number (< 1000-2000 genes per
# feature).
out.file = paste("QC-plots_", sample, ".pdf", sep = "")
pdf(out.file, onefile=TRUE,useDingbats=FALSE)

par(mfrow=c(1, 2))
par(pin = c(2, 2))


hist.genes = hist(colSums((exp.values) > 0), col = 'light blue', breaks = 100, xlab = 'Number of genes', ylab = 'Number of features', main = 'Genes per feature')
legend(0,-max(hist.genes$counts)/2, c(paste('mean = ', round(mean(colSums((exp.values) > 0)))), paste('SD = ', round(sd(colSums((exp.values) > 0))))), lty=c(2, 2), lwd=c(2, 2), col=c('red', 'black'), cex=0.75, xpd = TRUE, bty = 'n')

abline(v=mean(colSums((exp.values) > 0)), col='red', lwd=2, lty=2)
abline(v=(mean(colSums((exp.values) > 0)) - sd(colSums((exp.values) > 0))), col='black', lwd=2, lty=2)
abline(v=(mean(colSums((exp.values) > 0)) + sd(colSums((exp.values) > 0))), col='black', lwd=2, lty=2)


hist.trans = hist(colSums(exp.values), col = 'light green', breaks = 100, xlab = 'Number of transcripts', ylab = 'Number of features', main = 'Transcripts per feature')
legend(0, -max(hist.trans$counts)/2, c(paste('mean = ', round(mean(colSums(exp.values)))), paste('SD = ', round(sd(colSums(exp.values))))), lty=c(2, 2), lwd=c(2, 2), col=c('red', 'black'), cex=0.75, xpd = TRUE, bty = 'n')

abline(v=mean(colSums(exp.values)), col='red', lwd=2, lty=2)
abline(v=(mean(colSums(exp.values)) - sd(colSums(exp.values))), col='black', lwd=2, lty=2)
abline(v=(mean(colSums(exp.values)) + sd(colSums(exp.values))), col='black', lwd=2, lty=2)

dev.off()


# Step 3: Filtering your data. A. Start by removing features with less than a
# desired number of genes, in this example we use 1000 genes as a minimum but
# you might want to lower this number depending on sample quality.

# min.genes = 800
# 
# exp.values.cut = cbind(colSums(exp.values !=0) >= min.genes)
# exp.values = exp.values[,exp.values.cut]
# exp.values = exp.values[rowSums(exp.values) > 0,]
# 
# # B. Check how many genes and features you have left in your expression matrix.
# 
# dim(exp.values)

# C. Alternatively, you can remove features with few transcripts (or in combination with step B).

min.trans = 300

exp.values.cut = cbind(colSums(exp.values) >= min.trans)
exp.values = exp.values[,exp.values.cut]
exp.values = exp.values[rowSums(exp.values) > 0,]

dim(exp.values)

# D. Remove lowly expressed genes; this can be done in two main ways. Remove
# based on minimum average expression, in this example we use 0.2.

# min.mean.exp = 0.1

# exp.values = exp.values[rowMeans(exp.values) > min.mean.exp,]

# dim(exp.values)

# E. Alternatively, you can remove genes based on minimum expression in minimum
# number of features. In this example we remove genes that have less than 2
# count in less than 10 features.

min.exp = 2
min.features = 4

exp.values = exp.values[rowSums(exp.values >= min.exp) >= min.features,]

dim(exp.values)

# F. Save the filtered expression matrix for later use.

save(exp.values, file = paste(sample, '_exp_values', sep=''))
write.csv(exp.values, 'CN53_C1_expdata.csv')
