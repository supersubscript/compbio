setwd("/home/henrik/compbio/thesis/code/")
# Do the track analysis
library(RColorBrewer)
library(reshape)
library(scales) #imports alpha
library(stringr)
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicRanges")
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 3

# Put the quantified files in the right order
plant.no = 2 # 2, 4, 13, 15, 18 
quant.files = list.files(paste0("/home/henrik/compbio/thesis/data/msb_plants/plant", plant.no, "/results"), full.names = TRUE)
numbers = sapply(quant.files, function(x) strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])]) 
numbers = unname(sapply(numbers, function(x) as.integer(gsub("t|\\.txt", "", x))))
quant.files = quant.files[order(numbers)]
library(plyr)
# source("https://bioconductor.org/biocLite.R")
# biocLite("limma")
# biocLite("affy")
# library(limma)
# biocLite('preprocessCore')
#load package
library(preprocessCore)
library(affy)

# Read in data 
data = lapply(quant.files, function(x) read.table(x, header=T, sep = "\t")) 
all.data = lapply(2:ncol(data[[1]]), function(x) t(plyr::ldply(unname(sapply(data, "[", x)), rbind)))
names(all.data) = colnames(data[[1]][2:ncol(data[[1]])])

##############################################
### Get plant expression densities
##############################################
# parameter = "x"
# parameter = "y"
# parameter = "z"
# parameter = "Boa.volume"
parameter = "Mean.cell.intensity"
variable = unname(sapply(data, "[", parameter))
variable = t(plyr::ldply(variable, rbind))
normalised.variable = normalize.quantiles(variable)
par(mfrow=c(1,2))
plotDensity(variable, lwd = lw.s)
plotDensity(normalised.variable, lwd = lw.s)

##############################################
### Get plant data for all timepoints
##############################################
# Modify for stretching factor
stretching = readLines(paste0("/home/henrik/compbio/thesis/data/PNAS/plant", plant.no, "/plant", plant.no, "_tiff_resolutions.txt"))[-c(1:3)]
stretching = sapply(stretching, function(x) strsplit(x, " ")[[1]][5])
stretching = as.numeric(unname(sapply(stretching, function(x) gsub("\\[|\\]|,", "", x))))
all.data$z = sapply(1:ncol(all.data$z), function(time) all.data$z[,time] / stretching[time])

par(mfrow = c(1, 1))
par1 = "z"
par2 = "Mean.cell.intensity"
plot(c(all.data[[par1]]), c(all.data[[par2]]))
nas = which(is.na(c(all.data[[par1]])) | is.na(c(all.data[[par2]]))) 
cor(c(all.data[[par1]])[-nas], c(all.data[[par2]])[-nas], method="pearson")
# heatmap(all.data$z, na.rm = T) # Plot correlation in z-value between timepoints

### Plot z-values
plot(NA, xlim = c(0, 84), ylim = c(0, 12))
times = (0:(ncol(all.data[[1]]) - 1))*4
for(ii in 1:ncol(all.data[[1]])){
  points(rep(times[ii], nrow(all.data$z)), (all.data$z)[,ii], col = 1, pch = 19, type = "p")  
}

