#!/usr/bin/env Rscript
# Script to extract and plot the images resulting from
# multiplying the processed and the thresholded images.
setwd("~/compbio/src/assignments/bia1/code/")
library(RColorBrewer)
library(scales) #imports alpha
library(stats)
library(tiff)
library(raster)
palette(brewer.pal(n = 8, name = "Set1"))

# Adapt this accordingly
files = list.files(path="~/project_out", pattern="multiplied", full.names=T)

# Li
count = 1
par(mfrow=c(4,2),mar = c(0,0,0,0), oma = c(1, 1, 1, 1) + 0.0, mai = c(.0, .0, .0, .0))
for(file in files[1:8]){
  f = readTIFF(file, all=T, native = T)[[50]] # Take slice 50
  image(f, axes=F, col = c("black", "white"))  
  mtext(side=2, LETTERS[count], line = -3, las = 1, cex = 2, at = .9) 
  count = count + 1
}

# Otsu
count = 1
par(mfrow=c(4,2),mar = c(0,0,0,0), oma = c(1, 1, 1, 1) + 0.0, mai = c(.0, .0, .0, .0))
for(file in files[9:16]){
  f = readTIFF(file, all=T, native = T)[[50]] 
  image(f, axes=F, col = c("black", "white"))  
  mtext(side=2, LETTERS[count], line = -3, las = 1, cex = 2, at = .9)  
  count = count + 1
}

