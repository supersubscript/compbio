#!/usr/bin/env Rscript
# Plot all the maxima we have inferred from our Find Maxima approach
setwd("~/compbio/src/assignments/bia1/code/")
library(RColorBrewer)
library(scales) #imports alpha
library(stats)
library(plot3D)
palette(brewer.pal(n = 8, name = "Set1"))
color = brewer.pal(n = 8, name = "Set1")
lw.s = 3

# Load files and modify colours
files = list.files(path="~/project_out/", pattern="res", full.names=F)
alphs = c(seq(1,.1, length.out=50), rev(seq(1,.1, length.out=50))) 
colors = alpha(colorRampPalette(color[c(1,2,3,2,1)])(100), alphs)

# Do Li
par(mfrow=c(4,2),mar = c(0,0,0,0), oma = c(0, 0, 0, 0) + 0.0, mai = c(.0, .0, .0, .0))
count = 1
for (file in files[1:8]){
  data = read.csv(paste0("~/project_out/", file), sep = "\t")
  colnames(data) = c("Index", "x", "y", "z")
  cols = cut(data$z, 100) # Colour by z-value
  scatterplot3d(x=data$x*.5, y=data$y*.5, z=data$z*1.981818, color=colors[cols], xlab="", ylab="", zlab="", axis = F, asp = .8, pch = 20, angle=120, main = "", mar = c(0,0,0,0))
  mtext(side=2, LETTERS[count], line = -8, las = 1, at=7, cex = 1.5)  
  count = count + 1
}

# Do Otsu
for (file in files[9:16]){
  data = read.csv(paste0("~/project_out/", file), sep = "\t")
  colnames(data) = c("Index", "x", "y", "z")
  cols = cut(data$z, 100)
  scatterplot3d(x=data$x*.5, y=data$y*.5, z=data$z*1.981818, color=colors[cols], xlab="", ylab="", zlab="", axis = F, asp = .8, pch = 20, angle=120, main = "", mar = c(0,0,0,0))
  mtext(side=2, LETTERS[count], line = -8, las = 1, at=7, cex = 2)  
  count = count + 1
}
