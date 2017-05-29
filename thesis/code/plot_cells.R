#!/usr/bin/env Rscript
setwd("/home/henrik/projects/sum16_solver/Lisa_expPSF_tiffs_SF_Corrected/plant#4/results")
# Do the track analysis
library(RColorBrewer)
library(scales) #imports alpha
library(scatterplot3d)
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 3


data = read.table("t76.txt", skip=1) # Take in whichever file you want

colnames(data) = c("id", "x", "y", "z", "volume", "MeanIntensity")
scatterplot3d(data$x, data$y, data$z)
