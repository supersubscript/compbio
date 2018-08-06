#!/usr/bin/env Rscript
setwd("/home/henrik/projects/sum16_solver/MSB_WCK/code/")
# Do the track analysis
library(RColorBrewer)
library(scales) #imports alpha
library(scatterplot3d)
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 3

data = read.csv("../rk.test", skip=0, sep = "\t", header = F)
colnames(data) = c("x", "y", "z", "radius", "WUS", "wus", "CLV3", "clv3", "Y", "KAN1", "kan1", "A", "B", "AnchorG1", "AnchorG2", "AnchorS2a", "AnchorS2b", "nNeighbours")

nTimepoints = data[1,1]
nCells = data[2,1]
cells.trajectories = lapply(1:ncol(data), function(x)
  matrix(ncol = nCells, nrow = nTimepoints))
names(cells.trajectories) = colnames(data)
data = data[which(!is.na(data$nNeighbours)), ]

for(time in 1:nTimepoints) {
  for (species in 1:length(colnames(data))) {
    cells.trajectories[[species]][time, ] = data[((time - 1) * nCells + 1):(time * nCells), species]
  }
}

matplot((cells.trajectories[["CLV3"]][1:100,]), type = "l", col = 1:8)
matplot((cells.trajectories[["clv3"]][1:100,]), type = "l", col = 1:8)
