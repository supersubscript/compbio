#!/usr/bin/env Rscript
# Script to plot the total cell volume and the inferred
# cell size, as approximated by 
library(RColorBrewer)
library(scales) #imports alpha
library(stats)
library(plot3D)
palette(brewer.pal(n = 8, name = "Set1"))
color = brewer.pal(n = 8, name = "Set1")
lw.s = 3

# Read in files
files = list.files(path="~/project_out/", pattern="voxels", full.names=F)

# Get the total area
areas = vector("numeric", length=length(files))
count = 1
for(file in files){
  data = read.csv(paste0("~/project_out/", file), sep = ",")
  areas[count] = sum(data$Area)
  count = count + 1
}

# Get the right cell volume
voxel.volume = 4/3*pi*6**3 * 0.5**2 * 1.981818
xy.z.discrepancy = 4

# Plot!
par(mfrow=c(1,1),mar=c(5,5,5,5))
barplot(c(areas[1:8], 0, areas[9:16])/(voxel.volume*xy.z.discrepancy), 
        col = alpha(1:9, .7), names.arg = c(LETTERS[1:8], "", LETTERS[1:8]), yaxt ="n")
mtext("Li",   side=3, line = 1, at = 4.50, cex = 1)
mtext("Otsu", side=3, line = 1, at = 15.5, cex = 1)
axis(2)
axis(4, labels = c("0", "18e6"), at = c(0, 20000))
mtext("Number of cells", side=2, las = 3, line = 2.5)
mtext(paste0("Volume [mikron]"), side=4, las = 3, line = 2.5)
