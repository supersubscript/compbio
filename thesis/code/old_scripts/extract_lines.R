#!/usr/bin/env Rscript
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

args = commandArgs(trailingOnly=TRUE)
plant.no = args[1]
# plant.no =2  # 1, 2, 4, 13, 15, 18
plant.files = list.files(paste0("../data/PNAS/plant", plant.no, "/tracking_data"), full.names = TRUE)

# Put in right order
numbers = sapply(plant.files, function(x) strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])]) 
numbers = unname(sapply(numbers, function(x) as.integer(strsplit(x, "hrs")[[1]][1])))
plant.files = plant.files[order(numbers)]

# Go through all cells in output data and find who the mother is
data = readLines(textConnection(gsub(":", ",", readLines(plant.files[1]))))
data = lapply(data, function(x) gsub(",", "", x))
data = lapply(data, function(x) as.integer(str_split(x," ")[[1]]))

cell.lines = lapply(lapply(data, function(x) x[1]), function(x) list(x)) # These are all the original cell lines

for(ii in 1:length(data)){
  for(jj in 1:length(cell.lines)){
    # If mother is in last time point
    if(data[[ii]][1] %in% cell.lines[[jj]][[1]]){
      cell.lines[[jj]] = c(cell.lines[[jj]], list(data[[ii]][-1])) # Add daughter cells
      break
    }
  }
}
cell.lines[sapply(cell.lines, length) < 2] = lapply(cell.lines[sapply(cell.lines, length) < 2], function(x) c(x, list(NA)))

### Now do all but the first
for (time in 2:length(plant.files)) {
  data = readLines(textConnection(gsub(":", ",", readLines(plant.files[time]))))
  data = lapply(data, function(x) gsub(",", "", x))
  data = lapply(data, function(x) as.integer(str_split(x," ")[[1]]))

  # Now sort out the cell lines
  for (ii in 1:length(data)) {
    added = FALSE
    for (jj in 1:length(cell.lines)) {
      # If mother is in last time point
      if (data[[ii]][1] %in% cell.lines[[jj]][[time]]) {
        # If next timeframe doesn't already exist, add it
        if (length(cell.lines[[jj]]) != time + 1)
          cell.lines[[jj]] = c(cell.lines[[jj]], list(data[[ii]][-1])) # Add daughter cells to next time frame
        else
          # If it exists, add to it
          cell.lines[[jj]][[time + 1]] = c(cell.lines[[jj]][[time + 1]], data[[ii]][-1]) # Add daughter cells to next time frame
        added = TRUE
        break
      }
    }
    # If we didn't find it
    if(!added){
      new.line = list(c(sapply(1:(time-1), function(x) list(NA)), list(data[[ii]][1])))
      if(length(data[[ii]]) > 1)
        new.line[[1]] = c(new.line[[1]], list(data[[ii]][-1]))
      cell.lines = c(cell.lines, new.line)
    }
  }

  # Fill up space for the ones that didn't have daughters
  cell.lines[sapply(cell.lines, length) == time] = lapply(cell.lines[sapply(cell.lines, length) == time], function(x) c(x, list(NA)))
}

# Make it a data frame of lists instead
cell.lines = t(sapply(cell.lines, function(x) x))

# Make all integer(0) NAs
for(ii in 1:ncol(cell.lines)){
  cell.lines[which(!(sapply(cell.lines[,ii], length))), ii] = list(NA)
}


# Extract the tracks that are there for all timepoints
# there.throughout = lapply(1:nrow(cell.lines), function(x) list())
# for(ii in 1:nrow(cell.lines)){
#   for(jj in 1:length(cell.lines[ii, ])){
#     if(!(NA %in% cell.lines[ii, jj][[1]])){
#       there.throughout[[ii]] = c(there.throughout[[ii]], list(cell.lines[ii, jj][[1]]))
#       # cat(cell.lines[ii, jj][[1]], " ")
#     }
#   }
#   # cat("\n")
# }
# there.throughout[lapply(there.throughout, length) == 22]


# Put the quantified files in the right order
quant.files = list.files(paste0("/home/henrik/compbio/thesis/data/msb_plants/plant", plant.no, "/results"), full.names = TRUE)
numbers = sapply(quant.files, function(x) strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])]) 
numbers = unname(sapply(numbers, function(x) as.integer(gsub("t|\\.txt", "", x))))
quant.files = quant.files[order(numbers)]

# Read in data
data = lapply(quant.files, function(x) read.table(x, header=T, sep = "\t"))

# Find all the expression data, volume, etc. for the lineages
cell.lines.data =  vector("list", nrow(cell.lines))  # prealloc
for(ii in 1:nrow(cell.lines)){
  lineage.data = matrix(NA, nrow = 0, ncol = 1 + ncol(data[[1]]))
  for(jj in 1:length(quant.files)){
    cells.in.lineage = which(data[[jj]]$Cell.id %in% cell.lines[[ii, jj]])
    if(length(cells.in.lineage) > 0)
      for(kk in 1:length(cells.in.lineage))
        lineage.data = rbind(lineage.data, c(Time=(jj-1)*4, as.vector(data[[jj]][cells.in.lineage[kk], ])))
    else
      lineage.data = rbind(lineage.data, rep(NA, ncol(lineage.data)))
  }
  cell.lines.data[[ii]] = lineage.data
}

len = sapply(cell.lines.data, function(x) length(x[!is.na(x)]))
cell.lines.data = cell.lines.data[order(-len)]
cell.lines.data = cell.lines.data[sort(len, decreasing = T) > 0]

# plot(NA, xlim = c(0,84), ylim = c(0,100))
# parameter = "x"
# parameter = "y"
parameter = "z"
# parameter = "Boa.volume"
# parameter = "Mean.cell.intensity"
plot(NA, xlim = c(0, 84), ylim = c(0, 10))
for(ii in 1:length(cell.lines.data)){
  points(cell.lines.data[[ii]][, "Time"], unlist(cell.lines.data[[ii]][, parameter]), col = 1, pch = 19, type = "p")  
}

