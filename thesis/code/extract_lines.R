#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/")
# Do the track analysis
# library(RColorBrewer)
# library(reshape)
# library(scales) #imports alpha
# library(stringr)
# # source("https://bioconductor.org/biocLite.R")
# # biocLite("GenomicRanges")
# palette(brewer.pal(n = 8, name = "Set1"))
# lw.s = 3
# #lapply(plant.files, function(x) read.csv(textConnection(gsub(":", ",", readLines(x))), sep=",", header=FALSE))
# 
# plant.files = list.files("../data/PNAS/plant13/tracking_data/", full.names = TRUE)
# 
# # Go through all cells in output data and find who the mother is
# 
# data = readLines(textConnection(gsub(":", ",", readLines(plant.files[1]))))
# data = lapply(data, function(x) gsub(",", "", x))
# data = lapply(data, function(x) as.integer(str_split(x," ")[[1]]))
# 
# cell.lines = lapply(lapply(data, function(x) x[1]), function(x) list(x)) # These are all the original cell lines
# for(ii in 1:length(data)){
#   for(jj in 1:length(cell.lines)){
#     # If mother is in last time point
#     if(data[[ii]][1] %in% cell.lines[[jj]][[1]]){
#       cell.lines[[jj]] = c(cell.lines[[jj]], list(data[[ii]][-1])) # Add daughter cells
#       break
#     } 
#   }
# }
# cell.lines[sapply(cell.lines, length) < 2] = lapply(cell.lines[sapply(cell.lines, length) < 2], function(x) c(x, list(NA)))
# 
# ### Now do all but the first
# for (time in 2:length(plant.files)) {
#   data = readLines(textConnection(gsub(":", ",", readLines(plant.files[time]))))
#   data = lapply(data, function(x) gsub(",", "", x))
#   data = lapply(data, function(x) as.integer(str_split(x," ")[[1]]))
#   
#   # Now sort out the cell lines
#   for (ii in 1:length(data)) {
#     added = FALSE
#     for (jj in 1:length(cell.lines)) {
#       # If mother is in last time point
#       if (data[[ii]][1] %in% cell.lines[[jj]][[time]]) {
#         # If next timeframe doesn't already exist, add it
#         if (length(cell.lines[[jj]]) != time + 1)
#           cell.lines[[jj]] = c(cell.lines[[jj]], list(data[[ii]][-1])) # Add daughter cells to next time frame
#         else
#           # If it exists, add to it
#           cell.lines[[jj]][[time + 1]] = c(cell.lines[[jj]][[time + 1]], data[[ii]][-1]) # Add daughter cells to next time frame
#         added = TRUE
#         break
#       }
#     }
#     # If we didn't find it
#     if(!added){
#       new.line = list(c(sapply(1:(time-1), function(x) list(NA)), list(data[[ii]][1])))
#       if(length(data[[ii]]) > 1)
#         new.line[[1]] = c(new.line[[1]], list(data[[ii]][-1]))
#       cell.lines = c(cell.lines, new.line)
#     }
#   }
#   
#   # Fill up space for the ones that didn't have daughters
#   cell.lines[sapply(cell.lines, length) == time] = lapply(cell.lines[sapply(cell.lines, length) == time], function(x) c(x, list(NA)))
# }
# 
# # Make it a data frame of lists instead
# cell.lines = t(sapply(cell.lines, function(x) x))
# 
# # Make all integer(0) NAs
# for(ii in 1:ncol(cell.lines)){
#   cell.lines[which(!(sapply(cell.lines[,ii], length))), ii] = list(NA)
# }


# # The tracks that are there for all timepoints
# there.throughout = lapply(1:nrow(cell.lines), function(x) list())
# for(ii in 1:nrow(cell.lines)){
#   for(jj in 1:22){
#     if(!(NA %in% cell.lines[ii, jj][[1]])){
#       there.throughout[[ii]] = c(there.throughout[[ii]], list(cell.lines[ii, jj][[1]]))
#       cat(cell.lines[ii, jj][[1]], " ")
#     }
#   }
#   # cat("\n")
# }
# there.throughout[lapply(there.throughout, length) == 22]

