#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/")
# Do the track analysis
library(RColorBrewer)
library(reshape)
library(scales) #imports alpha
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicRanges")
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 3
#lapply(plant.files, function(x) read.csv(textConnection(gsub(":", ",", readLines(x))), sep=",", header=FALSE))

plant.files = list.files("../data/PNAS/plant13/tracking_data/", full.names = TRUE)




file = plant.files[4]
data = read.csv(textConnection(gsub(":", ",", readLines(file))), sep = ",", header = FALSE)
get.tracks = function(data) {
  # Load data
  # data = read.csv(textConnection(gsub(":", ",", readLines(file))), sep = ",", header =
  # FALSE)
  # 
  # # Take out duplicates, remove NAs
  # data = unname(split(data, seq(nrow(data))))
  # data = lapply(data, function(x)
  #   unique(unlist(unname(x[!is.na(x)]))))
  
  removal = c()
  for (ii in 1:length(data)) {
    for (jj in ii:length(data)) {
      if (ii == jj)
        next
      if (any(data[[ii]] %in% data[[jj]])) {
        data[[jj]] = union(data[[ii]], data[[jj]]) # propagate down
        removal = append(removal, ii)
        break
      }
    }
  }
  data = data[-removal]
  return(data)
}

# get.tracks(data)




# cell.lines = list()
# for(file in 1:length(plant.files)){
#   data = read.csv(textConnection(gsub(":", ",", readLines(plant.files[file]))), sep = ",", header = FALSE)
#   data = unname(split(data, seq(nrow(data))))
#   data = lapply(data, function(x)
#     unique(unlist(unname(x[!is.na(x)]))))
#   
#     # for(ii in data){
#   #   for(jj in cell.lines){
#   #     if(ii[length(ii)] %in% jj){
#   #       jj = append(jj, ii[-1])
#   #       break
#   #     }
#   #   }
#   #   cell.lines = c(cell.lines, list(jj))
#   # }
# }








# 
# file.names = sapply(list.files("../data/PNAS/plant13/tracking_data/", full.names = FALSE), function(x)
#   gsub(".txt", "", x, fixed = TRUE))
# track.data = get.tracks(plant.files)#lapply(plant.files, get.tracks)
# names(track.data) = file.names
# 
# ########################################
# results.files = list.files("/home/henrik/projects/sum16_solver/data/plant#13/results", full.names = TRUE)
# segm.data  = lapply(results.files, function(x) read.csv(x, sep="\t", header=TRUE))
# 
# file = 21
# cell.lineage = 1
# segm.data[[file]][which(segm.data[[file]][,1] %in% track.data[[file]][[cell.lineage]]), "Mean.cell.intensity"]
# 




