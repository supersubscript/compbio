#!/usr/bin/env Rscript
# Do the track analysis
setwd("~/compbio/src/assignments/bia1/code/")
library(RColorBrewer)
library(scales) #imports alpha
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 3

# Read in data
# 14 if lower & so
track.file = "../data/All Spots statistics14.csv"
data = read.csv(track.file, sep = ",")
data = data[,c(2,4,5,6,7,9,14,15,16,17,18)] # Take away columns we don't want

############################################
### Plot spots / time
############################################
no.cells.per.time = unname(table(data[,"POSITION_T"]))
plot(1:92, no.cells.per.time, type = "l", yaxt = "l", ylab = "Number of cells", xlab = "Time", lwd = lw.s, col = 1)
axis(2)

############################################
### Find the tracks and filter accordinly
############################################
no.tracks = length(unique(data[,"TRACK_ID"])) # How many tracks do we have? 
cat("no tracks throughout: ", no.tracks, "\n")

# Take away the ones that weren't there to begin with, as well as the ones
# marked as "None" -- these correspond to crappy, short-lived tracks that
# shouldn't be accounted for.
init.filter = which(data$POSITION_T == 0) 
ids = unique(data[init.filter, "TRACK_ID"]) # Retrieve the actual tracks present at t = 0.
ids = ids[which(ids != "None")] # Remove Nones.
no.tracks = length(unique(data[which(data$TRACK_ID %in% ids),"TRACK_ID"])) # How many tracks do we have? 
cat("no tracks in the beginning: ", no.tracks, "\n")

# Find all tracks where cells disappear during the course of the
# time series..
all.splits = c()
for(ii in ids) {
  cells.in.id = which(data$TRACK_ID == ii) 
  # If faulty, remove and go to next track.
  if (length(unname(table(data[cells.in.id, "POSITION_T"]))) != 92 |
      any(diff(table(data[cells.in.id, "POSITION_T"])) < 0))
  {
    data = data[-cells.in.id, ]
    next
  }
  # Collect all the splits we have left
  all.splits = c(all.splits, sum(diff(table(data[cells.in.id, "POSITION_T"]))))
}

output = c(length(all.splits), mean(all.splits), sd(all.splits) / sqrt(length(all.splits)))

#############################################
# Redo analysis above, but allow for disappearances.
#############################################
# Find all tracks where cells disappear throughout
data = read.csv(track.file, sep = ",")
data = data[,c(2,4,5,6,7,9,14,15,16,17,18)]
colnames(data)[c(2, 4:6)] = c("id", "x","y","t")
all.splits = c()

for(ii in ids) {
  cells.in.id = which(data$id == ii)
  # How many splits based on how many are there in the end? 
  l = table(data[cells.in.id, "t"])
  all.splits = c(all.splits, l[length(l)] - 1)
}
output = rbind(output, c(length(all.splits), mean(all.splits), sd(all.splits) / sqrt(length(all.splits))))

#####################################################################
# Redo analysis, but only remove cells which fall outside of the edge
#####################################################################
# Find all tracks where cells disappear throughout
data = read.csv(track.file, sep = ",")
data = data[, c(2, 4, 5, 6, 7, 9, 14, 15, 16, 17, 18)]
colnames(data)[c(2, 4:6)] = c("id", "x","y","t")
all.splits = c()
width = 1100
height = 700
pix.per.um = 1100/709.5 # ca 1.55 pixels / mikron
limit = 30 # pixels

for(ii in ids) {
  cells.in.id = which(data$id == ii)
  if(any(data[cells.in.id, "x"] < 0 + limit | data[cells.in.id, "x"] > width  - limit | 
         data[cells.in.id, "y"] < 0 + limit | data[cells.in.id, "y"] > height - limit)){
    data = data[-cells.in.id, ]
    next
  }
      
  # How many splits based on how many are there in the end? 
  all.splits = c(all.splits, (table(data[cells.in.id, "t"]))[length(table(data[cells.in.id, "t"]))] - 1)
}
output = rbind(output, c(length(all.splits), mean(all.splits), sd(all.splits) / sqrt(length(all.splits))))

output[c(2,1),] = output[c(1,2), ]
output[c(2,3),] = output[c(3,2), ]
print(output)

