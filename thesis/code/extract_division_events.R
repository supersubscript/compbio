#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/")
# Do the track analysis
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 3

get.plant.data = function(plant.no){
  library(reshape)
  library(scales) #imports alpha
  library(stringr)
  library(plyr)
  library(parallel)

# Read in the right files
plant.files = list.files(paste0("../data/PNAS/plant", plant.no, "/tracking_data"), full.names = TRUE)
quant.files = list.files(paste0("/home/henrik/compbio/thesis/data/msb_plants/plant", plant.no, "/results"), full.names = TRUE)

# Put in right order
order.plant = sapply(plant.files, function(x) strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])]) 
order.plant = unname(sapply(order.plant, function(x) as.integer(strsplit(x, "hrs")[[1]][1])))
plant.files = plant.files[order(order.plant)]
order.quant = sapply(quant.files, function(x) strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])]) 
order.quant = unname(sapply(order.quant, function(x) as.integer(gsub("t|\\.txt", "", x))))
quant.files = quant.files[order(order.quant)]

# Read in quantified data
quant.data = lapply(quant.files, function(x) read.table(x, header=T, sep = "\t"))
params = colnames(quant.data[[1]])
all.quant.data = lapply(1:ncol(quant.data[[1]]), function(x) t(plyr::ldply(unname(sapply(quant.data, "[", x)), rbind)))
names(all.quant.data) = params

output.data = matrix(NA, ncol=15, nrow=0)
colnames(output.data) = c("cell.id", "lineage.id", "mother.id", "daughter1.id", 
                          "daughter2.id", "age", "t","x","y","z", "vol", "expr", 
                          "mother.expr", "daughter1.expr", "daughter2.expr")

### Read in first generation
data = readLines(textConnection(gsub(":", ",", readLines(plant.files[1]))))
data = lapply(data, function(x) gsub(",", "", x))
data = lapply(data, function(x) as.integer(str_split(x," ")[[1]]))
cell.lines = lapply(lapply(data, function(x) x[1]), function(x) list(x)) # These are all the original cell lines

# Add daughter cells to time 2
for(ii in 1:length(data)){
  for(jj in 1:length(cell.lines)){
    if(data[[ii]][1] %in% cell.lines[[jj]][[1]]){
      cell.lines[[jj]] = c(cell.lines[[jj]], list(data[[ii]][-1])) # Add daughter cells
      break
    }
  }
}

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

# Read in data
data = lapply(quant.files, function(x) read.table(x, header=T, sep = "\t"))

# Find all the expression data, volume, etc. for the lineages
cell.lines.data =  vector("list", nrow(cell.lines))  # prealloc
for (ii in 1:nrow(cell.lines)) {
  lineage.data = matrix(NA, nrow = 0, ncol = 1 + ncol(data[[1]]))
  for (jj in 1:length(plant.files)) {
    cells.in.lineage = which(data[[jj]]$Cell.id %in% cell.lines[[ii, jj]])
    if (length(cells.in.lineage) > 0)
      for (kk in 1:length(cells.in.lineage))
        lineage.data = rbind(lineage.data, c(Time = (jj - 1) * 4, as.vector(data[[jj]][cells.in.lineage[kk],])))
      else
        lineage.data = rbind(lineage.data, rep(NA, ncol(lineage.data)))
  }
  cell.lines.data[[ii]] = lineage.data
}

len = sapply(cell.lines.data, function(x) length(x[!is.na(x)]))
cell.lines.data = cell.lines.data[order(-len)]
cell.lines.data = cell.lines.data[sort(len, decreasing = T) > 0]

############################
### Back to normal
############################
all.division.data = lapply(1:length(plant.files), function(x) {
  data = readLines(textConnection(gsub(":", ",", readLines(plant.files[x]))))
  data = lapply(data, function(y) gsub(",", "", y))
  data = lapply(data, function(y) as.integer(str_split(y," ")[[1]]))
  data
})

# Get the mother cell of all cells
mothers = lapply(1:length(all.division.data), function(x)
  lapply(all.division.data[[x]], function(y) {
    mother.id = y[1] #rep(y[1], length(cell.id))
    cell.id = y[-1]
    if (length(cell.id) == 1)
      return(matrix(c(t=x, mother.id, cell.id), ncol = 3, nrow = 1))
    
    # print(unname(cbind(t = x, mother.id, cell.id)))
    mums = rep(mother.id, length(cell.id))
    ts = rep(x, length(cell.id))
    return(cbind(ts, mums, cell.id))
  }))
mothers = do.call(rbind, unlist(mothers, recursive = F))
first.gen = cbind(rep(0, length(all.division.data[[1]])),
                  rep(NA, length(all.division.data[[1]])),
                  sapply(all.division.data[[1]], "[[", 1))
mothers = rbind(first.gen, mothers)
colnames(mothers) = c("t", "mother.id", "cell.id")

for(time in 1:length(plant.files)){
  division.data = all.division.data[[time]]
  # division.data = readLines(textConnection(gsub(":", ",", readLines(plant.files[time]))))
  # division.data = lapply(division.data, function(x) gsub(",", "", x))
  # division.data = lapply(division.data, function(x) as.integer(str_split(x," ")[[1]]))
  division.data = division.data[sapply(division.data, length) > 2] # Only keep division events

  # for every cell, take out the data FROM all.quant.data that is interesting for that and add it to a row. Append this to the output matrix.
  for(ii in 1:length(division.data)){
    this.cell.index = all.quant.data$Cell.id[, time] %in% division.data[[ii]][1] #
    daughter1.in = all.quant.data$Cell.id[, time] %in% division.data[[ii]][2]
    daughter2.in = all.quant.data$Cell.id[, time] %in% division.data[[ii]][3]
    xyzvi = rep(NA, 5)
    
    # Does this cell exist in the quantified data?
    if (any(this.cell.index)) {
      xyzvi = c(
        all.quant.data$x[, time][which(this.cell.index)],
        all.quant.data$y[, time][which(this.cell.index)],
        all.quant.data$z[, time][which(this.cell.index)],
        all.quant.data$Boa.volume[, time][which(this.cell.index)],
        all.quant.data$Mean.cell.intensity[, time][which(this.cell.index)])}

    # Do the daughter cells exist in the quantified data?
    d1 = ifelse(any(daughter1.in), all.quant.data$Mean.cell.intensity[, time][which(daughter1.in)], NA)
    d2 = ifelse(any(daughter2.in), all.quant.data$Mean.cell.intensity[, time][which(daughter2.in)], NA)

    # Get the identifiers
    lineage.id = which(sapply(1:nrow(cell.lines),  function(x) division.data[[ii]][1] %in% cell.lines[[x, time]]))
    this.cell.id = division.data[[ii]][1]
    mother.id = (mothers[mothers[, "cell.id"] == this.cell.id & mothers[,"t"] == time - 1, "mother.id"]) # Which is the mother id for this cell?
    mother.id = ifelse(length(mother.id) == 0, NA, mother.id)
    mother.expr = output.data[which(output.data[,"t"] == time - 4 & 
                                    output.data[,"cell.id"] == mother.id), "expr"]
    mother.expr = ifelse(length(mother.expr) == 0, NA, mother.expr)
    
    # Find how long the mother has been alive for
    age = 0
    latest = mother.id
    # reference.id = unname(mothers[mothers[, "cell.id"] == mother.id & mothers[,"t"] == time - 1, "mother.id"]) # grandmother
    for (kk in time:2) {
      # Check the mother. Does she have several offspring?
      reference.id = unname(mothers[mothers[, "cell.id"] == latest & mothers[,"t"] == kk - 2, "mother.id"])
      if(length(reference.id[!is.na(reference.id)]) == 0)
        break
      present.previous = mothers[,"t"] == kk - 1
      # print(reference.id)
      ancestor = mothers[, "mother.id"] == reference.id
      no.offspring = length(mothers[ancestor & present.previous, ])
      # print(no.offspring)
      latest = reference.id
      
      # If the mother to the reference has more than one offspring, then a division event has taken place
      if (no.offspring > 1) {
        break
      }
      age = age + 1
    }

    # Print output data
    output.data = rbind(output.data, c(this.cell.id, # the id of the cell in this timepoint
                                       lineage.id, # the cell lineage id
                                       mother.id, # mother id in previous timepoint
                                       division.data[[ii]][2], #daughter1 id
                                       division.data[[ii]][3], # daughter2 id
                                       age,
                                       (time - 1) * 4,
                                       xyzvi, mother.expr, d1, d2))
  }
}
output.data
}
plants = c(2, 4, 13, 15, 18)

library(parallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
results = parLapply(cl, plants, get.plant.data)
stopCluster(cl)


# plot(output.data[, "t"], output.data[, "expr"])
# 
# hist(output.data[,"expr"], breaks=30, col = 1)
# plot(density(output.data[,"expr"][!is.na(output.data[,"expr"])], kernel = "gaussian"), col = 1, lwd = lw.s)
# hist(output.data[,"t"], breaks=30, col = 1)
# hist(output.data[, "vol"], breaks=30)
# 
# cor(output.data[, "t"][!is.na(output.data[, "x"])], output.data[, "x"][!is.na(output.data[, "x"])])
