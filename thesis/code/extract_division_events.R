#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/")
library(parallel)
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 2
par(lwd=lw.s, cex = 1.5, ps = 15)
interval = 4
maxtime = 84

# TODO: Plant 18 has a gap in one of the time series 
# source("normalise_coordinates.R")

# get.plant.data = function(plant.no){
  # Load the necessary packages
  packages = c("reshape", "scales", "stringr", "plyr", "parallel")
  lapply(packages, library, character.only = TRUE)
  
  # Read in the right files
  mapping.files = list.files(paste0("../data/PNAS/plant", plant.no, "/tracking_data"), full.names = TRUE)
  quant.files = list.files(paste0("../data/clv3_complete/plant", plant.no, "/Results"), full.names = TRUE, pattern=".txt")
  timepoints = unname(sapply(sapply(sapply(mapping.files, function(x) tail(strsplit(x, "/")[[1]], 1)), function(x) gsub("hrs|\\.txt", "", x)), function(x) as.integer(tail(strsplit(x, "_")[[1]], 1))))
  
  no.track.files = length(mapping.files)
  no.quant.files = length(quant.files)
  
  # Put in right order
  order.plant = sapply(mapping.files, function(x) strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])]) 
  order.plant = unname(sapply(order.plant, function(x) as.integer(strsplit(x, "hrs")[[1]][1])))
  order.quant = sapply(quant.files, function(x) strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])]) 
  order.quant = unname(sapply(order.quant, function(x) as.integer(gsub("t|\\.txt", "", x))))
  
  mapping.files = mapping.files[order(order.plant)]
  quant.files = quant.files[order(order.quant)]
  timepoints = c(0, timepoints[order(order.plant)])
  
  # Read in the mapping data between timepoints
  all.mapping.data = lapply(1:length(mapping.files), function(x) {
    data = readLines(textConnection(gsub(":", ",", readLines(mapping.files[x]))))
    data = lapply(data, function(y) gsub(",", "", y))
    data = lapply(data, function(y) as.integer(str_split(y," ")[[1]]))
    data
  })
  # Read in quantified data
  quant.data = lapply(quant.files, function(x) read.table(x, header=T, sep = "\t"))
  if(length(quant.data) != 0){
    all.quant.data = lapply(1:ncol(quant.data[[1]]), function(x) t(ldply(unname(sapply(quant.data, "[", x)), rbind)))
    names(all.quant.data) = colnames(quant.data[[1]])
  }
  
  id.files = list.files(paste0("../data/plant", plant.no, "_correspondence"), full.names = T, pattern="correspondence")
  id.files.order = sapply(sapply(id.files, function(x) tail(strsplit(x, "_")[[1]], 1)), function(x) as.integer(gsub("\\.txt", "", x)))
  id.files = id.files[order(id.files.order)]
  id.maps = lapply(id.files, function(x) read.table(x, header = F, sep = ","))
  
  # TODO: Some cells have multiple nuclei
  for (ii in 1:length(all.mapping.data)) {
    for (jj in 1:length(all.mapping.data[[ii]])) {
      id = all.mapping.data[[ii]][[jj]][1]
      row = which(id.maps[[ii]][, 1] == id)[1]
      if (is.na(row)) {
        all.mapping.data[[ii]][[jj]][1] = id + 10000
      } else {
        all.mapping.data[[ii]][[jj]][1] = id.maps[[ii]][row, 2]
      }
      all.mapping.data[[ii]][[jj]][-1] = 
        sapply(all.mapping.data[[ii]][[jj]][-1], function(id) {
        row = which(id.maps[[ii + 1]][, 1] == id)[1]
        if (is.na(row))
          return(id + 10000)
        else
          return(id.maps[[ii + 1]][row, 2])
      })
    }
  }
  
  
  #############################################################################
  # Extract cell lines
  #############################################################################
  cell.lines = cbind(cbind(t(t(cbind(lapply(all.mapping.data[[1]], "[[", 1), 
                                     lapply(all.mapping.data[[1]], "[", -1))))), 
                     matrix(list(NA), ncol = no.track.files - 2 + 1, 
                            nrow = length(all.mapping.data[[1]])))
  # colnames(cell.lines) = timepoints
  for(time in 2:no.track.files){
    this.mapping.data = all.mapping.data[[time]]
    parents = lapply(this.mapping.data, "[[", 1)
    children = lapply(this.mapping.data, "[", -1)
    parents.found = which(sapply(parents, function(x) x %in% unlist(cell.lines[, time])))
    parents.not.found = !1:length(parents) %in% parents.found
    
    # Add the children of the parents' which are found
    for(ii in parents.found){
      for(line in 1:nrow(cell.lines)){
        if(parents[[ii]] %in% cell.lines[[line, time]]){
          already.there = cell.lines[[line, time + 1]]
          cell.lines[[line, time + 1]] = c(already.there[!is.na(already.there)], as.integer(children[[ii]]))
          break
        }
      }
    }
    
    # Create new cell lines for the ones which aren't found
    new.lines = cbind(cbind(
      matrix(list(NA), ncol = time - 1, nrow = length(all.mapping.data[[time]][-parents.found])), 
      t(t(cbind(lapply(all.mapping.data[[time]][-parents.found], "[[", 1), 
                lapply(all.mapping.data[[time]][-parents.found], "[", -1))))), 
      matrix(list(NA), ncol = no.track.files - time, nrow = length(all.mapping.data[[time]][-parents.found])))
    cell.lines = rbind(cell.lines, new.lines)
  }
  cell.lines = cell.lines[order(-apply(cell.lines, 1, function(x) sum(!is.na(unlist(x))))), ]
  
  #############################################################################
  # Extract cell line data
  #############################################################################
  # Find all the expression data, volume, etc. for the lineages
  if(plant.no != 1){
    cell.lines.data =  vector("list", nrow(cell.lines))  # prealloc
    for (ii in 1:nrow(cell.lines)) {
      lineage.data = matrix(NA, nrow = 0, ncol = 1 + ncol(quant.data[[1]]) + 1) # add time and dist2top
      for (time in 1:no.track.files) {
        cells.in.lineage = which(quant.data[[time]]$Cell.id %in% cell.lines[[ii, time]])
        topcell = quant.data[[time]][order(-quant.data[[time]]$z), 2:4][1, ] # Take top 3 cells
        top.coords = apply(topcell, 2, mean)
        top.coords[1] = mean(quant.data[[time]]$x)
        top.coords[2] = mean(quant.data[[time]]$y)
        
        if (length(cells.in.lineage) > 0){
          for (kk in 1:length(cells.in.lineage)){
            dist2tops = sqrt(as.vector(quant.data[[time]][cells.in.lineage[kk],]$x - top.coords[1]) ** 2 + 
                             as.vector(quant.data[[time]][cells.in.lineage[kk],]$y - top.coords[2]) ** 2 +
                             as.vector(quant.data[[time]][cells.in.lineage[kk],]$z - top.coords[3]) ** 2)
            lineage.data = rbind(lineage.data, c(Time = timepoints[time], as.vector(quant.data[[time]][cells.in.lineage[kk],]), dist2top=dist2tops))
          }
        } else
          lineage.data = rbind(lineage.data, rep(NA, ncol(lineage.data)))
      }
      cell.lines.data[[ii]] = lineage.data
    }
    # Sort cell lines data by number fo members
    cell.lines.members = sapply(cell.lines.data, function(x) length(x[!is.na(x)]))
    # cell.lines.data = cell.lines.data[order(-cell.lines.members)]
    # cell.lines.data = cell.lines.data[!sapply(cell.lines.data, function(x) sum(!is.na(unlist(x))) < 50)]
    cell.lines.data = lapply(cell.lines.data, function(x) apply(x, 1:2, function(y) unlist(y)))
    # cell.lines.data = cell.lines.data[order(sapply(cell.lines.data, function(x) median(x[,"dist2top"], na.rm = TRUE)))]
  }
  
  #############################################################################
  # Extract mother cells to all cells
  #############################################################################
  # Get the mother cell of all cells
  mothers = lapply(1:length(all.mapping.data), function(x)
    lapply(all.mapping.data[[x]], function(y) {
      mother.id = y[1]
      cell.id = y[-1]
      mums = rep(mother.id, length(cell.id))
      ts = rep(timepoints[x + 1], length(cell.id))
      return(cbind(ts, mums, cell.id))
    }))
  mothers = do.call(rbind, unlist(mothers, recursive = F))
  first.gen = cbind(rep(0, length(all.mapping.data[[1]])),
                    rep(NA, length(all.mapping.data[[1]])),
                    sapply(all.mapping.data[[1]], "[[", 1))
  mothers = rbind(first.gen, mothers)
  colnames(mothers) = c("t", "mother.id", "cell.id")
  output.data = matrix(NA, ncol = 17, nrow = 0)
  colnames(output.data) = c("cell.id", "lineage.id", "mother.id", "daughter1.id", 
                            "daughter2.id", "age", "t","x","y","z", "vol", "expr", 
                            "dist2top", "daughter1.vol", "daughter2.vol", "daughter1.expr", "daughter2.expr")
  
  #############################################################################
  # Get quantified data for all division events
  #############################################################################
  ### Loop through all the time points
  for(time in 1:no.track.files){
    division.data = all.mapping.data[[time]]
    division.data = division.data[sapply(division.data, length) > 2] # Only keep division events
    topcell       = quant.data[[time]][order(-quant.data[[time]]$z), 2:4][1:3, ] # Take top 3 cells
    top.coords    = apply(topcell, 2, mean)
    top.coords[1] = mean(quant.data[[time]]$x)
    top.coords[2] = mean(quant.data[[time]]$y)
       
    # for every cell, take out the data FROM all.quant.data that is interesting for that and add it to a row. Append this to the output matrix.
    for(ii in 1:length(division.data)){
      xyzvid          = rep(NA, 6)
      this.cell.index = NA
      daughter1.index = NA
      daughter2.index = NA
      mother.expr     = NA
      daughter1.vol   = NA
      daughter2.vol   = NA
      daughter1.expr  = NA
      daughter2.expr  = NA
      dist2top        = NA
      
      if(plant.no != 1){
        this.cell.index = all.quant.data$Cell.id[, time] %in% division.data[[ii]][1] 
        daughter1.index = all.quant.data$Cell.id[, time] %in% division.data[[ii]][2]
        daughter2.index = all.quant.data$Cell.id[, time] %in% division.data[[ii]][3]
        
        # Does this cell exist in the quantified data?
        if (any(this.cell.index)) {
          xyzvid = c(
            all.quant.data$x[, time][which(this.cell.index)],
            all.quant.data$y[, time][which(this.cell.index)],
            all.quant.data$z[, time][which(this.cell.index)],
            all.quant.data$Boa.volume[, time][which(this.cell.index)],
            all.quant.data$Mean.cell.intensity[, time][which(this.cell.index)],
            sqrt((all.quant.data$x[, time][which(this.cell.index)] - top.coords[1]) ** 2 +
                 (all.quant.data$y[, time][which(this.cell.index)] - top.coords[2]) ** 2 +
                 (all.quant.data$z[, time][which(this.cell.index)] - top.coords[3]) ** 2)
          )
        }
        
        # Do the daughter cells exist in the quantified data?
        # TODO: fix this again
        if (plant.no == 18 & time > 11){
          daughter1.expr = ifelse(any(daughter1.index), all.quant.data$Mean.cell.intensity[, time][which(daughter1.index)], NA)
          daughter2.expr = ifelse(any(daughter2.index), all.quant.data$Mean.cell.intensity[, time][which(daughter2.index)], NA)
          daughter1.vol  = ifelse(any(daughter1.index), all.quant.data$Boa.volume[, time][which(daughter1.index)], NA)
          daughter2.vol  = ifelse(any(daughter2.index), all.quant.data$Boa.volume[, time][which(daughter2.index)], NA)
          
        } else {
          daughter1.expr = ifelse(any(daughter1.index), all.quant.data$Mean.cell.intensity[, time + 1][which(daughter1.index)], NA)
          daughter2.expr = ifelse(any(daughter2.index), all.quant.data$Mean.cell.intensity[, time + 1][which(daughter2.index)], NA)
          daughter1.vol  = ifelse(any(daughter1.index), all.quant.data$Boa.volume[, time + 1][which(daughter1.index)], NA)
          daughter2.vol  = ifelse(any(daughter2.index), all.quant.data$Boa.volume[, time + 1][which(daughter2.index)], NA)
          
        }
      }
      
      # Get the identifiers
      lineage.id = which(sapply(1:nrow(cell.lines), function(x) division.data[[ii]][1] %in% cell.lines[[x, time]]))
      this.cell.id = division.data[[ii]][1] # remember, id and index are different
      mother.id = mothers[mothers[, "cell.id"] == this.cell.id & mothers[, "t"] == timepoints[time], "mother.id"] # Which is the mother id for this cell?
      mother.id = ifelse(length(mother.id) == 0, NA, mother.id)
 
      ### Find how long the mother has been alive for
      age = 0
      last.mother = mother.id
      for(kk in (time-1):1){
        mother.to.mother = mothers[mothers[, "cell.id"] == last.mother & mothers[,"t"] == timepoints[kk], "mother.id"][1]
        no.offspring = length(which(mothers[, "mother.id"] == mother.to.mother & mothers[,"t"] == timepoints[kk]))
        if(no.offspring > 1 || no.offspring == 0)
          break
        last.mother = mother.to.mother
        age = age + timepoints[kk + 1] - timepoints[kk]
      }
 
      # Print output data
      output.data = rbind(output.data, 
                          c(this.cell.id, 
                            lineage.id, 
                            mother.id, 
                            division.data[[ii]][2], 
                            division.data[[ii]][3], 
                            age,
                            timepoints[time],
                            xyzvid, 
                            daughter1.vol, 
                            daughter2.vol,
                            daughter1.expr, 
                            daughter2.expr))
    }
  }
  # return(output.data)
# }


# plants = c(1, 2, 4, 13, 15, 18)
# no.plants = length(plants)
# no_cores = detectCores() - 1
# cl = makeCluster(no_cores)
# results = parLapply(cl, plants, get.plant.data)
# stopCluster(cl)
# 
# all.data = do.call(rbind, results)
all.data = output.data
#################
# ### Cell age is exponentially distributed, i.e. cells seem to live for an amount of time that depends on some underlying average rate.
# ### Actual division event due to stochasticity? Exponential distribution memoryless -- doesn't make sense for us in the short time scale.
# ### 4 hours enough to reveal that? 
# ################
# hist.fit.exponential = function(data, title) {
#   histogram = hist(data[, "age"], breaks = 100, main = title)
#   is.0 = histogram$counts == 0
#   temp = data.frame(y = histogram$counts[!is.0], x = histogram$mids[!is.0])
#   mod = nls(y ~ exp(a + b * x), data = temp, start = list(a = .5, b = .5))
#   lines(c(0, temp$x), predict(mod, list(x = c(0, temp$x))), col = 1)
# }
# par(mfrow = c(4, 2))
# hist.fit.exponential(all.data, "all")
# for(ii in 1:no.plants) hist.fit.exponential(results[[ii]], plants[ii])
# 
# ################
# ### In big enough intervals, cellular divisions follow a normal distribution with the logged radius. Does the probability to divide 
# ### decay exponentially with the radius? Also radius squared? 
# ################
par(mfrow = c(1, 3))
no.breaks = 30
hist(log(sqrt(all.data[, "x"] ** 2 + all.data[, "y"] ** 2)), breaks = no.breaks, main = "ln(rad)")
hist(((all.data[, "x"] ** 2 + all.data[, "y"] ** 2)), breaks = no.breaks, main = "rad**2")
hist((sqrt(all.data[, "x"] ** 2 + all.data[, "y"] ** 2)), breaks = no.breaks, main = "rad")

# Plot all the quantification data
par(mfrow = c(4, 3))
no.breaks = 30
for(ii in 6:ncol(all.data)) if(!all(is.na(all.data[,ii]))) hist(all.data[, ii], breaks = no.breaks, main = colnames(all.data)[ii], freq = F)
for(ii in 6:ncol(all.data)) if(!all(is.na(all.data[,ii]))) plotDensity(data.frame(x=all.data[, ii]), main = colnames(all.data)[ii])

graphics.off(); plot.new()
par(mfrow = c(3, 3), mar = c(2,2,2,2))
plots = vector(mode="list", length=(ncol(all.data)-5));
interesting = (1:ncol(all.data))[-c(1:5, 8:10)] # Remove xyz and indices
counter = 1
for (ii in interesting) {
  not.ii = interesting[interesting != ii]
  for (jj in not.ii) {
    if (length(all.data[!is.na(all.data[, ii]), ii]) != 0 &
        length(all.data[!is.na(all.data[, jj]), jj]) != 0) {
      x.lab = colnames(all.data)[jj]
      y.lab = colnames(all.data)[ii]
      is.not.na = !is.na(all.data[, y.lab]) & !is.na(all.data[, x.lab])
      correlation = cor(all.data[is.not.na, x.lab] , all.data[is.not.na, y.lab])
      title = paste0(y.lab, " vs ", x.lab, "\n cor: ", round(correlation, 2))
      plot(all.data[, jj], all.data[, ii], main = title, xlab = x.lab, ylab = y.lab)
    }
  }
  plots[[counter]] = recordPlot()
  plot.new()
  counter = counter+1
}
plots[[1]]
plots[[2]]
plots[[3]]
plots[[4]]
plots[[5]]
plots[[6]]
plots[[7]]
plots[[8]]
plots[[9]]

scatter.hist(all.data[,"expr"], all.data[,"dist2top"], x.breaks=50, y.breaks=50)
scatter.hist(all.data[,"expr"], all.data[,"dist2top"], x.breaks=50, y.breaks=50)


# 
# par(mfrow = c(3, 1))
# hist(all.data[,"daughter1.expr"] / all.data[,"expr"], breaks = 50)
# hist(all.data[,"daughter2.expr"] / all.data[,"expr"], breaks = 50)
# hist(all.data[,"daughter1.expr"] / all.data[,"daughter2.expr"], breaks = 50)

###################################
###################################
plot(cell.lines.data[[1]][,"Time"])
par(mfrow = c(4, 3), mar = c(2, 2, 2, 2))
for(ii in 1:12){
  plot(cell.lines.data[[ii]][,"Time"], cell.lines.data[[ii]][, "Mean.cell.intensity"], ylim = c(0,160), xlim = c(0,84))
  par(new = TRUE)
  plot(cell.lines.data[[ii]][,"Time"], cell.lines.data[[ii]][, "dist2top"], col = 2, ylim=c(0,10), xaxt="n", yaxt="n", xlim = c(0, 84))
  axis(4)
}

# plot(cell.lines.data[[1]][,"x"], cell.lines.data[[1]][, "y"], col = 1:22)
# text(cell.lines.data[[1]][,"x"], cell.lines.data[[1]][, "y"], labels = 1:length(cell.lines.data[[1]][, "y"]))

# par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))
# divisions.per.time = (as.integer(table(all.data[,"t"])))
# expression.lvls = sapply(unique(all.data[,"t"]), function(x) mean(all.data[which(all.data[,"t"] == x), "vol"], na.rm=T))
# plot(expression.lvls)
# plot(divisions.per.time, expression.lvls, main=round(cor(divisions.per.time, expression.lvls, use="complete"), 2))
# cor.test(divisions.per.time, expression.lvls, use="complete")
