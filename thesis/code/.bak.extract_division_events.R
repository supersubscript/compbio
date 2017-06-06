#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/")
library(parallel)
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 2
par(lwd = lw.s, cex = 1.5, ps = 15)

# TODO: Plant 18 has a gap in one of the time series 
# TODO: We don't really need mother expression in output table

get.plant.data = function(plant.no){
  # Load the necessary packagesa
  packages = c("reshape", "scales", "stringr", "plyr", "parallel")
  lapply(packages, library, character.only = TRUE)
  
  # Read in the right files
  mapping.files = list.files(paste0("../data/PNAS/plant", plant.no, "/tracking_data"), full.names = TRUE)
  quant.files = list.files(paste0("../data/msb_plants/plant", plant.no, "/results"), full.names = TRUE)
  
  no.track.files = length(mapping.files)
  no.quant.files = length(quant.files)
  
  # Put in right order
  order.plant = sapply(mapping.files, function(x) strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])]) 
  order.plant = unname(sapply(order.plant, function(x) as.integer(strsplit(x, "hrs")[[1]][1])))
  order.quant = sapply(quant.files, function(x) strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])])
 
  time.interval.ends = sapply(mapping.files, function(x) tail(strsplit(x, "/")[[1]], 1))
  time.interval.ends = sapply(time.interval.ends, function(x) gsub("\\.txt|hrs", "", x))
  time.interval.ends = unname(as.numeric(sapply(time.interval.ends, function(x) strsplit(x, "_")[[1]][3])))
  time.interval.ends = c(0, time.interval.ends[order(order.plant)])
  
  order.quant = unname(sapply(order.quant, function(x) as.integer(gsub("t|\\.txt", "", x))))
  mapping.files = mapping.files[order(order.plant)]
  quant.files = quant.files[order(order.quant)]
  
  # Read in quantified data
  quant.data = lapply(quant.files, function(x) read.table(x, header=T, sep = "\t"))
  if(length(quant.data) != 0){
    all.quant.data = lapply(1:ncol(quant.data[[1]]), function(x) t(ldply(unname(sapply(quant.data, "[", x)), rbind)))
    names(all.quant.data) = colnames(quant.data[[1]])
  }
  
  # Read in the mapping data between timepoints
  all.mapping.data = lapply(1:length(mapping.files), function(x) {
    data = readLines(textConnection(gsub(":", ",", readLines(mapping.files[x]))))
    data = lapply(data, function(y) gsub(",", "", y))
    data = lapply(data, function(y) as.integer(str_split(y," ")[[1]]))
    data
  })
  
  #############################################################################
  # Extract cell lines
  #############################################################################
  ### Read in first generation and construct the cell line list
  first.mapping.data = all.mapping.data[[1]] 
  cell.lines = lapply(lapply(first.mapping.data, function(x) x[1]), function(x) list(x)) # These are all the original cell lines
  for (ii in 1:length(first.mapping.data)) {
    for (jj in 1:length(cell.lines)) {
      if (first.mapping.data[[ii]][1] %in% cell.lines[[jj]][[1]]) {
        cell.lines[[jj]] = c(cell.lines[[jj]], list(first.mapping.data[[ii]][-1])) # Add daughter cells
        break
      }
    }
  }
  
  ### Now do all but the first
  for (time in 2:length(mapping.files)) {
    mapping.data = all.mapping.data[[time]] 
    # time.interval.end = mapping.files[time]
    
    # Now sort out the cell lines
    for (ii in 1:length(mapping.data)) {
      added = FALSE
      for (jj in 1:length(cell.lines)) {
        # If mother is in last time point
        if (mapping.data[[ii]][1] %in% cell.lines[[jj]][[time]]) {
          # If next timeframe doesn't already exist, add it
          if (length(cell.lines[[jj]]) != time + 1)
            cell.lines[[jj]] = c(cell.lines[[jj]], list(mapping.data[[ii]][-1])) # Add daughter cells to next time frame
          else
            # If it exists, add to it
            cell.lines[[jj]][[time + 1]] = c(cell.lines[[jj]][[time + 1]], mapping.data[[ii]][-1]) # Add daughter cells to next time frame
          added = TRUE
          break
        }
      }
      # If we didn't find it
      if (!added) {
        new.cell.line = list(c(sapply(1:(time - 1), 
                                      function(x) list(NA)), list(mapping.data[[ii]][1])))
        if (length(mapping.data[[ii]]) > 1)
          new.cell.line[[1]] = c(new.cell.line[[1]], 
                                 list(mapping.data[[ii]][-1]))
        cell.lines = c(cell.lines, new.cell.line)
      }
    }
    # Fill up space for the ones that didn't have daughters
    cell.lines[sapply(cell.lines, length) == time] = lapply(cell.lines[sapply(cell.lines, length) == time], function(x) c(x, list(NA)))
  }
  
  # Make it a data frame of lists instead
  cell.lines = t(sapply(cell.lines, function(x) x))
  
  # Make all integer(0) NAs
  for(ii in 1:ncol(cell.lines)) cell.lines[which(!(sapply(cell.lines[,ii], length))), ii] = list(NA)
  
  #############################################################################
  # Extract cell line data
  #############################################################################
  # Find all the expression data, volume, etc. for the lineages
  if(plant.no != 1){
    cell.lines.data =  vector("list", nrow(cell.lines))  # prealloc
    for (ii in 1:nrow(cell.lines)) {
      lineage.data = matrix(NA, nrow = 0, ncol = 1 + ncol(quant.data[[1]]))
      for (jj in 1:no.track.files) {
        cells.in.lineage = which(quant.data[[jj]]$Cell.id %in% cell.lines[[ii, jj]])
        if (length(cells.in.lineage) > 0){
          for (kk in 1:length(cells.in.lineage))
            lineage.data = rbind(lineage.data, c(Time = (jj - 1) * 4, as.vector(quant.data[[jj]][cells.in.lineage[kk], ])))
        } else
          lineage.data = rbind(lineage.data, rep(NA, ncol(lineage.data)))
      }
      cell.lines.data[[ii]] = lineage.data
    }
    
    # Sort cell lines data by number fo members
    cell.lines.members = sapply(cell.lines.data, function(x) length(x[!is.na(x)]))
    cell.lines.data = cell.lines.data[order(-cell.lines.members)]
    cell.lines.data = cell.lines.data[sort(cell.lines.members, decreasing = T) > 0]
  }
  
  #############################################################################
  # Extract mother cells to all cells
  #############################################################################
  # Get the mother cell of all cells
  mothers = lapply(1:length(all.mapping.data), function(x)
    lapply(all.mapping.data[[x]], function(y) {
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
  first.gen = cbind(rep(0, length(all.mapping.data[[1]])),
                    rep(NA, length(all.mapping.data[[1]])),
                    sapply(all.mapping.data[[1]], "[[", 1))
  mothers = rbind(first.gen, mothers)
  colnames(mothers) = c("t", "mother.id", "cell.id")
  
  #############################################################################
  # Get quantified data for all division events
  #############################################################################
  ### Loop through all the time points
  output.data = matrix(NA, ncol = 16, nrow = 0)
  colnames(output.data) = c("cell.id", "lineage.id", "mother.id", "daughter1.id", 
                            "daughter2.id", "age", "t","x","y","z", "vol", "expr", "dist2top", "mother.expr", 
                            "daughter1.expr", "daughter2.expr")
  for(time in 1:length(mapping.files)){
    division.data = all.mapping.data[[time]]
    division.data = division.data[sapply(division.data, length) > 2] # Only keep division events

    top.coords = quant.data[[time]][which(quant.data[[time]]$z == max(quant.data[[time]]$z)), 2:4]
    topcell.coords = top.coords[order(-top.coords$z, top.coords$y**2 + top.coords$x**2), ][1, ]
    
    # for every cell, take out the data FROM all.quant.data that is interesting for that and add it to a row. Append this to the output matrix.
    for(ii in 1:length(division.data)){
      xyzvid = rep(NA, 6)
      this.cell.index = NA
      daughter1.index = NA
      daughter2.index = NA
      daughter1.expr = NA
      daughter2.expr = NA
      
      if(plant.no != 1){
        this.cell.index = all.quant.data$Cell.id[, time] %in% division.data[[ii]][1]
        daughter1.index = all.quant.data$Cell.id[, time] %in% division.data[[ii]][2]
        daughter2.index = all.quant.data$Cell.id[, time] %in% division.data[[ii]][3]
        
        # Does this cell exist in the quantified data?
        if (sum(this.cell.index) != 0) {
          xyzvid = c(
            all.quant.data$x[which(this.cell.index), time],
            all.quant.data$y[which(this.cell.index), time],
            all.quant.data$z[which(this.cell.index), time],
            all.quant.data$Boa.volume[which(this.cell.index), time],
            all.quant.data$Mean.cell.intensity[which(this.cell.index), time],
            sqrt((all.quant.data$x[which(this.cell.index), time] - topcell.coords[1])**2 +
             (all.quant.data$y[which(this.cell.index), time] - topcell.coords[2])**2 +
             (all.quant.data$z[which(this.cell.index), time] - topcell.coords[3])**2)
            )
        }
        
        # Do the daughter cells exist in the quantified data?
        # TODO: Generalise this
        # print(time)
        if ((plant.no == 4 && time > 6) || (plant.no == 15 && time > 3)) {
          daughter1.expr = ifelse(sum(daughter1.index) != 0, all.quant.data$Mean.cell.intensity[which(daughter1.index), time], NA)
          daughter2.expr = ifelse(sum(daughter2.index) != 0, all.quant.data$Mean.cell.intensity[which(daughter2.index), time], NA)
        } else {
          daughter1.expr = ifelse(sum(daughter1.index) != 0, all.quant.data$Mean.cell.intensity[which(daughter1.index), time + 1], NA)
          daughter2.expr = ifelse(sum(daughter2.index) != 0, all.quant.data$Mean.cell.intensity[which(daughter2.index), time + 1], NA)
        }
      }
      
      # Get the identifiers
      lineage.id = which(sapply(1:nrow(cell.lines), function(x) division.data[[ii]][1] %in% cell.lines[[x, time]]))
      this.cell.id = division.data[[ii]][1] # remember, id and index are different
      mother.id = mothers[mothers[, "cell.id"] == this.cell.id & mothers[, "t"] == time - 1, "mother.id"] # Which is the mother id for this cell?
      mother.id = ifelse(length(mother.id) == 0, NA, mother.id)
      mother.expr = output.data[which(output.data[,"t"] == time - 4 & output.data[,"cell.id"] == mother.id), "expr"]
      # mother.expr = output.data[which(output.data[,"t"] == times[time] - times[time - 1] & output.data[,"cell.id"] == mother.id), "expr"]
      mother.expr = ifelse(length(mother.expr) == 0, NA, mother.expr)
      
      ### Find how long the mother has been alive for
      age = 0
      last.mother = mother.id
      for (kk in time:1) {
        # Check the mother. Does she have several offspring?
        reference.id = mothers[mothers[, "cell.id"] == last.mother &
                               mothers[, "t"] == kk - 2, "mother.id"]
        is.not.na = reference.id[!is.na(reference.id)]
        if (length(is.not.na) == 0)
          break
        is.present.previous = mothers[, "t"] == kk - 1
        is.ancestor = mothers[, "mother.id"] == reference.id
        last.mother = reference.id

        # If the mother to the reference has more than one offspring, then a division event has taken place
        no.offspring = length(mothers[is.ancestor & is.present.previous, ])
        if (no.offspring > 1) {
          break
        }
        age = age + 1
      }
      
      # Print output data
      output.data = rbind(output.data, 
                          c(this.cell.id, 
                            lineage.id, 
                            mother.id, 
                            division.data[[ii]][2], 
                            division.data[[ii]][3], 
                            age,
                            (time - 1) * 4,
                            xyzvid, 
                            mother.expr, 
                            daughter1.expr, 
                            daughter2.expr))
    }
  }
  return(output.data)
}

###############################################################################
### Actually run the thing
###############################################################################
plants = c(2, 4, 13, 15, 18)
no.plants = length(plants)
no_cores = detectCores() - 1
cl = makeCluster(no_cores)
results = parLapply(cl, plants, get.plant.data)
stopCluster(cl)

all.data = do.call(rbind, results)

# ################
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
# par(mfrow = c(1, 3))
# no.breaks = 15
# hist(log(sqrt(all.data[, "x"] ** 2 + all.data[, "y"] ** 2)), breaks = no.breaks, main = "ln(rad)")
# hist(((all.data[, "x"] ** 2 + all.data[, "y"] ** 2)), breaks = no.breaks, main = "rad**2")
# ### !!! Radius! Biomodal!
# hist((sqrt(all.data[, "x"] ** 2 + all.data[, "y"] ** 2)), breaks = no.breaks, main = "rad")
# plotDensity(data.frame(x=(sqrt(all.data[, "x"] ** 2 + all.data[, "y"] ** 2))))
# plotDensity(data.frame(x=(sqrt(results[[1]][, "x"] ** 2 + results[[1]][, "y"] ** 2))))
# plotDensity(data.frame(x=(sqrt(results[[2]][, "x"] ** 2 + results[[2]][, "y"] ** 2))))
# plotDensity(data.frame(x=(sqrt(results[[3]][, "x"] ** 2 + results[[3]][, "y"] ** 2))))
# plotDensity(data.frame(x=(sqrt(results[[4]][, "x"] ** 2 + results[[4]][, "y"] ** 2))))
# plotDensity(data.frame(x=(sqrt(results[[5]][, "x"] ** 2 + results[[5]][, "y"] ** 2))))
# 
# # Plot all the quantification data
# # all.data = results[[5]]
# par(mfrow = c(3, 4))
# no.breaks = 1000
# for(ii in 6:ncol(all.data)) if(!all(is.na(all.data[,ii]))) hist(all.data[, ii], breaks = no.breaks, main = colnames(all.data)[ii], freq = T)
# 
# ############
# # Does the time seem to peak at certain intervals? 
# ############
# 
# lapply(1:no.plants, function(x) 
#   {
#   par(mfrow = c(3, 4))
#   for(ii in 6:ncol(results[[x]])) if(!all(is.na(results[[x]][,ii]))) hist(results[[x]][, ii], breaks = no.breaks, main = colnames(results[[x]])[ii], freq = T)
#   Sys.sleep(2)})
# 
# 
# ##################
# ### Do we have double gaussians?
# ##################
# par(mfrow = c(3, 4))
# library(affy)
# plotDensity(data.frame(x=all.data[,"expr"]), main = paste0("expr ", "all"))
# silent = lapply(1:no.plants, function(plant) plotDensity(data.frame(x=results[[plant]][,"expr"]), main = paste0("expr ", plant)))
# plotDensity(data.frame(x=all.data[,"dist2top"]), main = paste0("dist2top ", "all"))
# silent = lapply(1:no.plants, function(plant) plotDensity(data.frame(x=results[[plant]][,"dist2top"]), main = paste0("dist2top ", plant)))
# 
# ##################
# ### Plot all cross-analyses
# ##################
# graphics.off(); plot.new()
# par(mfrow = c(3, 4))
# plots = vector(mode = "list", length = (ncol(all.data) - 5))
# for (ii in 6:ncol(all.data)) {
#   not.ii = (6:ncol(all.data))[(6:ncol(all.data) != ii)] 
#   for (jj in not.ii) {
#     if (length(all.data[!is.na(all.data[, ii]), ii]) != 0 &
#         length(all.data[!is.na(all.data[, jj]), jj]) != 0) {
#       x.lab = colnames(all.data)[jj]
#       y.lab = colnames(all.data)[ii]
#       is.not.na = !is.na(all.data[, y.lab]) & !is.na(all.data[, x.lab])
#       correlation = cor(all.data[is.not.na, x.lab] , all.data[is.not.na, y.lab])
#       title = paste0(y.lab, " vs ", x.lab, "\n cor: ", round(correlation, 2))
#       plot(all.data[, jj], all.data[, ii], main = title, xlab = x.lab, ylab = y.lab)
#     }
#   }
#   plots[[ii - 5]] = recordPlot()
#   plot.new()
# }
# plots[[1]]
# plots[[2]]
# plots[[3]]
# plots[[4]]
# plots[[5]]
# plots[[6]]
# plots[[7]]
# plots[[8]]
# plots[[9]]
# plots[[10]]
# 
# #######################
# # No. divisions per time
# #######################
# 
# tmp.hist = hist(all.data[, "daughter1.expr"] / all.data[, "expr"], breaks = 50, xlim = c(0, 7))
# par(mfrow = c(3, 1))
# hist(all.data[, "daughter1.expr"] / all.data[, "expr"], main = "d1 / m")
# hist(all.data[, "daughter2.expr"] / all.data[, "expr"], main = "d2 / m")
# hist(all.data[, "daughter1.expr"] / all.data[, "daughter2.expr"], main = "d1 / d2", breaks = 50)
