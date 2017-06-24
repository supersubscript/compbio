#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/")
source("aux.R", local = TRUE)
source("read_data.R", local = TRUE)
source("process_data.R", local = TRUE)

# 4 and 15 lacks a timepoint in quant. 18 lacks a tracking timepoint.
plant.no = 2

# get.plant.data = function(plant.no){
  # Load the necessary packages
  packages = c("reshape", "scales", "stringr", "plyr", "parallel", "limma", "affy")
  try(lapply(packages, library, character.only = TRUE), silent = TRUE)
  
  # Read in data
  all.data     = read.data(plant.no, map.quant = TRUE) # Remember to check this
  mapping.data = all.data$mapping.data
  quant.data   = all.data$quant.data
  timepoints   = all.data$timepoints
  
  # Extract cell lines
  cell.lines = get.cell.lines(mapping.data, sort.order = "no.members")
  sublines   = get.sublines(mapping.data, timepoints)
  lineages   = get.lineages(sublines)
  
  lineage.data = lapply(lineages, 
                        get.lineage.data, 
                        quant.data = quant.data, 
                        timepoints = timepoints)
  lineage.data = lineage.data[lapply(lineage.data, length) > 0] # Only save ones with data
  
  lineage.sublines.data = lapply(lineages, 
                                 get.lineage.sublines.data, 
                                 quant.data = quant.data, 
                                 timepoints = timepoints)
  
  
  
  
  # subline.data = get.sublines.data(
  #   sublines,
  #   mapping.data,
  #   quant.data,
  #   timepoints,
  #   min.no.members = 1,
  #   param          = "Mean.cell.intensity",
  #   no.z           = 1,
  #   no.xy          = 1
  # )
  # 
  # subline.data = subline.data[order.na(subline.data)]
  
  
  # Find all the expression data, volume, etc. for the lineages
  # if (plant.no != 1) {
  #   cell.lines.data = get.cell.line.data(
  #     cell.lines,
  #     quant.data,
  #     timepoints,
  #     min.no.members = 1,
  #     param          = "Mean.cell.intensity",
  #     no.z           = 1,
  #     no.xy          = 1
  #   )
  # }
  
  # Extract mother cells to all cells
  mothers = get.mothers(mapping.data, timepoints)
  
  # Retrieve table with all division events and their corresponding information
  all.data = get.division.events(
    mapping.data,
    quant.data,
    cell.lines,
    param = "Mean.cell.intensity",
    no.z  = 1,
    no.xy = 1,
    plant.no)
  # return(all.data)
# }

# plants    = c(1, 2, 4, 13, 15, 18)
# plants    = 2
# no.plants = length(plants)
# no_cores  = detectCores() - 1
# cl        = makeCluster(no_cores)
# silent    = clusterExport(cl, ls())
# results   = parLapply(cl, plants, get.plant.data)
# stopCluster(cl)

# all.data = do.call(rbind, results)

################
### In big enough intervals, cellular divisions follow a normal distribution with the logged radius. Does the probability to divide
### decay exponentially with the radius? Also radius squared?
################

no.breaks = 30
par(mfrow = c(1, 3))
hist(log(sqrt(all.data[, "x"] ** 2 + all.data[, "y"] ** 2)), breaks = no.breaks, main = "ln(rad)")
hist(((all.data[, "x"] ** 2 + all.data[, "y"] ** 2)),        breaks = no.breaks, main = "rad**2")
hist((sqrt(all.data[, "x"] ** 2 + all.data[, "y"] ** 2)),    breaks = no.breaks, main = "rad")

# Plot all the quantification data
par(mfrow = c(4, 3))
no.breaks = 30
for(ii in 6:ncol(all.data)) if(!all(is.na(all.data[,ii]))) hist(all.data[, ii], breaks = no.breaks, main = colnames(all.data)[ii], freq = F)
for(ii in 6:ncol(all.data)) if(!all(is.na(all.data[,ii]))) plotDensity(data.frame(x=all.data[, ii]), main = colnames(all.data)[ii])
# plot(density(all.data[, "dist2top"], na.rm = T, bw = .2))

graphics.off(); plot.new()
par(mfrow = c(3, 3), mar = c(2, 2, 2, 2))
plots       = vector(mode = "list", length = (ncol(all.data) - 5))
interesting = (1:ncol(all.data))[-c(1:5, 8:10)] # Remove xyz and indices
counter     = 1
for (ii in interesting) {
  not.ii = interesting[interesting != ii]
  for (jj in not.ii) {
    if (length(all.data[!is.na(all.data[, ii]), ii]) != 0 &
        length(all.data[!is.na(all.data[, jj]), jj]) != 0) {
      x.lab       = colnames(all.data)[jj]
      y.lab       = colnames(all.data)[ii]
      is.not.na   = !is.na(all.data[, y.lab]) & !is.na(all.data[, x.lab])
      correlation = cor(all.data[is.not.na, x.lab] , all.data[is.not.na, y.lab])
      title       = paste0(y.lab, " vs ", x.lab, "\n cor: ", round(correlation, 2))
      plot(all.data[, jj], all.data[, ii], main = title, xlab = x.lab, ylab = y.lab)
    }
  }
  plots[[counter]] = recordPlot()
  plot.new()
  counter = counter + 1
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

##################################
# plot.12(cell.lines.data)


plot(subline.data[[1]][,"Time"], subline.data[[1]][,"Mean.cell.intensity"])

subline.data = subline.data[lapply(subline.data, function(x) sum(!is.na(unlist(x)))) > 0]
subline.data = subline.data[order(unlist(lapply(subline.data, function(x)
  mean(unlist(x[, 8]), na.rm = T))))]

plot.12(subline.data, batch = 1)

plot.lineage(lineage.data[["351"]])


# df = data.frame(ID = 1:12, do.call(rbind, lapply(subline.data[1:12], function(x) unlist(x[,"Mean.cell.intensity"]))))
# t  = data.frame(ID = 1:12, do.call(rbind, lapply(subline.data[1:12], function(x) unlist(x[,"Time"]))))


# library(traj)
# s1 = step1measures(df, t, ID = TRUE)
# s2 = step2factors(s1, discard = TRUE)
# s3 = step3clusters(s2, nclusters = 4)
