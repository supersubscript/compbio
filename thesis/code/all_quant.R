#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/")
# Do the track analysis
library(plyr)
library(RColorBrewer)
library(scales)
library(affy)
library(limma)
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 2
par(lwd=lw.s, cex = 1.5, ps = 15)
    
plants = c(2, 4, 13, 15, 18)
no.plants = length(plants)

quant.files = lapply(plants, function(plant.no) list.files(paste0("../data/clv3_complete/plant", plant.no, "/Results"), full.names = TRUE, pattern=".txt"))
order.quant = lapply(quant.files, function(y) sapply(y, function(x) strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])]))
order.quant = lapply(order.quant, function(y) unname(sapply(y, function(x) as.integer(gsub("t|\\.txt", "", x)))))
quant.files = lapply(1:length(quant.files), function(x) quant.files[[x]][order(order.quant[[x]])])

# Sort 
timepoints = lapply(quant.files, function(file) sapply(file, function(times) tail(strsplit(times, "/")[[1]], 1)))
timepoints = lapply(timepoints, function(file) unname(sapply(file, function(times) unname(as.integer(gsub("t|\\.txt", "", times))))))

quant.files = lapply(1:length(quant.files), function(x) quant.files[[x]][order(timepoints[[x]])])

all.quant.data =  lapply(1:no.plants, function(plant.no) {
  unlist(lapply(quant.files[[plant.no]], function(y)
         lapply(y, function(x) read.table(x, header=T, sep = "\t"))), recursive = F)
  })

series = lapply(all.quant.data, function(y) lapply(1:ncol(y[[1]]), function(x) unname(sapply(y, "[", x))))
all.all.quant.data = list()
for (ii in 1:5) all.all.quant.data[[ii]] = lapply(1:6, function(x) t(ldply(unname(series[[ii]][[x]]), rbind)))

no.clv3.nuclei = sapply(lapply(all.quant.data, function(x) sapply(x, nrow)), function(y) t(matrix(y)))
no.clv3.nuclei = rbind.fill.matrix(no.clv3.nuclei)

timepoints = rbind.fill.matrix(sapply(timepoints, function(x) t(matrix(x))))

for(plant in 1:nrow(timepoints)){
  for(time in 2:ncol(timepoints)){
    if(!is.na(timepoints[plant, time]) && !is.na(timepoints[plant, time-1])){
      if(timepoints[plant, time] != timepoints[plant, time - 1] + 4){
        timepoints[plant, (time + 1):(ncol(timepoints))] = timepoints[plant, (time):(ncol(timepoints) - 1)]
        no.clv3.nuclei[plant, (time + 1):(ncol(timepoints))] = no.clv3.nuclei[plant, (time):(ncol(timepoints) - 1)]
        timepoints[plant, time] = NA
        no.clv3.nuclei[plant, time] = NA
      }
    }
  }
}

mean.no.clv3.nuclei = apply(no.clv3.nuclei, 2, mean, na.rm=T)
sd.no.clv3.nuclei = apply(no.clv3.nuclei, 2, sd, na.rm = T)

clv3.mean = lapply(all.quant.data, function(x) sapply(x, function(y) mean(y[[6]])))
clv3.sd = lapply(all.quant.data, function(x) sapply(x, function(y) sd(y[[6]])))

vol.mean = lapply(all.quant.data, function(x) sapply(x, function(y) mean(y[[5]])))
vol.sd = lapply(all.quant.data, function(x) sapply(x, function(y) sd(y[[5]])))

################################
### PLOT
################################
par(mfrow = c(3,2))
for(ii in 1:nrow(no.clv3.nuclei)) plot(timepoints[ii,], no.clv3.nuclei[ii,], type = "b", main = "no.clv3.nuclei", xlim=c(0,84))
plot(seq(0,84,4), mean.no.clv3.nuclei, type = "b", main = "mean.no.clv3.nuclei")


par(mfrow = c(3,1))
plot(seq(0,84,4), mean.no.clv3.nuclei, type = "b", main = "mean.no.clv3.nuclei")
arrows(seq(0,84,4),
  mean.no.clv3.nuclei - sd.no.clv3.nuclei / sqrt(apply(no.clv3.nuclei, 2, function(x) sum(!is.na(x)))),
  seq(0,84,4), mean.no.clv3.nuclei + sd.no.clv3.nuclei / sqrt(apply(no.clv3.nuclei, 2, function(x) sum(!is.na(x)))),
  length = 0.05, angle = 90, code = 3)
plot(seq(0,84,4), sd.no.clv3.nuclei, type = "b", main = "sd.no.clv3.nuclei")
plot(seq(0,84,4), sd.no.clv3.nuclei / mean.no.clv3.nuclei, type = "b", main = "sd.no.clv3.nuclei / mean.no.clv3.nuclei")

par(mfrow = c(3, 2))
lapply(1:length(clv3.mean), function(x) plot(clv3.mean[[x]], type = "b", main = paste0("clv3.mean ", plants[x])))

par(mfrow = c(3, 2))
lapply(1:length(clv3.sd), function(x) plot(clv3.sd[[x]], type = "b", main = paste0("clv3.sd ", plants[x])))

par(mfrow = c(3, 2))
lapply(1:length(vol.mean), function(x) plot(vol.mean[[x]], type = "b", main = paste0("vol.mean ", plants[x])))

par(mfrow = c(3, 2))
lapply(1:length(vol.sd), function(x) plot(vol.sd[[x]], type = "b", main = paste0("vol.sd ", plants[x])))

################################
################################

plant.index = 1
data = all.all.quant.data[[plant.index]]
names(data) = c("id", "x","y","z","vol","int")
topcells = apply(data$z, 2, function(x) which(unlist(x) == max(x[!is.na(x)]))[1])
# apply(data$z, 2, function(x) which(unlist(x) == max(x[!is.na(x)]))[1])
top.x = sapply(1:length(topcells), function(x) data$x[topcells[x], x])
top.y = sapply(1:length(topcells), function(x) data$y[topcells[x], x])
top.z = sapply(1:length(topcells), function(x) data$z[topcells[x], x])

dist2top = sqrt(sweep(data$x, 2, top.x)**2 + sweep(data$y, 2, top.y)**2 + sweep(data$z, 2, top.z)**2)
dat = data.frame(x=dist2top, y = data$int)
plot(c(dist2top), c(data$int))
timepoint = 14
order.z = order(c(dist2top[, timepoint]))
plot(dist2top[, timepoint][order.z],  data$int[, timepoint][order.z]/(data$int[, timepoint][order.z][1]), xlim = c(0, 10), ylim = c(0,10))

################################
################################
par(mfrow=c(1,1))
all.dist2top = lapply(all.all.quant.data, function(data) {
  names(data) = c("id", "x","y","z","vol","int")
  topcells = apply(data$z, 2, function(x) which(unlist(x) == max(x[!is.na(x)]))[1])
  # apply(data$z, 2, function(x) which(unlist(x) == max(x[!is.na(x)]))[1])
  top.x = sapply(1:length(topcells), function(x) data$x[topcells[x], x])
  top.y = sapply(1:length(topcells), function(x) data$y[topcells[x], x])
  top.z = sapply(1:length(topcells), function(x) data$z[topcells[x], x])
  
  sqrt(sweep(data$x, 2, top.x)**2 + sweep(data$y, 2, top.y)**2 + sweep(data$z, 2, top.z)**2)
})

all.z = data.matrix(plyr::ldply(lapply(all.all.quant.data, "[[", 4), rbind))
all.vol = data.matrix(plyr::ldply(lapply(all.all.quant.data, "[[", 5), rbind))
all.int = data.matrix(plyr::ldply(lapply(all.all.quant.data, "[[", 6), rbind))
all.dist2top = data.matrix(plyr::ldply(all.dist2top, rbind))

plot(c(all.dist2top), c(all.int))

df = data.frame(dist2top=c(all.dist2top), int=c(all.int), z = c(all.z), vol=c(all.vol))
plot(df, cex = .5, col = alpha(1, .1))
  
# hist(c(all.dist2top)/c(all.int), breaks = 1000)
# plotDensity(data.frame(x=c(all.dist2top), y=c(all.int)))
plot(c(all.dist2top), c(all.int), col = alpha(1,.1))
cor.test(c(all.dist2top), c(all.int), use="complete")

# plot(c(all.dist2top[which(all.dist2top < 2)]), c(all.int[which(all.dist2top < 2)]))








###### Observation: all vol to int plots seem to be composed of two functions. One linear-ish and one sigmoidal. Cells about to divide?
# plot(dat)
plot(c(data$int), c(dist2top))
cor.test(c(dist2top), c(data$int), use="complete")
plot(c(data$vol), c(data$int)) # Sigmoidal?
cor.test(c(data$vol), c(data$int), use="complete")
plot(c(dist2top), c(data$vol))
cor.test(c(dist2top), c(data$vol), use="complete")

plotDensities(data$int)
plotDensities(dist2top)

par(mfrow=c(1,1))
for(ii in 1:ncol(data$vol)) {
  # plot(c(data$z[,ii]), c(data$int[,ii]), xlim = c(0,10), ylim = c(0,160))
  # plot(c(data$vol[,ii]), c(data$int[,ii]), xlim = c(0,2), ylim = c(0,160))
  # plotDensity(data.frame(x=data$z[,ii]))
  # plotDensity(data.frame(x=data$vol[,ii]))
  Sys.sleep(.5)
}
