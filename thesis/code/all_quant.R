#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/")
# Do the track analysis
library(plyr)
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 2
par(lwd=lw.s, cex = 1.5, ps = 15)
    
plants = c(2, 4, 13, 15, 18)

# quant.files = list.files(paste0("../data/msb_plants/plant", plant.no, "/results"), full.names = TRUE)
all.quant.data =  lapply(plants, function(plant.no) {
  unlist(lapply(list.files(paste0("../data/msb_plants/plant", plant.no, "/results"), full.names = TRUE), function(y)
         lapply(y, function(x) read.table(x, header=T, sep = "\t"))), recursive = F)
  })

# all.quant.data.merged = rbind.fill.matrix(t(t(all.quant.data)))

no.clv3.nuclei = sapply(lapply(all.quant.data, function(x) sapply(x, nrow)), function(y) t(matrix(y)))
no.clv3.nuclei = rbind.fill.matrix(no.clv3.nuclei)

mean.no.clv3.nuclei = apply(no.clv3.nuclei, 2, mean, na.rm = T)
sd.no.clv3.nuclei = apply(no.clv3.nuclei, 2, sd, na.rm = T)

clv3.mean = lapply(all.quant.data, function(x) sapply(x, function(y) mean(y[[6]])))
clv3.sd = lapply(all.quant.data, function(x) sapply(x, function(y) sd(y[[6]])))

vol.mean = lapply(all.quant.data, function(x) sapply(x, function(y) mean(y[[5]])))
vol.sd = lapply(all.quant.data, function(x) sapply(x, function(y) sd(y[[5]])))

################################
### PLOT
################################
par(mfrow = c(3,2))
apply(no.clv3.nuclei, 1, plot, type = "b", main = "no.clv3.nuclei")

par(mfrow = c(3,1))
plot(mean.no.clv3.nuclei, type = "b", main = "mean.no.clv3.nuclei")
arrows(1:22,
  mean.no.clv3.nuclei - sd.no.clv3.nuclei / sqrt(apply(no.clv3.nuclei, 2, function(x) sum(!is.na(x)))),
  1:22, mean.no.clv3.nuclei + sd.no.clv3.nuclei / sqrt(apply(no.clv3.nuclei, 2, function(x) sum(!is.na(x)))),
  length = 0.05, angle = 90, code = 3)
plot(sd.no.clv3.nuclei, type = "b", main = "sd.no.clv3.nuclei")
plot(sd.no.clv3.nuclei / mean.no.clv3.nuclei, type = "b", main = "sd.no.clv3.nuclei / mean.no.clv3.nuclei")

par(mfrow = c(3, 2))
lapply(1:length(clv3.mean), function(x) plot(clv3.mean[[x]], type = "b", main = paste0("clv3.mean ", plants[x])))

par(mfrow = c(3, 2))
lapply(1:length(clv3.sd), function(x) plot(clv3.sd[[x]], type = "b", main = paste0("clv3.sd ", plants[x])))

par(mfrow = c(3, 2))
lapply(1:length(vol.mean), function(x) plot(vol.mean[[x]], type = "b", main = paste0("vol.mean ", plants[x])))

par(mfrow = c(3, 2))
lapply(1:length(vol.sd), function(x) plot(vol.sd[[x]], type = "b", main = paste0("vol.sd ", plants[x])))
