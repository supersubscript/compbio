#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/")
source("aux.R")

plants    = c(2, 4, 13, 15)
no.plants = length(plants)

plant.no = 15
p.idx = all.quant.data[,"plant"] == plant.no
l1.data = lapply(plants, read.data, only.l1 = TRUE)
l2.data = lapply(plants, read.data, only.l2 = TRUE)

l1.quant.data   = lapply(l1.data, "[[", "quant.data") 
l1.mapping.data = lapply(l1.data, "[[", "mapping.data") 
l2.quant.data   = lapply(l2.data, "[[", "quant.data") 
l2.mapping.data = lapply(l2.data, "[[", "mapping.data") 
timepoints      = lapply(l1.data, "[[", "timepoints") 

# Add time and plant.no to quant.data
for(ii in 1:length(l1.quant.data)) for(jj in 1:length(l1.quant.data[[ii]])) {
  l1.quant.data[[ii]][[jj]] = cbind(t = rep(jj, nrow(l1.quant.data[[ii]][[jj]])), l1.quant.data[[ii]][[jj]])
  l1.quant.data[[ii]][[jj]] = cbind(plant = rep(plants[ii], nrow(l1.quant.data[[ii]][[jj]])), l1.quant.data[[ii]][[jj]])  
  l2.quant.data[[ii]][[jj]] = cbind(t = rep(jj, nrow(l2.quant.data[[ii]][[jj]])), l2.quant.data[[ii]][[jj]])
  l2.quant.data[[ii]][[jj]] = cbind(plant = rep(plants[ii], nrow(l2.quant.data[[ii]][[jj]])), l2.quant.data[[ii]][[jj]])
}

all.l1.quant.data   = do.call(rbind, unlist(l1.quant.data, recursive = F))
all.l2.quant.data   = do.call(rbind, unlist(l2.quant.data, recursive = F))

vars = c("plant", "t", "id","x","y","z","vol","expr","d2t")
colnames(all.l1.quant.data) = vars
colnames(all.l2.quant.data) = vars

par(mfcol = c(3, 2), mar = c(4, 4, 4, 4))
for(ii in 7:9) for (jj in (ii:9)[-which(ii:9 == ii)])
  plot(all.l1.quant.data[, ii], all.l1.quant.data[, jj], 
       xlab = colnames(all.l1.quant.data)[ii], ylab = colnames(all.l1.quant.data)[jj])
for(ii in 7:9) for (jj in (ii:9)[-which(ii:9 == ii)])
  plot(all.l2.quant.data[, ii], all.l2.quant.data[, jj], 
       xlab = colnames(all.l2.quant.data)[ii], ylab = colnames(all.l2.quant.data)[jj])

# Plot specific plant
par(mfcol = c(3, 2), mar = c(4, 4, 4, 4))
for(ii in 7:9) for (jj in (ii:9)[-which(ii:9 == ii)])
  plot(all.l1.quant.data[p.idx, ii], all.l1.quant.data[p.idx, jj], 
       xlab = colnames(all.l1.quant.data)[ii], ylab = colnames(all.l1.quant.data)[jj])
for(ii in 7:9) for (jj in (ii:9)[-which(ii:9 == ii)])
  plot(all.l2.quant.data[p.idx, ii], all.l2.quant.data[p.idx, jj], 
       xlab = colnames(all.l2.quant.data)[ii], ylab = colnames(all.l2.quant.data)[jj])


no.breaks = 50
par(mfcol = c(3, 2), mar = c(4, 4, 4, 4))
for(ii in 7:9) hist(all.l1.quant.data[, ii], main = colnames(all.l1.quant.data[ii]), breaks = no.breaks)
for(ii in 7:9) hist(all.l2.quant.data[, ii], main = colnames(all.l2.quant.data[ii]), breaks = no.breaks)

for(ii in 7:9) hist(all.l1.quant.data[p.idx, ii], main = colnames(all.l1.quant.data[ii]), breaks = no.breaks)
for(ii in 7:9) hist(all.l2.quant.data[p.idx, ii], main = colnames(all.l2.quant.data[ii]), breaks = no.breaks)

### Take out central zone
# all.l1.quant.data = all.l1.quant.data[all.l1.quant.data[,"d2t"] < 5, ]
# all.l2.quant.data = all.l2.quant.data[all.l2.quant.data[,"d2t"] < 5, ]

### Take out specific plant
# all.l1.quant.data = all.l1.quant.data[all.l1.quant.data[,"plant"] == 4, ]
# all.l2.quant.data = all.l2.quant.data[all.l2.quant.data[,"plant"] == 4, ]

p1 = plot.violin(all.l1.quant.data, "t", "expr")
p2 = plot.violin(all.l2.quant.data, "t", "expr")
p3 = plot.violin(all.l1.quant.data, "t", "vol")
p4 = plot.violin(all.l2.quant.data, "t", "vol")
grid.arrange(p1, p2, p3, p4, ncol = 2)



### Only l1
par(mfrow=c(2,1)); plot(NA, xlim=c(0,84), ylim=c(0,100))
no.clv3 = ldply(sapply(l1.quant.data, function(x) sapply(x, nrow)), rbind)
timepts = ldply(sapply(plants, get.q.timepoints), rbind)
sapply(1:nrow(timepts), function(x) lines(timepts[x, ], no.clv3[x, ], type="b", col = x))
plot(NA, xlim=c(0,84), ylim=c(0,100))
no.clv3 = ldply(sapply(l2.quant.data, function(x) sapply(x, nrow)), rbind)
timepts = ldply(sapply(plants, get.q.timepoints), rbind)
sapply(1:nrow(timepts), function(x) lines(timepts[x, ], no.clv3[x, ], type="b", col = x))
