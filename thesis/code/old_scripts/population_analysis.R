#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/")
source("aux.R")

plants          = c(2, 4, 13, 15)
no.plants       = length(plants)
quant.missing   = list(c(-1), c(24), c(-1), c(12))
mapping.missing = list(c(-1), c(-1), c(-1), c(-1))

plant.no = 4
no_cores  = detectCores() - 1
cl        = makeCluster(no_cores)
silent    = clusterExport(cl, ls())
all.data = parLapply(cl, 1:no.plants, function(x) 
  read.data(plants[x], quant.missing = quant.missing[[x]], map.quant = TRUE, only.l1 = FALSE))
l1.data = parLapply(cl, 1:no.plants, function(x) 
  read.data(plants[x], quant.missing = quant.missing[[x]], map.quant = TRUE, only.l1 = TRUE))
l2.data = parLapply(cl, 1:no.plants, function(x) 
  read.data(plants[x], quant.missing = quant.missing[[x]], map.quant = TRUE, only.l2 = TRUE))
silent = stopCluster(cl)

all.quant.data   = lapply(all.data, "[[", "quant.data") 
all.mapping.data = lapply(all.data, "[[", "mapping.data") 
l1.quant.data   = lapply(l1.data, "[[", "quant.data") 
l1.mapping.data = lapply(l1.data, "[[", "mapping.data") 
l2.quant.data   = lapply(l2.data, "[[", "quant.data") 
l2.mapping.data = lapply(l2.data, "[[", "mapping.data") 
timepoints      = lapply(l1.data, "[[", "timepoints") 
darktimes       = c(10:18, 34:42, 58:66, 82:84)

# Add time and plant.no to quant.data
for(ii in 1:length(l1.quant.data)) {
  for(jj in 1:length(l1.quant.data[[ii]])) {
    all.quant.data[[ii]][[jj]] = cbind(t = rep(jj, nrow(all.quant.data[[ii]][[jj]])), all.quant.data[[ii]][[jj]])
    all.quant.data[[ii]][[jj]] = cbind(plant = rep(plants[ii], nrow(all.quant.data[[ii]][[jj]])), all.quant.data[[ii]][[jj]]) 
    l1.quant.data[[ii]][[jj]] = cbind(t = rep(jj, nrow(l1.quant.data[[ii]][[jj]])), l1.quant.data[[ii]][[jj]])
    l1.quant.data[[ii]][[jj]] = cbind(plant = rep(plants[ii], nrow(l1.quant.data[[ii]][[jj]])), l1.quant.data[[ii]][[jj]])  
    l2.quant.data[[ii]][[jj]] = cbind(t = rep(jj, nrow(l2.quant.data[[ii]][[jj]])), l2.quant.data[[ii]][[jj]])
    l2.quant.data[[ii]][[jj]] = cbind(plant = rep(plants[ii], nrow(l2.quant.data[[ii]][[jj]])), l2.quant.data[[ii]][[jj]])
  }
}
all.all.quant.data = do.call(rbind, unlist(l1.quant.data, recursive = F))
all.l1.quant.data = do.call(rbind, unlist(l1.quant.data, recursive = F))
all.l2.quant.data = do.call(rbind, unlist(l2.quant.data, recursive = F))

vars = c("plant", "t", "id","x","y","z","vol","expr","d2t")
colnames(all.all.quant.data) = vars
colnames(all.l1.quant.data) = vars
colnames(all.l2.quant.data) = vars

cz.rad = 2
all.alla.quant.data = all.all.quant.data[which(all.all.quant.data[,"d2t"] < cz.rad), ]
all.l1.quant.data = all.l1.quant.data[which(all.l1.quant.data[,"d2t"] < cz.rad), ]
all.l2.quant.data = all.l2.quant.data[which(all.l2.quant.data[,"d2t"] < cz.rad), ]
all.p.idx = all.all.quant.data[,"plant"] == plant.no
l1.p.idx  = all.l1.quant.data[,"plant"] == plant.no
l2.p.idx  = all.l2.quant.data[,"plant"] == plant.no


par(mfcol = c(3, 2), mar = c(4, 4, 4, 4))
for(ii in 7:9) for (jj in (ii:9)[-which(ii:9 == ii)])
  plot(all.all.quant.data[, jj], all.all.quant.data[, ii], 
       xlab = colnames(all.all.quant.data)[jj], ylab = colnames(all.all.quant.data)[ii], main = "all")
par(mfcol = c(3, 2), mar = c(4, 4, 4, 4))
for(ii in 7:9) for (jj in (ii:9)[-which(ii:9 == ii)])
  plot(all.l1.quant.data[, jj], all.l1.quant.data[, ii], 
       xlab = colnames(all.l1.quant.data)[jj], ylab = colnames(all.l1.quant.data)[ii], main = "L1")
for(ii in 7:9) for (jj in (ii:9)[-which(ii:9 == ii)])
  plot(all.l2.quant.data[, jj], all.l2.quant.data[, ii], 
       xlab = colnames(all.l2.quant.data)[jj], ylab = colnames(all.l2.quant.data)[ii], main = "L2")

par(mfcol = c(3, 2), mar = c(4, 4, 4, 4))
for(ii in 7:9) for (jj in (ii:9)[-which(ii:9 == ii)])
  plot(all.l1.quant.data[l1.p.idx, jj], all.l1.quant.data[l1.p.idx, ii], 
       xlab = colnames(all.l1.quant.data)[jj], ylab = colnames(all.l1.quant.data)[ii])
for(ii in 7:9) for (jj in (ii:9)[-which(ii:9 == ii)])
  plot(all.l2.quant.data[l2.p.idx, jj], all.l2.quant.data[l2.p.idx, ii], 
       xlab = colnames(all.l2.quant.data)[jj], ylab = colnames(all.l2.quant.data)[ii])

no.breaks = 30
par(mfcol = c(3, 2), mar = c(4, 4, 4, 4))
for(ii in 7:9) hist(all.l1.quant.data[, ii], main = colnames(all.l1.quant.data[ii]), breaks = no.breaks)
for(ii in 7:9) hist(all.l2.quant.data[, ii], main = colnames(all.l2.quant.data[ii]), breaks = no.breaks)

for(ii in 7:9) hist(all.l1.quant.data[l1.p.idx, ii], main = colnames(all.l1.quant.data[ii]), breaks = no.breaks)
for(ii in 7:9) hist(all.l2.quant.data[l2.p.idx, ii], main = colnames(all.l2.quant.data[ii]), breaks = no.breaks)

par(mfrow=c(2,1))
t = all.l1.quant.data[,"t"]; expr = all.l1.quant.data[,"expr"]; there = nna(t) & nna(expr)
p = ggplot(data.frame(t=t[there], expr=expr[there]), aes(factor(t), expr)); p1 = p + geom_violin()
t = all.l2.quant.data[,"t"]; expr = all.l2.quant.data[,"expr"]; there = nna(t) & nna(expr)
p = ggplot(data.frame(t=t[there], expr=expr[there]), aes(factor(t), expr)); p2 = p + geom_violin()
grid.arrange(p1, p2 , ncol = 1)

t = all.l1.quant.data[,"t"]; vol = all.l1.quant.data[,"vol"]; there = nna(t) & nna(vol)
p = ggplot(data.frame(t=t[there], vol=vol[there]), aes(factor(t), vol)); p1 = p + geom_violin()
t = all.l2.quant.data[,"t"]; vol = all.l2.quant.data[,"vol"]; there = nna(t) & nna(vol)
p = ggplot(data.frame(t=t[there], vol=vol[there]), aes(factor(t), vol)); p2 = p + geom_violin()
grid.arrange(p1, p2 , ncol = 1)
