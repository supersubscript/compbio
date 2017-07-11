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

# all.l1.quant.data = all.l1.quant.data[all.l1.quant.data[,"d2t"] < 5, ]
# all.l2.quant.data = all.l2.quant.data[all.l2.quant.data[,"d2t"] < 5, ]

# all.l1.quant.data = all.l1.quant.data[all.l1.quant.data[,"plant"] == 4, ]
# all.l2.quant.data = all.l2.quant.data[all.l2.quant.data[,"plant"] == 4, ]

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



# t = all.quant.data[,"t"]; expr = all.quant.data[,"expr"]; there = nna(t) & nna(expr)
# p = ggplot(data.frame(t=t[there], expr=expr[there]), aes(factor(t), expr)); p + geom_violin()
# t = all.quant.data[,"t"]; d2t = all.quant.data[,"d2t"]; there = nna(t) & nna(d2t)
# p = ggplot(data.frame(t=t[there], d2t=d2t[there]), aes(factor(t), d2t)); p + geom_violin()


### Only l1
# no.clv3 = ldply(sapply(quant.data, function(x) sapply(x, nrow)), rbind)
# timepts = ldply(sapply(plants, get.q.timepoints), rbind)
# 
# par(mfrow=c(1,1)); plot(NA, xlim=c(0,84), ylim=c(0,110))
# sapply(1:nrow(timepts), function(x) lines(timepts[x, ], no.clv3[x, ], type="b", col = x))
