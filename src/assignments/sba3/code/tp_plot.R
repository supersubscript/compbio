#!/usr/bin/env Rscript
setwd("/local/data/public/hpa22/compbio/src/assignments/sba3/code/")
# Do the track analysis
library(plyr)
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 2
par(lwd=lw.s, cex = 1.5, ps = 15)
library(prodlim)
library(lattice)
library(gridExtra)

real.files = list.files("../figures", pattern = "crystal_contacts.csv", full.names = T)
real.files = grep("DYR", real.files, value = T)
pred.files = list.files("../figures", pattern = "prediction_contacts.csv", full.names = T)
pred.files = grep("DYR", pred.files, value = T)

order.files  = as.numeric(sapply(sapply(real.files, function(x) strsplit(x,"_theta_")[[1]][2]), function(x) strsplit(x, "_")[[1]][1]))
order.files2 = as.numeric(sapply(sapply(real.files, function(x) strsplit(x,"pc_weight_")[[1]][2]), function(x) strsplit(x, "_")[[1]][1]))

real.files = real.files[order(order.files, order.files2)]
pred.files = pred.files[order(order.files, order.files2)]

a = expand.grid(unique(order.files), unique(order.files2))
a = a[order(a[,1], a[,2]), ]
name = paste(a[,1], a[,2], sep = "_")

results = sapply(1:length(real.files), function(x){
  pred = read.csv(pred.files[x], sep = ",")
  real = read.csv(real.files[x], sep = ",")
  # real = read.csv("../figures/DYR_ECOLI_e3_n2_m40_theta_0.3_pc_weight_0.5_crystal_contacts.csv", sep = ",")
  # pred = read.csv("../figures/DYR_ECOLI_e3_n2_m40_theta_0.3_pc_weight_0.5_prediction_contacts.csv", sep = ",")
  TP = sum(!is.na(row.match(pred, real))) # Is in both
  FP = sum(is.na(row.match(pred, real)))  # Is in pred but not in real
  FN = sum(is.na(row.match(real, pred)))  # Is in real but not in pred
    
  sens = TP / (TP + FN)
  prec = TP / (TP + FP)
  list(TP = TP, FP = FP, FN = FN, sens=sens, prec=prec)    
})

colnames(results) = name

# x is PC-weight, y is 
prec = results["prec",]
sens = results["sens",]
prec.m = matrix(unlist(prec), ncol=11, nrow=11, byrow=T)
sens.m = matrix(unlist(sens), ncol=11, nrow=11, byrow=T)
par(mfrow=c(2,1))
plot1 = levelplot(prec.m,  col.regions = viridis(21), xlab = "Pseudocount weight", ylab = expression(theta), main = "Precision", panel = panel.levelplot.raster)
plot2 = levelplot(sens.m,  col.regions = viridis(21), xlab = "Pseudocount weight", ylab = expression(theta), main = "Sensitivity")
grid.arrange(plot1, plot2, ncol = 1)



