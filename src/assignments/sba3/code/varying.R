#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/src/assignments/sba3/code/")
# Do the track analysis
library(plyr)
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 2
par(lwd=lw.s, cex = 1.5, ps = 15)
library(prodlim)
real.files   = list.files("../figures", pattern = "DIScoresC", full.names = T)
real.files   = grep("DYR", real.files, value = T)
order.files  = as.numeric(sapply(sapply(real.files, function(x) strsplit(x,"_theta_")[[1]][2]), function(x) strsplit(x, "_")[[1]][1]))
order.files2 = as.numeric(sapply(sapply(real.files, function(x) strsplit(x,"pc_weight_")[[1]][2]), function(x) strsplit(x, "_")[[1]][1]))
real.files   = real.files[order(order.files, order.files2)]

a = expand.grid(unique(order.files), unique(order.files2))
a = a[order(a[,1], a[,2]), ]
name = paste(a[,1], a[,2], sep = "_")


theta.0.3.files = which(a[,1] == 0.3)[c(2,4,6,8,10)]
pc.0.5.files    = which(a[,2] == 0.5)[c(2,4,6,8,10)]

TPs = function(nIncluded, scores) {
  TP = sum(scores[1:nIncluded, 4] < 5.0)
  TP / nIncluded
}

results = list()
for(ii in real.files[theta.0.3.files]){
  scores = read.csv(ii, sep = ",", header = F)
  fract = sapply(1:nrow(scores), function(x) TPs(x, scores))
  results = append(results, list(fract))
}

par(mfrow = c(2,1))
plot(NA, ylim= c(0,1), xlim = c(0,1000), xlab = "TP + FP", ylab = "TP / (TP + FP)", main = expression(" "*theta*" = 0.3, varying pc weight"))
for(ii in 1:length(results)) lines(results[[ii]], col = ii)
legend("topright", legend=a[theta.0.3.files, 2], col = 1:5, lty = 1, lwd = lw.s, inset = c(0.02,0.02), bg = "white")

#########################################################################
results = list()
for(ii in real.files[pc.0.5.files]){
  scores = read.csv(ii, sep = ",", header = F)
  fract = sapply(1:nrow(scores), function(x) TPs(x, scores))
  results = append(results, list(fract))
}
plot(NA, ylim= c(0,1), xlim = c(0,1000), xlab = "TP + FP", ylab = "TP / (TP + FP)", main = expression("pc weight = 0.5, varying "*theta*" "))
for(ii in 1:length(results)) lines(results[[ii]], col = ii)
legend("topright", legend=a[pc.0.5.files, 1], col = 1:5, lty = 1, lwd = lw.s, inset = c(0.02,0.02), bg = "white")
