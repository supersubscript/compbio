#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/src/assignments/sba3/code/")
# Do the track analysis
library(plyr)
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 2
par(lwd=lw.s, cex = 1.5, ps = 15)
library(prodlim)

scores =
read.csv("../figures/DYR_ECOLI_e3_n2_m40_theta_0.3_pc_weight_0.5_DIScoresCompared.csv",
         sep = ",", header = F)


TPs = function(nIncluded) {
    TP = sum(scores[1:nIncluded, 4] < 5.0)
  TP / nIncluded
}

fract = sapply(1:nrow(scores), TPs)
par(mfrow = c(1, 1))

plot(fract, type = "l", ylab = "TP / (TP + FP)", xlab = "TP + FP")



