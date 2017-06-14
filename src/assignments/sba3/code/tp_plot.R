#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/src/assignments/sba3/code/")
# Do the track analysis
library(plyr)
library(RColorBrewer)
library(scales)
library(affy)
library(limma)
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 2
par(lwd=lw.s, cex = 1.5, ps = 15)

real = read.csv("../data/crystal_contacts.csv", sep = ",")
pred = read.csv("../data/prediction_contacts.csv", sep = ",")

TP = 0
FP = 0
for(ii in 1:nrow(pred)){
  added = FALSE
  for(jj in 1:nrow(real)){
    if(all(unlist(pred[ii, ]) == unlist(real[jj,]))) {
      TP = TP +1
      added = TRUE
      break
    }
  }
  if(!added)
    FP = FP + 1
}

