#!/usr/bin/env Rscript
setwd("~/compbio/src/assignments/cnsa2/code/")
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))

no.inputs = 1000
no.datap = 10 # N_S
no.replicates = 100
lw.s = 2
train.data = replicate(no.datap, sample(c(1, -1), no.inputs, replace = T))
test.data  = replicate(no.datap, sample(c(1, -1), no.inputs, replace = T))
weights = runif(no.inputs, min = 0, max = 1)

sum.classify = function(input) {
  ifelse(sum(input) < 0,-1, 1)
}
prod.classify = function(input) {
  ifelse(prod(input) < 0,-1, 1)
}

train.hebb = function(weights, ...) {
  targets = apply(train.data, 2, function(x)
    sum.classify(x))
  weights = sapply(1:no.inputs, function(x)
    sum(targets * train.data[x,]))
  1 / no.inputs * weights
}

train.delta = function(weights, learning.rate, ...) {
  targets = apply(train.data, 2, function(x)
    sum.classify(x))
  for(ii in 1:no.datap){
    outputs = apply(train.data, 2, function(x)
      sum.classify(weights*x))
    weights = weights + learning.rate * (outputs - targets)*train.data[, ii]
  }
}

get.stat = function(no.inputs, no.datap) {
  no.inputs <<- no.inputs
  no.datap <<- no.datap
  train.data <<-
    replicate(no.datap, sample(c(1, -1), no.inputs, replace = T))
  test.data  <<-
    replicate(no.datap, sample(c(1, -1), no.inputs, replace = T))
  weights <<- train.hebb(weights)
  
  train.s = 0
  test.s = 0
  for (ii in 1:no.datap) {
    if (sum.classify(train.data[, ii] * weights) == sum.classify(train.data[, ii]))
      train.s = train.s + 1
    if (sum.classify(test.data[, ii] * weights) == sum.classify(test.data[, ii]))
      test.s = test.s + 1
  }
  return(list(train = train.s, test = test.s))
}

# Run! ... and get statistics.
results = lapply(ceiling(2 ^ (seq(1, 10, 1))), function(x)
  replicate(no.replicates, get.stat(x, 10)))
means = sapply(results, function(x)
  apply(x, 1, function(y)
    mean(unlist(y)))) 
stds  = sapply(results, function(x)
  apply(x, 1, function(y)
    sd(unlist(y)) / sqrt(no.replicates))) 

plot(ceiling(2^(seq(1,10,1))), means[1,], type="l", ylim=c(0,10), lwd = lw.s)
arrows(ceiling(2^(seq(1,10,1))), means[1,]-stds[1,], ceiling(2^(seq(1,10,1))), means[1,]+stds[1,], length=0.05, angle=90, code=3, cex = 2, lwd = lw.s)
lines(ceiling(2^(seq(1,10,1))), means[2,], col=2, lwd=lw.s)
arrows(ceiling(2^(seq(1,10,1))), means[2,]-stds[2,], ceiling(2^(seq(1,10,1))), means[2,]+stds[2,], length=0.05, angle=90, code=3, cex = 2, col = 2, lwd = lw.s)
legend("bottomright", c("Train", "Test"), col=1:2, lwd=lw.s, lty=1, inset=c(.01,.05))
lines(0:1024, rep(5, 1025), lwd = 2, lty = 2)




