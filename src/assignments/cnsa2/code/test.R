#!/usr/bin/env Rscript
setwd("~/compbio/src/assignments/cnsa2/code/")
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))

# no.inputs = 100
# no.datap = 10 # N_S
no.replicates = 100
# learning.rate = .2
lw.s = 2
# threshold = 0
# train.data = replicate(no.datap, sample(c(1, -1), no.inputs, replace = T))
# test.data  = replicate(no.datap, sample(c(1, -1), no.inputs, replace = T))
# weights = runif(no.inputs, min = 0, max = 1)

sum.classify = function(input, threshold) {
  ifelse(sum(input) < threshold, -1, 1)
}
prod.classify = function(input, threshold) {
  ifelse(prod(input) < threshold, -1, 1)
}

train.hebb = function(no.inputs, no.datap, train.data, test.data, threshold, ...) {
  targets = apply(train.data, 2, function(x)
    sum.classify(x, threshold))
  weights = sapply(1:no.inputs, function(x)
    sum(targets * train.data[x,]))
  1 / no.inputs * weights
}

train.delta = function(no.inputs, no.datap, train.data, test.data, threshold, learning.rate,  ...) {
  targets = apply(train.data, 2, function(x)
    sum.classify(x, threshold))
  a = matrix(runif(no.inputs * no.datap), no.inputs, no.datap)
  for (ii in 1:no.datap) {
    rand = sample(1:no.datap, 1)
    output = sum.classify(train.data[, rand] * a, threshold)
    a =
      a + learning.rate * (targets[rand] - output) * train.data[, rand]
  }
  a
}

train.perceptron = function(no.inputs, no.datap, train.data, test.data, threshold, learning.rate, ...) {
  targets = apply(train.data, 2, function(x)
    sum.classify(x, threshold))
  a = matrix(0, no.inputs, no.datap)
  for (ii in 1:no.datap) {
    rand = sample(1:no.datap, 1)
    output = sum.classify(train.data[, rand] * a, threshold)
    a =
      a + learning.rate / 2 * (targets[rand] - output) * train.data[, rand]
    threshold =
      threshold - learning.rate / 2 * (targets[rand] - output)
  }
  list(a, threshold)
}

get.stat = function(no.inputs, no.datap, learning.rate, ...) {
  threshold = 0
  train.data = replicate(no.datap, sample(c(1, -1), no.inputs, replace = T))
  test.data  = replicate(no.datap, sample(c(1, -1), no.inputs, replace = T))
  tr = train.delta(no.inputs, no.datap, train.data, test.data, threshold, learning.rate)
  
  w = tr[[1]]
  if(length(tr) > 1)
    threshold = tr[[2]]
  train.s = 0
  test.s = 0
  for (ii in 1:no.datap) {
    if (sum.classify(train.data[, ii] * w, threshold) == sum.classify(train.data[, ii], threshold))
      train.s = train.s + 1
    if (sum.classify(test.data[, ii] * w, threshold) == sum.classify(test.data[, ii], threshold))
      test.s = test.s + 1
  }
  return(list(train = train.s, test = test.s))
}

# Run! ... and get statistics.
# results = get.stat(10, 10)
# print(results)

# What parameters should we try?
rates = seq(0,1,0.2)
nodes = ceiling(2 ^ (seq(1, 10, 1)))
# nodes = 100

results = lapply(nodes, function(x)
  lapply(rates, function(z)
    replicate(no.replicates, get.stat(no.inputs=x, no.datap=10, learning.rate=z))))

# OUTPUT FOR MULTIPLE no.inputs VALUES. learning.rate == 5.
res = lapply(results, function(x) x[[2]])
means = sapply(res, function(x)
  apply(x, 1, function(y)
    mean(unlist(y))))
stds  = sapply(res, function(x)
  apply(x, 1, function(y)
    sd(unlist(y)) / sqrt(no.replicates)))
plot(nodes, means[1,], type="l", ylim=c(0,10), lwd = lw.s)
arrows(nodes, means[1,]-stds[1,], nodes, means[1,]+stds[1,], length=0.05, angle=90, code=3, cex = 2, lwd = lw.s)
lines(nodes, means[2,], col=2, lwd=lw.s)
arrows(nodes, means[2,]-stds[2,], nodes, means[2,]+stds[2,], length=0.05, angle=90, code=3, cex = 2, col = 2, lwd = lw.s)
legend("bottomright", c("Train", "Test"), col=1:2, lwd=lw.s, lty=1, inset=c(.01,.05))
lines(0:1024, rep(5, 1025), lwd = 2, lty = 2)

# OUTPUT FOR A SINGLE no.inputs VALUE
# results = results[[1]]
# means = sapply(results, function(x)
#   apply(x, 1, function(y)
#     mean(unlist(y))))
# stds  = sapply(results, function(x)
#   apply(x, 1, function(y)
#     sd(unlist(y)) / sqrt(no.replicates)))
# 
# 
# # PLOTS FOR A SINGLE no.inputs VALUE
# plot(rates, means[1,], type="l", ylim=c(0,10), lwd = lw.s)
# arrows(rates, means[1,]-stds[1,], rates, means[1,]+stds[1,], length=0.05, angle=90, code=3, cex = 2, lwd = lw.s)
# lines(rates, means[2,], col=2, lwd=lw.s)
# arrows(rates, means[2,]-stds[2,], rates, means[2,]+stds[2,], length=0.05, angle=90, code=3, cex = 2, col = 2, lwd = lw.s)
# legend("bottomright", c("Train", "Test"), col=1:2, lwd=lw.s, lty=1, inset=c(.01,.05))
# lines(0:1024, rep(5, 1025), lwd = 2, lty = 2)

