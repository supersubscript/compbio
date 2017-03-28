#!/usr/bin/env Rscript
setwd("~/compbio/src/assignments/cnsa2/code/")
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))

no.inputs = 100
no.datap = 10 # N_S
no.replicates = 100
learning.rate = .2
threshold = 0
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

train.delta = function(no.inputs, no.datap, learning.rate, ...) {
  targets = apply(train.data, 2, function(x)
    sum.classify(x))
  # runif(no.inputs*no.datap)
  weights <<- matrix(0, no.inputs, no.datap)
  # print(learning.rate)
  for (ii in 1:no.datap) {
    rand = sample(1:no.datap, 1)
    outputs = sum.classify(train.data[, rand] * weights)
    weights <<-
      weights + learning.rate/2 * (targets[rand] - outputs) * train.data[, rand]
  }
  weights
}
get.stat = function(no.inputs, no.datap, learning.rate) {
  no.inputs <<- no.inputs
  no.datap <<- no.datap
  learning.rate <<- learning.rate
  train.data <<-
    replicate(no.datap, sample(c(1, -1), no.inputs, replace = T))
  test.data  <<-
    replicate(no.datap, sample(c(1, -1), no.inputs, replace = T))
  weights <<- train.delta(no.inputs, no.datap, learning.rate)
  
  train.s = 0
  test.s = 0
  # print(train.data)
  # print(weights)# sapply(1:no.inputs, function(x)
  # print(no.datap)
  # print(no.inputs)
  # print(dim(train.data))
  # print(length(weights))
  for (ii in 1:no.datap) {
    if (sum.classify(train.data[, ii] * weights) == sum.classify(train.data[, ii]))
      train.s = train.s + 1
    if (sum.classify(test.data[, ii] * weights) == sum.classify(test.data[, ii]))
      test.s = test.s + 1
  }
  return(list(train = train.s, test = test.s))
}

# Run! ... and get statistics.
# results = get.stat(10, 10)
# print(results)
rates = seq(0,.2,0.01)
results = lapply(ceiling(2 ^ (seq(1, 10, 1))), function(x)
  lapply(rates, function(z)
    replicate(
      no.replicates, get.stat(x, 10, learning.rate = z)
    )))


# OUTPUT FOR MULTIPLE learning.rate VALUES. no.inputs == 5.
# res = lapply(results, function(x) x[[5]])
# means = sapply(res, function(x)
#   apply(x, 1, function(y)
#     mean(unlist(y))))
# stds  = sapply(res, function(x)
#   apply(x, 1, function(y)
#     sd(unlist(y)) / sqrt(no.replicates)))
# https://www.cs.auckland.ac.nz/~pat/706_98/ln/node148.html

# OUTPUT FOR A SINGLE no.inputs VALUE
results = results[[8]]
means = sapply(results, function(x)
  apply(x, 1, function(y)
    mean(unlist(y))))
stds  = sapply(results, function(x)
  apply(x, 1, function(y)
    sd(unlist(y)) / sqrt(no.replicates)))

# plot(ceiling(2^(seq(1,10,1))), means[1,], type="l", ylim=c(0,10), lwd = lw.s)
# arrows(ceiling(2^(seq(1,10,1))), means[1,]-stds[1,], ceiling(2^(seq(1,10,1))), means[1,]+stds[1,], length=0.05, angle=90, code=3, cex = 2, lwd = lw.s)
# lines(ceiling(2^(seq(1,10,1))), means[2,], col=2, lwd=lw.s)
# arrows(ceiling(2^(seq(1,10,1))), means[2,]-stds[2,], ceiling(2^(seq(1,10,1))), means[2,]+stds[2,], length=0.05, angle=90, code=3, cex = 2, col = 2, lwd = lw.s)
# legend("bottomright", c("Train", "Test"), col=1:2, lwd=lw.s, lty=1, inset=c(.01,.05))
# lines(0:1024, rep(5, 1025), lwd = 2, lty = 2)

# PLOTS FOR A SINGLE no.inputs VALUE
plot(rates, means[1,], type="l", ylim=c(0,10), lwd = lw.s)
arrows(rates, means[1,]-stds[1,], rates, means[1,]+stds[1,], length=0.05, angle=90, code=3, cex = 2, lwd = lw.s)
lines(rates, means[2,], col=2, lwd=lw.s)
arrows(rates, means[2,]-stds[2,], rates, means[2,]+stds[2,], length=0.05, angle=90, code=3, cex = 2, col = 2, lwd = lw.s)
legend("bottomright", c("Train", "Test"), col=1:2, lwd=lw.s, lty=1, inset=c(.01,.05))
# lines(0:1024, rep(5, 1025), lwd = 2, lty = 2)


