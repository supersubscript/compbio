#!/usr/bin/env Rscript
setwd("~/compbio/src/assignments/cnsa2/code/")
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))

no.weights = 10
no.replicates = 200
lw.s = 2
threshold = 0

classify = function(input, threshold = 0, fct = sum, ...) {
  ifelse(fct(input) < threshold, -1, 1)
}

train.hebb = function(data, no.patterns, ...) {
  targets = apply(data, 2, function(x) classify(x))
  weights = apply(data, 1, function(x) sum(targets * x))
  1/no.weights * weights
}

train.delta = function(data, no.patterns, learning.rate, no.iterations, ...) {
  targets = apply(data, 2, function(x) classify(x)) 
  weights = 0 # Initialise matrix
  for (ii in 1:no.iterations) {
    rand = sample(1:no.patterns, 1)
    output = classify(data[, rand] %*% weights)
    weights = weights + learning.rate * (targets[rand] - output) * data[, rand]
  }
  weights
}
train.delta.time = function(data, no.patterns, learning.rate, ...) {
  targets = apply(data, 2, function(x) classify(x)) 
  weights = 0 # Initialise matrix
  outputs = apply(data, 2, function(x) classify(weights*x)) 
  iters = 1
  while(!all.equal(targets, outputs)) {
    rand = sample(1:no.patterns, 1)
    # output =  classify(data[, rand] %*% weights)
    weights = weights + learning.rate * (targets[rand] - outputs[rand]) * data[, rand]
    outputs = apply(data, 2, function(x) classify(weights*x)) 
    iters = iters + 1
  }
  list(weights=weights, iters=iters)
}

##############################################################################
### Simulate for given number of patterns, with a specified 
### training function and classification function.
##############################################################################
run = function(no.patterns, train.fct, fct = sum, learning.rate = .2, ...) {
  # Generate data and set weights accordingly
  train.data = replicate(no.patterns, sample(c(1, -1), no.weights, replace = T))
  test.data = replicate(no.patterns, sample(c(1, -1), no.weights, replace = T))
  weights = train.fct(train.data, learning.rate)
  
  # Check how many patterns are recalled
  train.s = 0
  test.s = 0
  for (ii in 1:no.patterns) {
    if (classify(train.data[, ii] * weights) == classify(train.data[, ii])) 
      train.s = train.s + 1
    if (classify(test.data[, ii] * weights) == classify(test.data[, ii])) 
      test.s = test.s + 1
  }
  return(list(train = train.s, test = test.s))
}

### RUN SIMULATION
patterns.stored = ceiling(2^(seq(1, 10, 1)))  # How many patterns stored? 
results = lapply(patterns.stored, function(x) lapply(1:no.replicates, function(y) run(x, train.delta, fct = sum)))

### TRAIN DATA
train.means = sapply(results, function(x) mean(sapply(x, function(y) y$train))) / patterns.stored
train.stds = sapply(results, function(x) sd(sapply(x, function(y) y$train)) / sqrt(length(no.replicates))) / patterns.stored

plot(1, cex = 0, xlim = c(1.9, 1025), ylim = c(0, 1), log = "x", xlab = "Patterns stored", ylab = "Fraction correct", bty = "n")
lines(patterns.stored, train.means, lwd = lw.s, col = 1)
suppressWarnings(arrows(patterns.stored, train.means - train.stds, patterns.stored, train.means + train.stds, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 1))

### TEST DATA
test.means = sapply(results, function(x) mean(sapply(x, function(y) y$test)))/patterns.stored
test.stds = sapply(results, function(x) sd(sapply(x, function(y) y$test))/sqrt(length(no.replicates)))/patterns.stored

lines(patterns.stored, test.means, lwd = lw.s, col = 2)
suppressWarnings(arrows(patterns.stored, test.means - test.stds, patterns.stored, test.means + test.stds, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 2))
legend("bottomright", c("Train", "Test"), col = 1:2, lty = 1, lwd = lw.s, inset = c(0.05, 0.05))
