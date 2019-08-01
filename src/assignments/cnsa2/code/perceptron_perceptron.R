#!/usr/bin/env Rscript
setwd("~/compbio/src/assignments/cnsa2/code/")
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))
bias = TRUE
no.weights = 2 + ifelse(bias, 1, 0)
no.replicates = 100
lw.s = 2

classify = function(input, fct, ...) {
  ifelse(fct(input) < 0, -1, 1)
}

train.hebb = function(data, no.patterns, fct, ...) {
  targets = apply(data, 2, function(x) classify(x, fct=fct))
  weights = apply(data, 1, function(x) sum(targets * x))
  1/no.weights * weights
}

train.delta = function(data, no.patterns, learning.rate, no.iterations, fct, delta.fct, ...) {
  targets = apply(data, 2, function(x) classify(x, fct=fct)) # Find the correct answers
  weights = runif(no.weights, min=-0.01, max = 0.01) # Initialise weights
  for (ii in 1:no.iterations) {
    rand = sample(1:no.patterns, 1)
    output = classify(data[, rand] %*% weights, fct=fct)
    weights = weights + learning.rate * (targets[rand] - output) * delta.fct(data[, rand])
  }
  weights
}

train.delta.time = function(data, no.patterns, learning.rate, fct, delta.fct, ...) {
  targets = apply(data, 2, function(x) classify(x, fct=fct)) 
  weights = runif(no.weights, max = 0.01, min = -0.01)#0 # Initialise matrix
  outputs = apply(data, 2, function(x) classify(weights*x, fct=fct)) 
  iters = 1
  while (!all(targets == outputs)) {
    rand = sample(1:no.patterns, 1)
    weights = weights + learning.rate * (targets[rand] - classify(weights %*% delta.fct(data[,rand]), fct)) * delta.fct(data[, rand])
    outputs = apply(data, 2, function(x) classify(weights %*% x, fct=fct)) 
    # print(weights)
    iters = iters + 1
  }
  iters
}

train.perceptron = function(data, no.patterns, learning.rate, no.iterations, fct, ...) {
  targets = apply(data, 2, function(x) classify(x,  fct=fct)) 
  weights = runif(no.weights, min=-0.01, max = 0.01) # Initialise matrix
  for (ii in 1:no.iterations) {
    rand = sample(1:no.patterns, 1)
    output = classify(data[, rand] %*% weights, fct=fct)
    weights = weights + learning.rate / 2 * (targets[rand] - output) * data[, rand]
  }
  weights
}
train.perceptron.time = function(data, no.patterns, learning.rate, no.iterations, fct, ...) {
  targets = apply(data, 2, function(x) classify(x,  fct=fct)) 
  weights = runif(no.weights, min=-0.01, max = 0.01) # Initialise matrix
  outputs = apply(data, 2, function(x) classify(weights*x, fct=fct)) 
  iters = 1
  while (!all(targets == outputs)) {
    rand = sample(1:no.patterns, 1)
    weights = weights + learning.rate / 2 * (targets[rand] - outputs[rand]) * data[, rand]
    
    outputs = apply(data, 2, function(x) classify(weights %*% x, fct=fct)) 
    iters = iters + 1
  }
  iters
}
  
##############################################################################
### Simulate for given number of patterns, with a specified 
### training function and classification function.
##############################################################################
run = function(no.patterns, train.fct, fct, learning.rate, delta.fct = function(x) x, ...) {
  # Generate data and set weights accordingly
  train.data = replicate(no.patterns, sample(c(1, -1), no.weights, replace = T))
  train.data[1,] = rep(1, no.patterns) #rbind(train.data, rep(1, no.patterns))
  test.data = replicate(no.patterns, sample(c(1, -1), no.weights, replace = T))
  test.data[1,] = rep(1, no.patterns) #rbind(test.data, rep(1, no.patterns))
  weights = train.fct(data=train.data, no.patterns=no.patterns, no.iterations = no.patterns*10, learning.rate=learning.rate, fct=fct, delta.fct=delta.fct)
  
  # Check how many patterns are recalled
  train.s = 0
  test.s = 0
  for (ii in 1:no.patterns) {
    if (classify(train.data[, ii] * weights, fct=fct) == classify(train.data[, ii], fct=fct)) 
      train.s = train.s + 1
    if (classify(test.data[, ii] * weights, fct=fct) == classify(test.data[, ii], fct=fct)) 
      test.s = test.s + 1
  }
  return(list(train = train.s, test = test.s))
}

run.time = function(no.patterns, train.fct, fct, learning.rate, delta.fct = function(x) x,  ...) {
  # Generate data and set weights accordingly
  train.data = replicate(no.patterns, sample(c(1, -1), no.weights, replace = TRUE))
  iters = train.fct(data=train.data, no.patterns=no.patterns, learning.rate=learning.rate, fct=fct, delta.fct=delta.fct)
  iters 
}

##############################################################################
### RUN SIMULATIONS
##############################################################################
patterns.stored = ceiling(2^(seq(1, 10, 1)))  # How many patterns stored?
pat.time = c(5, 10, 15, 20)
learning = 0.01
erf =  function(x) 2 * pnorm(x * sqrt(2)) - 1
delta.rule = function(x) 100*x #erf #1/(1 - exp(-x))

### SUM NORMAL
results.hebb = lapply(patterns.stored, function(x)
  lapply(1:no.replicates, function(y)
    run(
      no.patterns = x,
      train.fct = train.hebb,
      fct = sum,
      learning.rate = learning
    )))
results.delta = lapply(patterns.stored, function(x)
  lapply(1:no.replicates, function(y)
    run(
      no.patterns = x,
      train.fct = train.delta,
      fct = sum,
      learning.rate = learning,
      delta.fct = delta.rule # x # log(x + 10)
    )))
results.perceptron = lapply(patterns.stored, function(x)
  lapply(1:no.replicates, function(y)
    run(
      no.patterns = x,
      train.fct = train.perceptron,
      fct = sum,
      learning.rate = learning
    )))

### PROD NORMAL
results.prod.hebb = lapply(patterns.stored, function(x)
  lapply(1:no.replicates, function(y)
    run(
      no.patterns = x,
      train.fct = train.hebb,
      fct = prod,
      learning.rate = learning
    )))
results.prod.delta = lapply(patterns.stored, function(x)
  lapply(1:no.replicates, function(y)
    run(
      no.patterns = x,
      train.fct = train.delta,
      fct = prod,
      learning.rate = learning,
      delta.fct = delta.rule #function(x) log(x + 10)
    )))
results.prod.perceptron = lapply(patterns.stored, function(x)
  lapply(1:no.replicates, function(y)
    run(
      no.patterns = x,
      train.fct = train.perceptron,
      fct = prod,
      learning.rate = learning
    )))

### SUM TIME
time.results.delta = lapply(pat.time, function(x)
  lapply(1:no.replicates, function(y)
    run.time(
      no.patterns = x,
      train.fct = train.delta.time,
      fct = sum,
      learning.rate = learning,
      delta.fct = delta.rule #function(x) x #log(x + 10)
    )))
time.results.perceptron = lapply(pat.time, function(x)
  lapply(1:no.replicates, function(y)
    run.time(
      no.patterns = x,
      train.fct = train.perceptron.time,
      fct = sum,
      learning.rate = learning
    )))

# ### PROD DELTA TIME
# time.results.prod.delta = lapply(pat.time, function(x)
#   lapply(1:no.replicates, function(y)
#     run.time(
#       no.patterns = x,
#       train.fct = train.delta.time,
#       fct = prod,
#       learning.rate = learning,
#       delta.fct = delta.rule #function(x) 1/(1 + exp(-x))
#     )))

# if(exists("../data/learning_raesults.RData")) {
#   load("../data/learning_results.RData")
# }

# save(results.hebb, results.delta, results.perceptron, time.results.delta, time.results.perceptron, results.prod.hebb, results.prod.delta, results.prod.perceptron,  file = "../data/learning_results.RData")
# time.results.prod.delta
##############################################################################
### CALCULATE MEANS + STDDEV
##############################################################################
# Sum
hebb.train.means = sapply(results.hebb, function(x) mean(sapply(x, function(y) y$train))) / patterns.stored
hebb.train.stds = sapply(results.hebb, function(x) sd(sapply(x, function(y) y$train)) / sqrt(length(no.replicates))) / patterns.stored
hebb.test.means = sapply(results.hebb, function(x) mean(sapply(x, function(y) y$test)))/patterns.stored
hebb.test.stds = sapply(results.hebb, function(x) sd(sapply(x, function(y) y$test))/sqrt(length(no.replicates)))/patterns.stored
delta.train.means = sapply(results.delta, function(x) mean(sapply(x, function(y) y$train))) / patterns.stored
delta.train.stds = sapply(results.delta, function(x) sd(sapply(x, function(y) y$train)) / sqrt(length(no.replicates))) / patterns.stored
delta.test.means = sapply(results.delta, function(x) mean(sapply(x, function(y) y$test)))/patterns.stored
delta.test.stds = sapply(results.delta, function(x) sd(sapply(x, function(y) y$test))/sqrt(length(no.replicates)))/patterns.stored
perceptron.train.means = sapply(results.perceptron, function(x) mean(sapply(x, function(y) y$train))) / patterns.stored
perceptron.train.stds = sapply(results.perceptron, function(x) sd(sapply(x, function(y) y$train)) / sqrt(length(no.replicates))) / patterns.stored
perceptron.test.means = sapply(results.perceptron, function(x) mean(sapply(x, function(y) y$test)))/patterns.stored
perceptron.test.stds = sapply(results.perceptron, function(x) sd(sapply(x, function(y) y$test))/sqrt(length(no.replicates)))/patterns.stored

# Sum time
delta.train.means.time = sapply(time.results.delta, function(x) mean(sapply(x, function(y) y))) / pat.time
delta.train.stds.time = sapply(time.results.delta, function(x) sd(sapply(x, function(y) y)) / sqrt(length(no.replicates))) / pat.time
perceptron.train.means.time = sapply(time.results.perceptron, function(x) mean(sapply(x, function(y) y))) / pat.time
perceptron.train.stds.time = sapply(time.results.perceptron, function(x) sd(sapply(x, function(y) y)) / sqrt(length(no.replicates))) / pat.time

# Prod
hebb.train.prod.means = sapply(results.prod.hebb, function(x) mean(sapply(x, function(y) y$train))) / patterns.stored
hebb.train.prod.stds = sapply(results.prod.hebb, function(x) sd(sapply(x, function(y) y$train)) / sqrt(length(no.replicates))) / patterns.stored
hebb.test.prod.means = sapply(results.prod.hebb, function(x) mean(sapply(x, function(y) y$test)))/patterns.stored
hebb.test.prod.stds = sapply(results.prod.hebb, function(x) sd(sapply(x, function(y) y$test))/sqrt(length(no.replicates)))/patterns.stored
delta.train.prod.means = sapply(results.prod.delta, function(x) mean(sapply(x, function(y) y$train))) / patterns.stored
delta.train.prod.stds = sapply(results.prod.delta, function(x) sd(sapply(x, function(y) y$train)) / sqrt(length(no.replicates))) / patterns.stored
delta.test.prod.means = sapply(results.prod.delta, function(x) mean(sapply(x, function(y) y$test)))/patterns.stored
delta.test.prod.stds = sapply(results.prod.delta, function(x) sd(sapply(x, function(y) y$test))/sqrt(length(no.replicates)))/patterns.stored
perceptron.train.prod.means = sapply(results.prod.perceptron, function(x) mean(sapply(x, function(y) y$train))) / patterns.stored
perceptron.train.prod.stds = sapply(results.prod.perceptron, function(x) sd(sapply(x, function(y) y$train)) / sqrt(length(no.replicates))) / patterns.stored
perceptron.test.prod.means = sapply(results.prod.perceptron, function(x) mean(sapply(x, function(y) y$test)))/patterns.stored
perceptron.test.prod.stds = sapply(results.prod.perceptron, function(x) sd(sapply(x, function(y) y$test))/sqrt(length(no.replicates)))/patterns.stored

# Prod time
# delta.train.prod.means.time = sapply(time.results.prod.delta, function(x) mean(sapply(x, function(y) y))) / pat.time
# delta.train.prod.stds.time = sapply(time.results.prod.delta, function(x) sd(sapply(x, function(y) y)) / sqrt(length(no.replicates))) / pat.time

##############################################################################
### PLOT
##############################################################################
par(mfcol = c(3, 2), oma = c(10, 5, 4, 2) + 0.0, mai = c(.0, 0, .0, .0))
### SUM HEBB
plot(1, cex = 0, xlim = c(min(patterns.stored), max(patterns.stored)), ylim = c(0, 1), log = "x", xlab = "", ylab = "", bty = "n", xaxt = "n")
lines(patterns.stored, hebb.train.means, lwd = lw.s, col = 1)
suppressWarnings(arrows(patterns.stored, hebb.train.means - hebb.train.stds, patterns.stored, hebb.train.means + hebb.train.stds, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 1))
lines(patterns.stored, hebb.test.means, lwd = lw.s, col = 2)
suppressWarnings(arrows(patterns.stored, hebb.test.means - hebb.test.stds, patterns.stored, hebb.test.means + hebb.test.stds, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 2))
# legend("bottomright", c("Train", "Test"), col = 1:2, lty = 1, lwd = lw.s, inset = c(0.05, 0.05))

### SUM DELTA
plot(1, cex = 0, xlim = c(min(patterns.stored), max(patterns.stored)), ylim = c(0, 1), log = "x", xlab = "", ylab = "", bty = "n", xaxt = "n")
lines(patterns.stored, delta.train.means, lwd = lw.s, col = 1)
suppressWarnings(arrows(patterns.stored, delta.train.means - delta.train.stds, patterns.stored, delta.train.means + delta.train.stds, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 1))
lines(patterns.stored, delta.test.means, lwd = lw.s, col = 2)
suppressWarnings(arrows(patterns.stored, delta.test.means - delta.test.stds, patterns.stored, delta.test.means + delta.test.stds, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 2))
mtext("Fraction correct", side=2, line=3,las=3, col="black", cex = 1)
# legend("bottomright", c("Train", "Test"), col = 1:2, lty = 1, lwd = lw.s, inset = c(0.05, 0.05))

### SUM PERCEPTRON
plot(1, cex = 0, xlim = c(min(patterns.stored), max(patterns.stored)), ylim = c(0, 1), log = "x", xlab = "", ylab = "", bty = "n")
lines(patterns.stored, perceptron.train.means, lwd = lw.s, col = 1)
suppressWarnings(arrows(patterns.stored, perceptron.train.means - perceptron.train.stds, patterns.stored, perceptron.train.means + perceptron.train.stds, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 1))
lines(patterns.stored, perceptron.test.means, lwd = lw.s, col = 2)
suppressWarnings(arrows(patterns.stored, perceptron.test.means - perceptron.test.stds, patterns.stored, perceptron.test.means + perceptron.test.stds, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 2))
# mtext("Patterns stored", side=1, line=3,las=1, col="black", cex = 1)
legend("bottomright", c("Train", "Test"), col = 1:2, lty = 1, lwd = lw.s, inset = c(0.05, 0.05))

### PROD HEBB
plot(1, cex = 0, xlim = c(min(patterns.stored), max(patterns.stored)), ylim = c(0, 1), log = "x", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n")
lines(patterns.stored, hebb.train.prod.means, lwd = lw.s, col = 1)
suppressWarnings(arrows(patterns.stored, hebb.train.prod.means - hebb.train.prod.stds, patterns.stored, hebb.train.prod.means + hebb.train.prod.stds, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 1))
lines(patterns.stored, hebb.test.prod.means, lwd = lw.s, col = 2)
suppressWarnings(arrows(patterns.stored, hebb.test.prod.means - hebb.test.prod.stds, patterns.stored, hebb.test.prod.means + hebb.test.prod.stds, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 2))
# legend("bottomright", c("Train", "Test"), col = 1:2, lty = 1, lwd = lw.s, inset = c(0.05, 0.05))

### PROD DELTA
plot(1, cex = 0, xlim = c(min(patterns.stored), max(patterns.stored)), ylim = c(0, 1), log = "x", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n")
lines(patterns.stored, delta.train.prod.means, lwd = lw.s, col = 1)
suppressWarnings(arrows(patterns.stored, delta.train.prod.means - delta.train.prod.stds, patterns.stored, delta.train.prod.means + delta.train.prod.stds, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 1))
lines(patterns.stored, delta.test.prod.means, lwd = lw.s, col = 2)
suppressWarnings(arrows(patterns.stored, delta.test.prod.means - delta.test.prod.stds, patterns.stored, delta.test.prod.means + delta.test.prod.stds, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 2))
# mtext("Fraction correct", side=2, line=3,las=3, col="black", cex = 1)
# legend("bottomright", c("Train", "Test"), col = 1:2, lty = 1, lwd = lw.s, inset = c(0.05, 0.05))

### PROD PERCEPTRON
plot(1, cex = 0, xlim = c(min(patterns.stored), max(patterns.stored)), ylim = c(0, 1), log = "x", xlab = "", ylab = "", bty = "n", yaxt = "n")
lines(patterns.stored, perceptron.train.prod.means, lwd = lw.s, col = 1)
suppressWarnings(arrows(patterns.stored, perceptron.train.prod.means - perceptron.train.prod.stds, patterns.stored, perceptron.train.prod.means + perceptron.train.prod.stds, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 1))
lines(patterns.stored, perceptron.test.prod.means, lwd = lw.s, col = 2)
suppressWarnings(arrows(patterns.stored, perceptron.test.prod.means - perceptron.test.prod.stds, patterns.stored, perceptron.test.prod.means + perceptron.test.prod.stds, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 2))
mtext("Patterns stored", side=1, line=3,las=1, col="black", cex = 1, at = 1.5)
legend("bottomright", c("Train", "Test"), col = 1:2, lty = 1, lwd = lw.s, inset = c(0.05, 0.05))

##############################################################################
### PLOT TIME
##############################################################################
# par(mfrow = c(1, 1), oma = c(10, 5, 4, 2) + 0.0, mai = c(.0, .6, .0, .5))
# par(mfrow=c(1,1))
# # SUM DELTA
# plot(1, cex = 0, xlim = c(min(pat.time), max(pat.time)), ylim = c(0, 10), log = "", xlab = "", ylab = "", bty = "n")
# lines(pat.time, delta.train.means.time, lwd = lw.s, col = 1)
# suppressWarnings(arrows(pat.time, delta.train.means.time - delta.train.stds.time, pat.time, delta.train.means.time + delta.train.stds.time, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 1))
# mtext("Learning time", side=2, line=3,las=3, col="black", cex = 1)
# 
# # SUM PERCEPTRON
# lines(pat.time, perceptron.train.means.time, lwd = lw.s, col = 2)
# suppressWarnings(arrows(pat.time, perceptron.train.means.time - perceptron.train.stds.time, pat.time, perceptron.train.means.time + perceptron.train.stds.time, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 2))

# # PROD DELTA
# lines(pat.time, delta.train.prod.means.time, lwd = lw.s, col = 3)
# suppressWarnings(arrows(pat.time, delta.train.prod.means.time - delta.train.prod.stds.time, pat.time, delta.train.prod.means.time + delta.train.prod.stds.time, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 3))
# mtext("Learning time", side=2, line=3,las=3, col="black", cex = 3)
# legend("bottomright", c("SDelta", "SPerceptron", "PDelta"), col = 1:3, lty = 1, lwd = lw.s, inset = c(0.05, 0.05))