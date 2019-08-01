#!/usr/bin/env Rscript
allow.inverse = FALSE
threshold = 0
pattern.length = 100
no.iterations = 5000
technical.replicates = 1
no.nodes = pattern.length + ifelse(allow.inverse, 1, 0)
no.patterns.stored = seq(5, 25)
distortion.fracs = seq(0, 0, by = 0.1)
fraction.ones = seq(0.5, 0.5, by = 0.1)
weights.removed = seq(0, 0, by = 0.1)

################################################
# Functions we use
################################################
calc.energy = function(input) {
  energy = 0
  for (ii in 1:ncol(weights)) {
    for (jj in 1:nrow(weights)) {
      energy = energy - 0.5 * weights[ii, jj] * input[ii] * input[jj]
    }
    energy = energy + threshold * input[ii]
  }
  return(energy)
}

hebb.learn.patterns = function(patterns,
                               weights,as
                               this.no.patterns.stored) {
  for (ii in 1:ncol(weights))
    for (jj in 1:nrow(weights)) {
      weights[ii, jj] <<-
        1.0 / this.no.patterns.stored * sum(patterns[ii, ] * patterns[jj, ])
    }
  diag(weights) <<- 0
}

get.h = function(ii, jj, nu, prev.weights, patterns) {
  sum = 0
  for (ll in (1:ncol(prev.weights))[-c(ii, jj)])
    sum = sum + prev.weights[ii, ll] * patterns[ll, nu]
  return(sum)
}

storkey.learn.patterns = function(patterns,
                                  weights,
                                  this.no.patterns.stored) {
  prev.w = matrix(0, no.nodes, no.nodes)
  for (kk in 1:ncol(patterns)) {
    for (ii in 1:ncol(weights))
      for (jj in 1:nrow(weights)) {
        weights[ii, jj] <<- prev.w + 1.0 / this.no.patterns.stored * (
          patterns[ii, kk] * patterns[jj, kk] - patterns[ii, kk] * get.h(ii, jj, kk, prev.w, patterns) - patterns[jj, kk] *
            get.h(jj, ii, kk, prev.w, patterns)
        )
      }
    prev.w = weights
  }
  diag(weights) <<- 0
}



update.async.random = function(input) {
  rand = sample(1:no.nodes, 1)
  input[rand] = ifelse(sum(weights[rand,] * input) > threshold, 1, -1)
  return(input)
}

update.neuron = function(input, index) {
  input[index] = ifelse(sum(weights[index,] * input) > threshold, 1, -1)
  return(input)
}

output.pattern = function(output, pattern) {
  all(output == pattern)
}

################################################
# Actually run stuff
################################################
patterns.stored <- vector("list", length(fraction.ones))
names(patterns.stored) = fraction.ones
for (ii in 1:length(patterns.stored)) {
  patterns.stored[[ii]] <- vector("list", length(weights.removed))
  names(patterns.stored[[ii]]) = weights.removed
  for (jj in 1:length(patterns.stored[[ii]])) {
    patterns.stored[[ii]][[jj]] <-
      vector("list", length(distortion.fracs))
    names(patterns.stored[[ii]][[jj]]) = distortion.fracs
    for (kk in 1:length(patterns.stored[[ii]][[jj]])) {
      patterns.stored[[ii]][[jj]][[kk]] <-
        vector("list", length(no.patterns.stored))
      names(patterns.stored[[ii]][[jj]][[kk]]) = no.patterns.stored
    }
  }
}

# Run over the different numbers of fraction ones.
data = list()
for (iter in 1:20){
for (this.no.patterns.stored in no.patterns.stored) {
  cat("Simulating with ",
      this.no.patterns.stored,
      "patterns stored.\n")
  for (this.fraction.ones in fraction.ones) {
    patterns = replicate(this.no.patterns.stored,
                         sample(c(
                           rep(1, pattern.length * this.fraction.ones),
                           rep(-1, pattern.length * (1 - this.fraction.ones))
                         ), pattern.length, replace = FALSE))
    nodes = sample(c(1, -1), no.nodes, replace = TRUE)
    
    for (this.weights.removed in weights.removed) {
      # Set the weights
      weights = matrix(runif(pattern.length ** 2,  min = -1, max = 1),
                       pattern.length,
                       pattern.length)
      diag(weights) = 0
      hebb.learn.patterns(patterns, weights, this.no.patterns.stored)
      to.remove = sample(1:no.nodes, this.weights.removed * no.nodes, replace = FALSE)
      weights[to.remove, to.remove] = 0
      
      for (this.distortion.frac in distortion.fracs) {
        patterns.found = 0 # This is what we measure.
        
        # Simulate over the patterns
        for (this.pattern in 1:ncol(patterns)) {
          # Do some replicates for data
          for (rep in 1:technical.replicates) {
            # Distort the given fraction of bits.
            output = patterns[, this.pattern]
            flipbits = sample(1:pattern.length,
                              this.distortion.frac * pattern.length,
                              replace = FALSE)
            if (length(flipbits) > 0)
              output[flipbits] = -1 * output[flipbits]
            
            # See if convergence.
            for (ii in 1:no.iterations) {
              output = update.async.random(output)
            }
            if (any(apply(patterns, 2, function(x)
              output.pattern(output, x)))) {
              patterns.found = patterns.found + 1
            }
          }
        }
        patterns.stored[[as.character(this.fraction.ones)]][[as.character(this.weights.removed)]][[as.character(this.distortion.frac)]][[as.character(this.no.patterns.stored)]] = patterns.found / technical.replicates
      }
    }
  }
}
  data = append(data, patterns.stored)
}

avg = apply(sapply(data, function(x) unlist(x)), 1, mean)
serr = apply(sapply(data, function(x) unlist(x)), 1, function(x) sd(x, na.rm = TRUE) / sqrt(length(x)))
plot(no.patterns.stored, avg, pch = 19, xlab = "Patterns stored", ylab = "Patterns retrieved")
lines(no.patterns.stored, col = "red")
arrows(no.patterns.stored, avg-serr, no.patterns.stored, avg+serr, length=0.05, angle=90, code=3, cex = 2)

# print(calc.energy(output))
# cat("patterns found:\t", patterns.found, "\n")
# image(matrix(input, nrow = sqrt(pattern.length - ifelse(allow.inverse, 1, 0)), ncol = sqrt(pattern.length - ifelse(allow.inverse, 1, 0))))
# image(matrix(output, nrow = sqrt(pattern.length - ifelse(allow.inverse, 1, 0)), ncol = sqrt(pattern.length - ifelse(allow.inverse, 1, 0))))
# output = update.neuron(output, ii %% pattern.length + 1)
