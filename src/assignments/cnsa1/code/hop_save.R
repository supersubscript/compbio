#!/usr/bin/env Rscript
allow.inverse = FALSE
threshold = 0
pattern.length = 100
no.iterations = 5000
technical.replicates = 1
no.nodes = pattern.length
no.patterns.stored = seq(10, 10)
distortion.fracs = seq(0, 0, by = 0.1)
fraction.ones = seq(0.5, 0.5, by = 0.1)
weights.removed = seq(0, .9, by = 0.1)
training.method = "hebb"

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

hebb.learn.patterns = function(patterns) {
  inp.w = matrix(0, no.nodes, no.nodes)
  for (ii in 1:ncol(inp.w))
    for (jj in 1:nrow(inp.w)) {
      inp.w[ii, jj] =
        1.0 / no.nodes * sum(patterns[ii,] * patterns[jj,])
    }
  diag(inp.w) = 0
  return(inp.w)
}

get.h = function(ii, jj, nu, prev.w, patterns) {
  h.sum = 0
  for (ll in (1:ncol(prev.w))[-c(ii, jj)])
    h.sum = h.sum + prev.w[ii, ll] * patterns[ll, nu]
  return(h.sum)
}

storkey.learn.patterns = function(patterns) {
  inp.w = matrix(0, no.nodes, no.nodes)
  prev.w = matrix(0, ncol(inp.w), nrow(inp.w))
  for (nu in 1:ncol(patterns)) {
    for (ii in 1:ncol(inp.w)) {
      for (jj in 1:nrow(inp.w)) {
        inp.w[ii, jj] = prev.w[ii, jj] + 1.0 / no.nodes * (
          patterns[ii, nu] * patterns[jj, nu] - patterns[ii, nu] * get.h(jj, ii, nu, prev.w, patterns) - patterns[jj, nu] *
            get.h(ii, jj, nu, prev.w, patterns)
        )
      }
    }
    prev.w = inp.w
  }
  diag(inp.w) = 0
  return(inp.w)
}

update.async.random = function(input) {
  rand = sample(1:no.nodes, 1)
  input[rand] = ifelse(sum(weights[rand, ] * input) > threshold, 1,-1)
  return(input)
}

update.async.ordered = function(input, index) {
  input[index] = ifelse(sum(weights[index, ] * input) > threshold, 1,-1)
  return(input)
}

output.pattern = function(output, pattern) {
  all(output == pattern)
}

################################################
# Actually run stuff
################################################
patterns.stored = vector("list", length(fraction.ones))
names(patterns.stored) = fraction.ones
for (ii in 1:length(patterns.stored)) {
  patterns.stored[[ii]] = vector("list", length(weights.removed))
  names(patterns.stored[[ii]]) = weights.removed
  for (jj in 1:length(patterns.stored[[ii]])) {
    patterns.stored[[ii]][[jj]] =
      vector("list", length(distortion.fracs))
    names(patterns.stored[[ii]][[jj]]) = distortion.fracs
    for (kk in 1:length(patterns.stored[[ii]][[jj]])) {
      patterns.stored[[ii]][[jj]][[kk]] =
        vector("list", length(no.patterns.stored))
      names(patterns.stored[[ii]][[jj]][[kk]]) = no.patterns.stored
    }
  }
}

# Set training method
if (tolower(training.method) == "hebb") {
  learn.patterns = hebb.learn.patterns
} else if (tolower(training.method) == "storkey") {
  learn.patterns = storkey.learn.patterns
}

data = list() # This we'll treat separately
for (iter in 1:technical.replicates) {
  for (this.no.patterns.stored in no.patterns.stored) {
    for (this.fraction.ones in fraction.ones) {
      # Generate patterns
      patterns = replicate(this.no.patterns.stored,
                           sample(c(
                             rep(1, pattern.length * this.fraction.ones),
                             rep(-1, pattern.length * (1 - this.fraction.ones))
                           ), pattern.length, replace = FALSE))
      
      for (this.weights.removed in weights.removed) {
        # Set the weights
        weights = learn.patterns(patterns)
        to.remove = sample(1:no.nodes, this.weights.removed * no.nodes, replace = FALSE)
        weights[to.remove, to.remove] = 0
        
        for (this.distortion.frac in distortion.fracs) {
          patterns.found = 0 # This is what we measure.
          # Simulate over the patterns
          for (this.pattern in 1:ncol(patterns)) {
            # Do some replicates for data
            for (rep in 1:technical.replicates) {
              output = patterns[, this.pattern]
              
              # Distort the given fraction of bits.
              flipbits = sample(1:pattern.length,
                                this.distortion.frac * pattern.length,
                                replace = FALSE)
              if (length(flipbits) > 0)
                output[flipbits] = -1 * output[flipbits]
              
              # See if convergence.
              for (ii in 1:no.iterations) {
                #print(calc.energy(output))
                output = update.async.random(output) #update.async.ordered(output, ii %% no.nodes + 1)
              }
              # Have we converged? To the right one?
              if (any(apply(patterns, 2, function(x)
                output.pattern(output, x))) &
                output == patterns[, this.pattern]) {
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
