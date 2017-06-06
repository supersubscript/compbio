#!/usr/bin/env Rscript
setwd("/local/data/public/hpa22/assignments/syba1")
rm(list = ls())

# Import arguments from BASH script (choice of tau1)
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
  stop("This script needs 1 argument: tau0 index.")
}

tau0.idx <- args[1]

# Simulation parameters
no.iterations = 20e6
no.reactions  = 6
no.species    = 3
epsilon       = 10e-3
reaction.sign = c(1, -1, 1, -1, 1, -1)
reactions = rep(NA, no.reactions)
species.index = c(1, 1, 2, 2, 3, 3)

# Iteration independent parameters
gamma    = 29 # Birthday
avg.x0   = 10 
avg.x1   = 10
avg.x2   = 1000
start.x0 = round(avg.x0) 
start.x1 = round(avg.x1) 
start.x2 = round(avg.x2)
start.t  = 0
tau1     = round((2 + 8 * gamma / 31) * 60) 
tau2     = 180 * 60
tau0     = seq(tau2 / 1000, tau1, length.out = 100)
beta1    = 1 / tau1
beta2    = 1 / tau2
lambda1  = avg.x1 * beta1
lambda2  = avg.x2 * beta2 / (avg.x0 * avg.x1)

# Calculate a weighted (co)variance
weighted.cov = function(var1, var2 = var1, weights, w.mean1, w.mean2 = w.mean1){
  1 / sum(weights) * sum(weights * (var1 - w.mean1) * (var2 - w.mean2))
}

step.size = 100000
start.pt = 100000
extras = matrix(NA, ncol = 15, nrow = (no.iterations-start.pt) / step.size)
count = 1


# Main function
# tau0 = tau2 / 1000 # Test with this
gillespie = function(tau0 = tau0, no.iterations = no.iterations){
  # System parameters
  # This time we change only for tau0
  beta0 = 1 / tau0
  lambda0 = avg.x0 * beta0
  
  # Preallocation and definition of output matrix
  probabilities = rep(NA, no.reactions) 
  out.values = matrix(NA, ncol = 1 + no.species, nrow = no.iterations)
  out.values[1, ] = c(start.t, start.x0, start.x1, start.x2)
  ii = 2
  colnames(out.values) = c("t", "x0","x1","x2")
  
  while(ii < no.iterations + 1) { 
    # Define reaction probabilities
    reactions = c(
      lambda0,
      beta0 * out.values[ii - 1, 2],
      lambda1,
      beta1   * out.values[ii - 1, 3],
      lambda2 * out.values[ii - 1, 2] * out.values[ii - 1, 3],
      beta2   * out.values[ii - 1, 4])
    sum.reactions = sum(reactions) # This is our "a_0" constant.
    
    # Increment / decrement reacting species, choosing one of the 6 possible reactions
    reaction.index = sample(1:no.reactions, 1, prob = reactions / sum.reactions) 
    reaction.species = 1:no.species == ceiling(reaction.index / 2) 
    
    # Update the time and chosen spec ies
    out.values[ii, 1:4] = out.values[ii - 1, 1:4] +
      c(rexp(1, sum.reactions), # Increment time (rexp is faster than runif trickery)
      reaction.sign[reaction.index] * reaction.species) # Update the correct species
    
    # Compute running average of each species, every 100000 iterations
    if ((ii - 1) %% step.size == 0 && ii > start.pt) {
      cat("iteration ", ii, " of ", no.iterations, "\n")
      wts = diff(out.values[start.pt:ii,  "t"])
      x0s = out.values[start.pt:(ii - 1), "x0"]
      x1s = out.values[start.pt:(ii - 1), "x1"]
      x2s = out.values[start.pt:(ii - 1), "x2"]

      # Weighted mean depending on different reaction times
      x0.mean = weighted.mean(x0s, wts)
      x1.mean = weighted.mean(x1s, wts)
      x2.mean = weighted.mean(x2s, wts)
      
      # Flux balance
      birth.x0 = lambda0
      death.x0 = beta0   * x0.mean
      birth.x1 = lambda1
      death.x1 = beta1   * x1.mean
      birth.x2 = lambda2 * x0.mean * x1.mean
      death.x2 = beta2   * x2.mean
      
      # Return the relative errors in balance
      extras[count, 1:3] = c((birth.x0 - death.x0) / birth.x0,
                             (birth.x1 - death.x1) / birth.x1,
                             (birth.x2 - death.x2) / birth.x2)
      
      # Fluctuation balance (from simulation and expected)
      n.00 = weighted.cov(var1 = x0s, var2 = x0s, weights = wts, w.mean1 = x0.mean, w.mean2 = x0.mean) / (x0.mean ** 2)
      n.01 = weighted.cov(var1 = x0s, var2 = x1s, weights = wts, w.mean1 = x0.mean, w.mean2 = x1.mean) / (x0.mean * x1.mean)
      n.02 = weighted.cov(var1 = x0s, var2 = x2s, weights = wts, w.mean1 = x0.mean, w.mean2 = x2.mean) / (x0.mean * x2.mean)
      n.11 = weighted.cov(var1 = x1s, var2 = x1s, weights = wts, w.mean1 = x1.mean, w.mean2 = x1.mean) / (x1.mean ** 2)
      n.12 = weighted.cov(var1 = x1s, var2 = x2s, weights = wts, w.mean1 = x1.mean, w.mean2 = x2.mean) / (x1.mean * x2.mean)
      n.22 = weighted.cov(var1 = x2s, var2 = x2s, weights = wts, w.mean1 = x2.mean, w.mean2 = x2.mean) / (x2.mean ** 2)
      
      # Theoretical answers
      n.00.exp = 1 / x0.mean
      n.01.exp = 0
      n.02.exp = tau0 / (tau0 + tau2) * n.00.exp
      n.11.exp = 1 / x1.mean
      n.12.exp = tau1 / (tau1 + tau2) * n.11.exp
      n.22.exp = 1 / x2.mean + n.02.exp + n.12.exp
      
      # Return the relative percent difference in balance
      extras[count, 4:9] = c((n.00 - n.00.exp) / n.00.exp,
                             (n.01 - n.01.exp),
                             (n.02 - n.02.exp) / n.02.exp,
                             (n.11 - n.11.exp) / n.11.exp,
                             (n.12 - n.12.exp) / n.12.exp,
                             (n.22 - n.22.exp) / n.22.exp)
      
      # Return normalised covariances
      extras[count, 10:15] = c(n.00, n.01, n.02, n.11, n.12, n.22)
    
      # Check if the errors are smaller than a threshold, and stop the simulation if true
      if(all(abs(extras[count, 1:9]) < epsilon)) {
        return(list(out.values, extras))
      }
      count = count + 1
    }
    ii = ii + 1
  }
  return(list(out.values, extras))
}

# Run the simulation
# a <- gillespie(tau1[as.numeric(tau1.idx)], no.iterations=no.iterations)
results = gillespie(tau0[as.numeric(tau0.idx)], no.iterations=no.iterations)
# results = gillespie(tau0, no.iterations=no.iterations)

# Store the resulting matrix in an .RData object
save(results, file = paste0(tau0.idx, "_result.RData"))

# All RData files save a matrix named "results", so loading any file 
# will import an "a" matrix (watch out for overwriting!)

