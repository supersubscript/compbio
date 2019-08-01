#!/usr/bin/env Rscript

# Params
x.init = 0.5
no.generations = 10000
pop.size = 1000
no.replicates = 1000

# Data
x = vector(mode = "numeric", no.generations)
x[1] = x.init

# Run!
means = vector(mode="numeric", no.replicates)
vars = vector(mode="numeric", no.replicates)

for (jj in 1:no.replicates) {
  for (ii in 2:no.generations) {
    dx = runif(1, min = -0.1, max = 0.1)
    df = abs(x[ii - 1] + dx) - abs(x[ii - 1])
    p.fix = (1 - exp(-2 * df)) / (1 - exp(-2 * pop.size * df))
    x[ii] = ifelse(runif(1) > p.fix, x[ii - 1] + dx, x[ii - 1])
  }
  if(jj == 1){
    plot(x, type = "l", ylim = c(-3, 3), bty = "n")
  }else
    lines(x, col = sample(rainbow(10)))

}

