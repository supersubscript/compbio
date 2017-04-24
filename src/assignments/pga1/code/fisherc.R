#!/usr/bin/env Rscript
setwd("~/compbio/src/assignments/pga1/code/")
library(RColorBrewer)
library(scales) #imports alpha
library(stats)
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 3

# Params
x.init = 0.5
no.generations = 5000
pop.size = 0
no.replicates = 1
lambda = 0.001
linew = 3

# Data
x = vector(mode = "numeric", no.generations)
sum.df = vector(mode = "numeric", no.generations)
x[1] = x.init

# Run!
run.sim = function(poppy.size) {
  target = 0
  for (t in 2:no.generations) {
    target = sin(2 * pi * lambda * t) / 2 + 1 / 2
    dx = runif(1, min = -0.1, max = 0.1)
    if (x[t - 1] + dx < 0) {
      if (x[t - 1] == 0) {
        x[t] = x[t - 1]
        sum.df[t] = sum.df[t - 1]
        next
      }
      dx = -x[t - 1]
    }
    
    df = -(abs(x[t - 1] + dx - target) - abs(x[t - 1] - target))
    p.fix = (1 - exp(-2 * df)) / (1 - exp(-2 * poppy.size * df))
    
    if (runif(1) < p.fix) {
      x[t] = x[t - 1] + dx
      sum.df[t] = sum.df[t - 1] + df
    } else {
      x[t] = x[t - 1]
      sum.df[t] = sum.df[t - 1]
    }
  }
  return(list(x, sum.df))
}

data10 = replicate(no.replicates, run.sim(10))
data10k = replicate(no.replicates, run.sim(1000))
xs.10  = matrix(unlist(data10[1,]), ncol = no.replicates, byrow = FALSE)
fs.10  = matrix(unlist(data10[2,]), ncol = no.replicates, byrow = FALSE)
xs.10k = matrix(unlist(data10k[1,]), ncol = no.replicates, byrow = FALSE)
fs.10k = matrix(unlist(data10k[2,]), ncol = no.replicates, byrow = FALSE)

mean10 = apply(xs.10, 1, mean)
mean10k = apply(xs.10k, 1, mean)
stdev10 = apply(xs.10, 1, function(x)
  sd(x) / sqrt(length(mean10)))
stdev10k = apply(xs.10k, 1, function(x)
  sd(x) / sqrt(length(mean10k)))

f.mean10 = apply(fs.10, 1, mean)
f.mean10k = apply(fs.10k, 1, mean)
f.stdev10 = apply(fs.10, 1, function(x)
  sd(x) / sqrt(length(mean10)))
f.stdev10k = apply(fs.10k, 1, function(x)
  sd(x) / sqrt(length(mean10k)))

par(
  mfrow = c(2, 1),
  oma = c(5, 2, 2, 0) + 0.0,
  mai = c(.0, 1, .0, 1)
)
plot(
  mean10,
  type = "l",
  #  ylim = c(-0, 1),
  ylim = c(-0.05, 1.1),
  col = 1,
  xaxt = "n",
  xlab = "",
  ylab = "x",
  lwd = linew
)
lines(mean10k,
      type = "l",
      col = 2,
      lwd = linew)
lines(
  1:no.generations,
  sin(2 * pi * lambda * 1:no.generations) / 2 + 1 / 2,
  lty = 2,
  lwd = linew
)
skip = 250
arrows((2:no.generations)[seq(2, no.generations, skip)],
       mean10[seq(1, no.generations, skip)] - stdev10[seq(2, no.generations, skip)],
       (2:no.generations)[seq(2, no.generations, skip)],
       mean10[seq(2, no.generations, skip)] + stdev10[seq(2, no.generations, skip)],
       length = 0.05,
       angle = 90,
       code = 3,
       col = 1,
       lwd = linew
)
arrows((2:no.generations)[seq(2, no.generations, skip)],
       mean10k[seq(1, no.generations, skip)] - stdev10k[seq(2, no.generations, skip)],
       (2:no.generations)[seq(2, no.generations, skip)],
       mean10k[seq(2, no.generations, skip)] + stdev10k[seq(2, no.generations, skip)],
       length = 0.05,
       angle = 90,
       code = 3,
       col = 2,
       lwd = linew
)

# Plot those f-stats
plot(
  f.mean10,
  col = 1,
  type = "l",
  ylab = "Summed change \n in fitness",
  xlab = "Generation",
  ylim=c(0,11),
  lwd = linew
)
mtext("Generation", side = 1, line = 2.5)
lines(f.mean10k, col = 2, lwd = 2)
arrows((2:no.generations)[seq(2, no.generations, skip)],
       f.mean10[seq(1, no.generations, skip)] - f.stdev10[seq(2, no.generations, skip)],
       (2:no.generations)[seq(2, no.generations, skip)],
       f.mean10[seq(2, no.generations, skip)] + f.stdev10[seq(2, no.generations, skip)],
       length = 0.05,
       angle = 90,
       code = 3,
       col = 1,
       lwd = linew
)
arrows((2:no.generations)[seq(2, no.generations, skip)],
       f.mean10k[seq(1, no.generations, skip)] - f.stdev10k[seq(2, no.generations, skip)],
       (2:no.generations)[seq(2, no.generations, skip)],
       f.mean10k[seq(2, no.generations, skip)] + f.stdev10k[seq(2, no.generations, skip)],
       length = 0.05,
       angle = 90,
       code = 3,
       col = 2,
       lwd = linew
)

legend(
  "topleft",
  c("N = 10", "N = 1000"),
  inset = 0.02,
  cex = 1.5,
  col = c(1, 2),
  lty = c(1, 1),
  bg = "white",
  lwd = linew
)


# plot(apply(xs.10, 1, sd), type = "l", ylim=c(-2.5,2.5))
# lines(apply(xs.10k, 1, sd), type = "l", col = "red")
