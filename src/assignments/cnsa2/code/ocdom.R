#!/usr/bin/env Rscript
setwd("~/compbio/src/assignments/cnsa2/code/")
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))
par(
  mfrow = c(1, 1),
  oma = c(6, 5, 4, 2) + 0.0,
  mai = c(.0, .6, .0, .5)
)

no.datapoints <- 2000
no.iter.ocdom = 1000
no.cells = 512
rho <-  0
m = 2
x <- rnorm(no.datapoints, mean = m, sd = 1)
y <- rnorm(no.datapoints, mean = rho * x, sd = sqrt(1 - rho ** 2))
data = cbind(x, y)
learning.rate = 0.05
no.weights = 2
# parint(var(x,y))
unit = function(len) {
  rep(1, len)
}

wminus = matrix(NA, no.cells, no.iter.ocdom)
train = function(inputs, covorcor, method = "sub", oja.c = 1, ...) {
    weights = matrix(runif(no.cells*2, min = 0, max = 0.1),no.cells,2)#runif(no.weights, min = 0, max = 1)
    K = matrix(0, no.cells, no.cells)
    sigma = 0.066
    cor.mat = cor(inputs)
    learning.rate = 0.01
    for(ii in 1:no.cells){
      for(jj in 1:no.cells){
        K[ii,jj] = exp(-(ii-jj)**2/(2*sigma**2)) - 1/9 *exp(-(ii-jj)**2/(18*sigma**2))
      }
    }
    for (iter in 1:no.iter.ocdom) {
      weights = weights + learning.rate*K%*%weights%*%cor.mat - learning.rate * sum(K%*%weights%*%cor.mat) / (2*no.cells)
      weights[which(weights > 1)] = 1
      weights[which(weights < 0)] = 0
      wminus[,iter] <<- weights[,2] - weights[,1]
    }
  return(weights)
}


############################
### Run simulation
############################
# plot(x, y, xlim = c(-3, 3) + m, ylim = c(-4, 4) + m * rho, col = 3, cex = .5)
# wminus[which(wminus > 1)] = 1
# wminus[which(wminus < 0)] = 0
plot(1, xlim=c(1, no.iter.ocdom), ylim =c(-1,1), cex = 0)
# print(any(weights != 1))
weights = train(data, cov, method = "ocdom")
tmp = sapply(1:no.cells, function(cell) lines(1:no.iter.ocdom, wminus[cell, ], cex = 0.1))
cat("sub+cov unnorm:\t", weights, "\n")
# weights = weights / sqrt(weights %*% weights)
cat("sub+cov norm:\t", weights, "\n")
print(any(weights != 1))
# lines(-3:3 * weights[1] + m,-3:3 * weights[2] + m * rho, col = 1, lwd = 5, lty = 3)
