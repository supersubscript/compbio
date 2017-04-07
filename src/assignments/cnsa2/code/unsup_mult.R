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
rho <-  .7

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

train = function(inputs, covorcor, method = "sub", oja.c = 1, ...) {
  weights = runif(no.weights, min = 0, max = 1)
  cor.mat = covorcor(inputs)
  if (method == "mult") {
    for (iter in 1:no.datapoints) {
      weights = weights + c(learning.rate * cor.mat %*% weights) - learning.rate * weights * (weights %*% (cor.mat %*% weights)) # eq. 11.30
    }
  } else if (method == "sub") {
    for (iter in 1:no.datapoints) {
      # weights = weights + learning.rate*(c(cor.mat %*% weights) - 1/no.weights * c(rowSums(cor.mat) %*% weights)) # Preserves weights
      # weights = weights + learning.rate*(c(cor.mat %*% weights) - 1/no.weights * c(weights * c(cor.mat %*% unit(2)))) # Doesn't preserve weights
      
      weights = weights + learning.rate*c(cor.mat %*% weights) - learning.rate/no.weights * sum(cor.mat%*%weights) # Doesn't preserve weights
      # change = learning.rate*(c(cor.mat %*% weights) - 1/no.weights * sum(cor.mat%*%weights))
      # if(any(weights + change > 1)){
      #   newunit = !(weights + change > 1)
      #   print(weights)
      #   print(newunit)
      #   change = learning.rate*(c(cor.mat %*% weights) - 1/sum(newunit) * sum((cor.mat*newunit)%*%weights))#sum(cor.mat%*%weights))
      # }
      # 
      # if(any(weights + change < 0)){
      #   newunit = !(weights + change < 0)
      #   change = learning.rate * ((c(cor.mat%*%newunit) %*% weights) - 1/sum(newunit) * sum((cor.mat*newunit)%*%weights))#sum(cor.mat%*%weights))
        # change = learning.rate*(c(cor.mat %*% weights) - 1/sum(newunit) * asum((cor.mat*newunit)%*%weights))#sum(cor.mat%*%weights))
        # print(change)
      weights[weights > 1] = 1
      weights[weights < 0] = 0
      # print(sum(weights))
      }
      # weights = weights + change
      
      # weights = weights + learning.rate*(c(cor.mat %*% weights) - 1/no.weights * sum(cor.mat%*%weights)) # Preserves weights
      # print(sum(weights))
      # weights[weights > 1] = 1
      # weights[weights < 0] = 0
  }else if (method == "hebb") {
    for (iter in 1:no.datapoints) {
      output = c(weights %*% inputs[iter,])
      weights = weights + learning.rate*(output*inputs[iter,] - output*sum(inputs[iter,] %*% unit(2))*unit(2)/no.weights) # Preserves weights
      # weights = weights / sqrt(weights%*%weights)
    }
  }else if (method == "oja"){
      for (iter in 1:(no.datapoints)) {
        output = c(weights %*% inputs[iter,])
        weights = weights + learning.rate*(output*inputs[iter,] - oja.c*output**2*weights) # Preserves weights
        # weights = weights / sqrt(weights%*%weights)
      }
    }
  return(weights)
}


############################
### Run simulation
############################
plot(
  x,
  y,
  xlim = c(-3, 3) + m,
  ylim = c(-4, 4) + m * rho,
  col = 3, cex = .5
)

# Covariance
weights = train(data, cov, method = "sub")
cat("sub+cov unnorm:\t", weights, "\n")
weights = weights / sqrt(weights %*% weights)
cat("sub+cov norm:\t", weights, "\n")
lines(
   -3:3 * weights[1] + m,
  -3:3 * weights[2] + m * rho,
  col = 1,
  lwd = 5,
  lty = 3
)
# weights = train(data, cov, method = "mult")
# cat("mul+cov unnorm:\t", weights, "\n")
# weights = weights / sqrt(weights %*% weights)
# cat("mul+cov norm:\t", weights, "\n")
# lines(
#   -3:3 * weights[1] + m,
#   -3:3 * weights[2] + m * rho,
#   col = 2,
#   lwd = 5,
#   lty = 3
# )

# Correlation
weights = train(data, cor, method = "sub")
cat("sub+cor unnorm:\t", weights, "\n")
weights = weights / sqrt(weights %*% weights)
cat("sub+cor norm:\t", weights, "\n")
lines(
  -3:3 * weights[1] + m,
  -3:3 * weights[2] + m * rho,
  col = 4,
  lwd = 5,
  lty = 3
)

# weights = train(data, cor, method = "mult")
# cat("mul+cor unnorm:\t", weights, "\n")
# weights = weights / sqrt(weights %*% weights)
# cat("mul+cor norm:\t", weights, "\n")
# print(prcomp(data, scale = T, center = T))
# lines(
#   -3:3 * weights[1] + m,
#   -3:3 * weights[2] + m * rho,
#   col = 5,
#   lwd = 5,
#   lty = 2
# )
print(prcomp(data.frame(data)))