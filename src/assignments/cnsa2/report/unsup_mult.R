#!/usr/bin/env Rscript
setwd("~/compbio/src/assignments/cnsa2/code/")
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))
learning.rate = 0.001
no.weights = 2
no.datapoints = 2000

# Get the Dayan/Abbott specified correlation matrix
correlation.matrix = function(data){
  t(data)%*%data / nrow(data)
}
covariance.matrix = function(data){
  new.vector = data - matrix(rep(colMeans(data),nrow(data)), byrow = T, ncol=2)
  t(new.vector)%*%(new.vector) / nrow(new.vector)
}

all.weights = matrix(NA, no.weights, no.datapoints + 1)
train = function(inputs, covorcor, method = "mult", ...) {
  weights = runif(no.weights, min = ifelse(method == "sub", 0, 1), max = 1)
  all.weights[,1] = weights
    init.weights = weights
    # weights = weights / sqrt(weights %*% weights)
  cor.mat = covorcor(inputs)
  if (method == "mult") {
    for (iter in 1:no.datapoints) {

        if(sum(weights) != 0)
          weights = weights + learning.rate * (c(cor.mat %*% weights) - sum(t(weights) %*% cor.mat) / sum(weights)  * weights)
        all.weights[,iter] = weights
      # print(sum(weights))
    }
  } else if (method == "sub") {
    for (iter in 1:no.datapoints) {
      weights = weights + learning.rate*(c(cor.mat%*%weights) - 1/no.weights * sum(cor.mat %*% weights))
      all.weights[,iter+1] = weights
      if(any(weights > 1 | weights < 0)){
        weights[weights > 1] = 1
        weights[weights < 0] = 0
        all.weights[,iter:no.datapoints] = weights
        break
      }
      all.weights[,iter + 1] = weights
    }
    }
  list(final=weights, init=init.weights, all.weights=all.weights)
}

train2 = function(inputs, covorcor, method = "mult", no.iterations = 10000, weights, ...) {
  init.weights = weights
  all.weights = matrix(NA, no.weights, no.iterations + 1)
  all.weights[,1] = weights
  # weights = weights / sqrt(weights %*% weights)
  cor.mat = covorcor(inputs)
  if (method == "mult") {
    for (iter in 1:no.iterations) {
      # weights = weights + c(learning.rate * cor.mat %*% weights) - learning.rate * weights * c(weights %*% c(cor.mat %*% weights)) # eq. 11.30, works fine
      if(sum(weights) != 0)
        weights = weights + learning.rate * (c(cor.mat %*% weights) - sum(t(weights) %*% cor.mat) / sum(weights)  * weights)
      all.weights[,iter+1] = weights
    }
  } else if (method == "sub") {
    for (iter in 1:no.iterations) {
      weights = weights + learning.rate*(c(cor.mat%*%weights) - 1/no.weights * sum(cor.mat %*% weights))
      if(any(weights > 1 | weights < 0)){
        weights[weights > 1] = 1
        weights[weights < 0] = 0
        all.weights[,iter:no.datapoints] = weights
        break
      }
      all.weights[,iter+1] = weights
    }
  }
  list(final=weights / sqrt(weights[1]**2  + weights[2]**2), init=init.weights, all.weights = all.weights)
}

####################################
### SIMULATE
####################################
global.rho = .5
generate.data = function(rho, m){
  x = rnorm(no.datapoints, mean = m, sd = 1)
  y = rnorm(no.datapoints, mean = rho*x, sd = sqrt(1 - rho ** 2))
  data = cbind(x, y)
  data
}

m = 0
rho = global.rho
dat = generate.data(rho, m)
x = dat[,1]
y = dat[,2]

############################
### Run simulation
############################
par(mfrow = c(2, 2), oma = c(6, 5, 4, 2) + 0.0, mai = c(.0, .6, .2, .5))
plot(x, y, xlim = c(-3, 3) + mean(x), ylim = c(-4, 4) + mean(y), col = alpha(3, 1), cex = .5, bty ="n", xaxt = "n", ylab = "")
weights = train2(dat, correlation.matrix, method = "mult", weights = c(-1,-1))$final
# weights = weights / sqrt(weights[1]**2  + weights[2]**2)
lines(
   -3:3  + m,
   -3:3 * weights[2] / weights[1] + m*rho, #weights[2] + m * rho,
   col = 1,
  lwd = 5,
  lty = 1
)
weights = train(dat, covariance.matrix, method = "mult", weights = c(-1,-1))$final
# weights = weights / sqrt(weights %*% weights)
weights = weights / sqrt(weights[1]**2  + weights[2]**2)
lines(
  -3:3  + m,
  -3:3 * weights[2] / weights[1] + m*rho, #weights[2] + m * rho,
  col = 2,
  lwd = 6,
  lty = 2
)

# Non-zero
m = 10
rho = global.rho
dat = generate.data(rho, m)
x = dat[,1]
y = dat[,2]

plot(x, y, xlim = c(-3, 3) + mean(x), ylim = c(-4, 4) + mean(y), col = 3, cex = .5, bty ="n", xaxt = "n", ylab = "")
weights = train(dat, correlation.matrix, method = "mult")$final
# weights = weights / sqrt(weights %*% weights)
weights = weights / sqrt(weights[1]**2  + weights[2]**2)

lines(
  -3:3  + m,
  -3:3 * weights[2] / weights[1] + mean(y), #weights[2] + m * rho,
  col = 1,
  lwd = 5,
  lty = 2
)

weights = train(dat, covariance.matrix, method = "mult")$final
# weights = weights / sqrt(weights %*% weights)
weights = weights / sqrt(weights[1]**2  + weights[2]**2)
lines(
  -3:3  + m,
  -3:3 * weights[2] / weights[1] + mean(y), #weights[2] + m * rho,
  col = 2,
  lwd = 5,
  lty = 3, ylab = ""
)
# legend("bottomright", c("Cor", "Cov"), lty = c(1,2), lwd = c(6,5), col = 1:2, inset = c(0.02,0.02))


m = 0
rho = global.rho
dat = generate.data(rho, m)
x = dat[,1]
y = dat[,2]

############################
### Run simulation
############################
plot(x, y, xlim = c(-3, 3) + mean(x), ylim = c(-4, 4) + mean(y), col = alpha(3, 1), cex = .5, bty ="n", ylab = "")
weights = train(dat, correlation.matrix, method = "sub")$final
# weights = weights / sqrt(weights[1]**2  + weights[2]**2)
print(weights)

lines(
  -3:3  + m,
  -3:3 * weights[2] / ifelse(weights[1]==0, 0.00000001, weights[1]) + m*rho, #weights[2] + m * rho,
  col = 1,
  lwd = 5,
  lty = 2
)
weights = train(dat, covariance.matrix, method = "sub")$final
# weights = weights / sqrt(weights %*% weights)
# weights = weights / sqrt(weights[1]**2  + weights[2]**2)
lines(
  -3:3  + m,
  -3:3 * weights[2] / ifelse(weights[1]==0, 0.00000001, weights[1]) + m*rho, #weights[2] + m * rho,
  col = 2,
  lwd = 6,
  lty = 3
)
mtext(side = 2, line = 3, at = 5, las = 3, "y")
mtext(side = 1, line = 3, at = 5, las = 1, "x")
mtext(side = 2, line = 4, at = 10, las = 3, "multiplicative")
mtext(side = 2, line = 4, at = 0, las = 3, "subtractive")
# print(weights)



# Non-zeor
m = 10
rho = global.rho
dat = generate.data(rho, m)
x = dat[,1]
y = dat[,2]

plot(x, y, xlim = c(-3, 3) + mean(x), ylim = c(-4, 4) + mean(y), col = 3, cex = .5, bty ="n", ylab = "")
weights = train(dat, correlation.matrix, method = "sub")$final
# weights = weights / sqrt(weights %*% weights)
# weights = weights / sqrt(weights[1]**2  + weights[2]**2)
print(weights)

lines(
  -3:3  + m,
  -3:3 * weights[2] / ifelse(weights[1]==0, 0.00000001, weights[1]) + mean(y), #weights[2] + m * rho,
  col = 1,
  lwd = 5,
  lty = 2
)

weights = train(dat, covariance.matrix, method = "sub")$final
# weights = weights / sqrt(weights %*% weights)
# weights = weights / sqrt(weights[1]**2  + weights[2]**2)
lines(
  -3:3  + m,
  -3:3 * weights[2] / ifelse(weights[1]==0, 0.00000001, weights[1]) + mean(y), #weights[2] + m * rho,
  col = 2,
  lwd = 5,
  lty = 3
)
print(weights)
legend("bottomright", c("Cor", "Cov"), lty = c(1,2), lwd = c(6,5), col = 1:2, inset = c(0.02,0.02))

print(prcomp(data.frame(dat)))



#################
par(mfrow=c(2,2))
w.seq = seq(-1,1, 0.2)
dat = generate.data(.5, 10)
w.mat = sapply(w.seq, function(x) sapply(w.seq, function(y) prod(train2(dat, correlation.matrix, method = "mult", weights =c(x,y))$final)))
image(w.mat, xaxt = "n", yaxt = "n", col = heat.colors(12))
mtext(side=3, line = 1, "correlation", las = 1)
mtext(side=2, line = 3, "multiplicative", las = 3)
axis(2, at=seq(0,1,.25), labels = seq(-1,1, 0.5))
w.mat = sapply(w.seq, function(x) sapply(w.seq, function(y) prod(train2(dat, covariance.matrix, method = "mult", weights =c(x,y))$final)))
image(w.mat, xaxt = "n", yaxt = "n", col = heat.colors(12))
mtext(side=3, line = 1, "covariance", las = 1)

w.seq = seq(0,1, 0.1)
w.mat = sapply(w.seq, function(x) sapply(w.seq, function(y) train2(dat, correlation.matrix, method = "sub", weights =c(x,y))$final[1]))
image(w.mat, xaxt = "n", yaxt = "n", col = heat.colors(12))
axis(2, at=seq(0,1,.25), labels = seq(0,1, 0.25))
axis(1, at=seq(0,1,.25), labels = seq(-1,1, 0.5))
mtext(side=2, line = 3, "subtractive", las = 3)
w.mat = sapply(w.seq, function(x) sapply(w.seq, function(y) train2(dat, covariance.matrix, method = "sub", weights =c(x,y))$final[1]))
image(w.mat, xaxt = "n", yaxt = "n", col = heat.colors(12))
axis(1, at=seq(0,1,.25), labels = seq(-1,1, 0.5))
# 

lw.s = 3
par(mfrow = c(2, 2), oma = c(6, 5, 2, 4) + 0.0, mai = c(.2, .2, .2, .2))
xmax = 30
w.weights = train2(dat, correlation.matrix, method = "mult", weights =c(.9, .1))$all.weights
plot(w.weights[1,], type = "l", ylim = c(-1, 1), col = 1, xlim = c(0, xmax), lwd = lw.s, xaxt = "n")
lines(w.weights[2,], col = 1, lwd = lw.s)
mtext(side = 2, at = .0, "multiplicative", las = 3, line = 4)
w.weights = train2(dat, correlation.matrix, method = "mult", weights =c(.9, .1))$all.weights
lines(w.weights[1,], type = "l", ylim = c(-1, 1), col = 2, xlim = c(0, xmax), lwd = lw.s, xaxt = "n", yaxt = "n", ylab = "")
lines(w.weights[2,], col = 2, lwd = lw.s)

w.weights = train2(dat, correlation.matrix, method = "mult", weights =c(-.6, .5))$all.weights
plot(w.weights[1,], type = "l", ylim = c(-1, 1), col = 1, xlim = c(0, xmax), lwd = lw.s, xaxt = "n", yaxt = "n", ylab = "")
lines(w.weights[2,], col = 1, lwd = lw.s)
w.weights = train2(dat, correlation.matrix, method = "mult", weights =c(-.5, .6))$all.weights
lines(w.weights[1,], type = "l", ylim = c(-1, 1), col = 2, xlim = c(0, xmax), lwd = lw.s, xaxt = "n", yaxt = "n", ylab = "")
lines(w.weights[2,], col = 2, lwd = lw.s)

w.weights = train2(dat, correlation.matrix, method = "sub", weights =c(.1,.9))$all.weights
plot(w.weights[1,], type = "l", ylim = c(0, 1), col = 1, xlim = c(0, xmax), lwd = lw.s)
lines(w.weights[2,], col = 1, lwd = lw.s)
mtext(side = 2, at = 1, "Weight", las = 3, line = 3)
mtext(side = 2, at = .5, "subtractive", las = 3, line = 4)

w.weights = train2(dat, correlation.matrix, method = "sub", weights =c(.4,.3))$all.weights
plot(w.weights[1,], type = "l", ylim = c(0, 1), col = 1, xlim = c(0, xmax), lwd = lw.s, yaxt = "n", ylab = "")
lines(w.weights[2,], col = 1, lwd = lw.s)
mtext(side = 1, at = -1.5, "Iteration", las = 1, line = 3)




# image(w.mat, xaxt = "n", yaxt = "n", col = 1:8)
# mtext(side=3, line = 1, "correlation", las = 1)
# mtext(side=2, line = 3, "multiplicative", las = 3)
# axis(2, at=seq(0,1,.25), labels = seq(-1,1, 0.5))
# w.mat = sapply(w.seq, function(x) sapply(w.seq, function(y) sum(train2(dat, covariance.matrix, method = "mult", weights =c(x,y))$final)))
# image(w.mat, xaxt = "n", yaxt = "n", col = 1:8)
# mtext(side=3, line = 1, "covariance", las = 1)

# image(w.mat, xaxt = "n", yaxt = "n", col = 1:8)
# axis(2, at=seq(0,1,.25), labels = seq(0,1, 0.25))
# axis(1, at=seq(0,1,.25), labels = seq(-1,1, 0.5))
# mtext(side=2, line = 3, "subtractive", las = 3)
# w.mat = sapply(w.seq, function(x) sapply(w.seq, function(y) train2(dat, covariance.matrix, method = "sub", weights =c(x,y))$final[1]))
# image(w.mat, xaxt = "n", yaxt = "n", col = 1:8)
# axis(1, at=seq(0,1,.25), labels = seq(-1,1, 0.5))

