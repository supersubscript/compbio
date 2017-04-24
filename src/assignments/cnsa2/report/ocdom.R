#!/usr/bin/env Rscript
setwd("~/compbio/src/assignments/cnsa2/code/")
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))
par(mfrow = c(1, 1), oma = c(6, 5, 4, 2) + 0.0, mai = c(.0, .6, .0, .5))

# Parameters
no.datapoints = 2000
no.iter.ocdom = 1000
no.cells = 512
rho =  .5
m = 0
sigma = 0.066
learning.rate = 0.01

# Circumference and inter-neuronal distances
a = 10/(no.cells)
total.circ = 10 #+ a

# Preallocate
wminus = matrix(NA, no.cells, no.iter.ocdom)
K = matrix(0, no.cells, no.cells)
K = sapply(1:no.cells, function(ii)
  sapply(1:no.cells, function(jj)
    exp(
      -ifelse(
        abs(ii - jj) < no.cells / 2,
        a * abs(ii - jj),
        total.circ - a * abs(ii - jj)
      ) ** 2 / (2 * sigma ** 2)
    ) - 1 / 9 * exp(
      -ifelse(
        abs(ii - jj) < no.cells / 2,
        a * abs(ii - jj),
        total.circ - a * abs(ii - jj)
      ) ** 2 / (18 * sigma ** 2)
    )))

# Train the weights!
train = function(weights = matrix(runif(no.cells * 2, min = 0, max = 0.01), no.cells, 2),
                 data = cbind(rnorm(no.datapoints, mean = m, sd = 1),
                              rnorm(no.datapoints, mean = rho * x, sd = sqrt(1 - rho ** 2))),
                 ...) {
  Q = cor(data) 
  
  # Update weights
  for (iter in 1:no.iter.ocdom) {
    change = learning.rate * K%*%weights%*%Q # Just calculate this once
    weights = weights +  change - rowMeans(change)
    weights[(weights > 1)] = 1
    weights[(weights < 0)] = 0
    wminus[, iter] = weights[, 2] - weights[, 1]
  }
  list(weights=weights, wminus=wminus)
}

no.replicates = 100
# results = replicate(no.replicates, train())
# save(results, file="../data/ocdom_results.RData")
# load("../data/ocdom_results.RData")

# Perform fourier transform on the 100 iterations.
last.wminus = sapply((results["wminus",]), function(x) x[, no.iter.ocdom]) # Get the final wminus
ffts = apply(last.wminus, 2, fft)
ffts = apply(ffts, 1:2, function(x) Re(sqrt(x**2)))
means = apply(ffts, 1, mean)
stds = apply(ffts, 1, function(x) sd(x) / sqrt(no.replicates))

lw.s = 2
par(mfrow=c(2,1))
plot(means, type = "l", col = 1, lwd = lw.s, xaxt = "n")#, ylim = c(-.1,1))
plot(K[,1], col = 2, lwd = lw.s, type = "l", xlab = "Neuron index", ylab = "Magnitude")
lines(K[,128], col = 3, lwd = lw.s, type = "l", xlab = "Neuron index", ylab = "Magnitude")
lines(K[,256], col = 4, lwd = lw.s, type = "l", xlab = "Neuron index", ylab = "Magnitude")
lines(K[,384], col = 5, lwd = lw.s, type = "l", xlab = "Neuron index", ylab = "Magnitude")
legend("bottomright", col=1:5, bg = "white", lwd = lw.s, c(expression("fft(w"["-"]*")"), expression("K"[1]), expression("K"[128]),expression("K"[256]),expression("K"[384])), lty = c(1), inset = c(0.02,0.02))
# suppressWarnings(arrows(1:no.iter.ocdom, means - stds, 1:no.iter.ocdom, means + stds, length = 0.05, angle = 90, code = 3, cex = 2, lwd = lw.s, col = 1))
mtext(side=1, las = 1, line = 3, "Iteration")

############################
### Run simulation and plot
############################
par(mfrow = c(2, 1), oma = c(6, 5, 0, 2) + 0.0, mai = c(.0, .2, .3, .4))
plot(1, xlim=c(1, no.iter.ocdom), ylim =c(-1,1), cex = 0, bty = "n", xaxt = "n")
mtext(side=2, las = 3, line = 3, "Weight")
train.out = train()
weights = train.out[[1]]
wminus = train.out[[2]]
tmp = sapply(1:no.cells, function(cell) lines(1:no.iter.ocdom, wminus[cell, ], cex = 0.1, col = sample(1:10)))
image(t(wminus), ylab = "", yaxt="n", xlab = "Iterations", xaxt = "n", bty= "n")
mtext(side=2, las = 3, line = 3, "Index")
axis(1, at=seq(0, 1,.20), labels = seq(0, no.iter.ocdom, no.iter.ocdom/5))
axis(2, at=seq(0, 1,.25), labels=c(1, seq(512/4, 512, 512/4)))
mtext(side=1, las = 1, line = 3, "Iteration")





### Save .GIF
# setwd("~/gif")
# frames = 512
# for(i in 1:frames){
#   # creating a name for each plot file with leading zeros
#   if (i < 10) {
#     name = paste('000', i, 'plot.png', sep = '')
#   }
#   else if (i < 100 && i >= 10) {
#     name = paste('00', i, 'plot.png', sep = '')
#   }
#   else if (i >= 100) {
#     name = paste('0', i, 'plot.png', sep = '')
#   }
#   
#   #saves the plot as a .png file in the working directory
#   png(name)
#   plot(fft(K[i,]), type = "l")
#   dev.off()
# }
# cat("sub+cov unnorm:\t", weights, "\n")
# weights = weights / sqrt(weights %*% weights)
# cat("sub+cov norm:\t", weights, "\n")
# print(any(weights != 1))
# lines(-3:3 * weights[1] + m,-3:3 * weights[2] + m * rho, col = 1, lwd = 5, lty = 3)
