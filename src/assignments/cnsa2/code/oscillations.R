#!/usr/bin/env Rscript
setwd("~/compbio/src/assignments/cnsa2/code/")
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))

# Dynamical parameters
method = "rk4"
no.replicates = 1
epsilon = 1
Ie = -1
# start.v = -0.6823278 + runif(1, max= 0.0001, min = -0.0001)
# start.u = -2*0.6823278 + runif(1, max= 0.0001, min = -0.0001)#runif(1, max = 1, min = -1)
start.v = runif(1, max = 1, min = -1)
start.u = runif(1, max = 1, min = -1)

print.ps = 4
# vol = c(3, 999999999, 9999999999)

# Simulation parameters
start.t = 0.0
end.t = 500.0
h = 0.1
no.iterations = (end.t - start.t) / h
output = data.frame(v = rep(NA, (end.t - start.t)*print.ps), u = rep(NA, (end.t - start.t)*print.ps))

# Init
initialise = function(v, u, epsilon, Ie){
  v <<- start.v
  u <<- start.u
  epsilon <<- epsilon
  Ie <<- Ie
  t <<- start.t
  output[1, ] <<- c(v, u)
  last.t <<- start.t
}

# Define derivatives
v.derivs = function(t, v, u) {
  v * (1 - v ** 2) - u + Ie
}

u.derivs = function(t, v, u) {
  epsilon * (v - 0.5 * u)
}

# Simple Euler integration
euler = function(t, v, u, h, ...) {
  v.k1 = v.derivs(t, v, u)
  u.k1 = u.derivs(t, v, u)
  v <<- v + v.k1 * h
  u <<- u + u.k1 * h
  t <<- t + h
}

# Fourth-order Runge-Kutta
rk4 = function(t, v, u, h, vol) {
  v.k1 = v.derivs(t, v, u)
  u.k1 = u.derivs(t, v, u)
  v.k2 = v.derivs(t + h / 2.0, v + h / 2.0 * v.k1, u + h / 2.0 * u.k1)
  u.k2 = u.derivs(t + h / 2.0, v + h / 2.0 * v.k1, u + h / 2.0 * u.k1)
  v.k3 = v.derivs(t + h / 2.0, v + h / 2.0 * v.k2, u + h / 2.0 * u.k2)
  u.k3 = u.derivs(t + h / 2.0, v + h / 2.0 * v.k2, u + h / 2.0 * u.k2)
  v.k4 = v.derivs(t + h, v + h * v.k3, u + h * u.k3)
  u.k4 = u.derivs(t + h, v + h * v.k3, u + h * u.k3)
  v <<- v + h / 6.0 * (v.k1 + 2 * v.k2 + 2 * v.k3 + v.k4)
  u <<- u + h / 6.0 * (u.k1 + 2 * u.k2 + 2 * u.k3 + u.k4)
  t <<- t + h
}

# Set numerical method
if (tolower(method) == "milstein") {
  fct = milstein
} else if (tolower(method) == "euler") {
  fct = euler
} else if (tolower(method) == "rk4") {
  fct = rk4
} else {
  stop("Error: Method not found.")
}

runSimulation = function(v.init, u.init, epsilon, Ie) {
  initialise(v.init, u.init, epsilon, Ie)
  for (iter in 1:(no.iterations-1)) {
  fct(t, v, u, h, vol)
  if (iter * h > last.t + 1/print.ps) {
    last.t = last.t + 1/print.ps
    output[last.t * print.ps+1, ] = c(v, u)
  }
  }
  return(output)
}
# Run single-valued simulations
# results = replicate(no.replicates, runSimulation(start.v, start.u, epsilon, Ie))

get.period = function(data){
  # Get minima
  indices = which(data[2:(length(data)-1)] < data[3:length(data)] & data[2:(length(data)-1)] < data[1:(length(data)-2)])
  times = indices + 1 # +1 due to offset
  periods = times[2:length(times)] - times[1:(length(times)-1)]
  return(periods)
}


# # Run simulations for multiple parameters
# start.vs = seq(0,1,.25)
# start.vs = c(start.vs, seq(0,1,.25))
# start.vs = c(start.vs, seq(0,1,.25))
# start.us = seq(0,1,.25)
# start.us = c(start.us, seq(0,1,.25))
# start.us = c(start.us, seq(0,1,.25))
# epsilons = c(rep(0.1, length(start.vs) / 3), rep(0.3, length(start.vs) / 3), rep(1, length(start.vs) / 3))
# results = sapply(1:no.replicates, function(x) runSimulation(start.vs[x], start.us[x], epsilons[x], Ies[x]))

results = replicate(no.replicates, runSimulation(start.v, start.u, .1, 0))
results.2 = replicate(no.replicates, runSimulation(start.v, start.u, .3, 0))
results.3 = replicate(no.replicates, runSimulation(start.v, start.u, 1, 0))
par(
  mfrow = c(3, 1),
  oma = c(10, 5, 4, 2) + 0.0,
  mai = c(.0, .6, .0, .5)
)
linewidth = 2.5
# plot(seq(start.t,end.t - 1/print.ps, 1 / print.ps), results[,1]$u, type = "l", xaxt = "n", xlab = "n", ylab = "u", lwd = linewidth, col = 1) # Plot u
# plot(seq(start.t,end.t - 1/print.ps, 1 / print.ps), results[,1]$v, type = "l", xaxt = "n", xlab = "n", ylab = "v", lwd = linewidth, col = 2) # Plot v
# plot(results[,1]$u, results[,1]$v, type = "l", ylab = "v", xlab = "u", lwd = linewidth, col = 1, xaxt = "n") # Plot u vs v
plot(results[,1]$u, results[,1]$v, type = "l", ylab = "v", xlab = "u", lwd = linewidth, col = 1, xlim = c(-1.1,1.1), ylim = c(-1.1,1.1), xaxt = "n") # Plot u vs v
legend("topleft", c(expression(epsilon == 0.1), expression(epsilon == 0.3), expression(epsilon == 1.0)), col = c(1,2,3), inset = c(0.01, 0.05), lty = 1, lwd = linewidth, bg="white", cex = 1.3, box.lwd = 0)
plot(results.2[,1]$u, results.2[,1]$v, type = "l", ylab = "v", xlab = "u", lwd = linewidth, col = 2, xaxt = "n") # Plot u vs v
plot(results.3[,1]$u, results.3[,1]$v, type = "l", ylab = "v", xlab = "u", lwd = linewidth, col = 3) # Plot u vs v
mtext("u", side=1, line=3,las=1, col="black", cex = .75)

# means = apply(results, 1, function(x) mean(unlist(lapply(x, function(y) y[end.t])))) # Get means over many simulations
# stds = apply(results, 1, function(x) sd(unlist(lapply(x, function(y) y[end.t])))/sqrt(no.replicates)) # Ibid (stds)

######## Ie = -1
results = replicate(no.replicates, runSimulation(start.v, start.u, .1, -1))
results.2 = replicate(no.replicates, runSimulation(start.v, start.u, .3, -1))
results.3 = replicate(no.replicates, runSimulation(start.v, start.u, 1, -1))
par(
  mfrow = c(3, 1),
  oma = c(10, 5, 4, 2) + 0.0,
  mai = c(.0, .6, .0, .5)
)
linewidth = 2.5
# plot(seq(start.t,end.t - 1/print.ps, 1 / print.ps), results[,1]$u, type = "l", xaxt = "n", xlab = "n", ylab = "u", lwd = linewidth, col = 1) # Plot u
# plot(seq(start.t,end.t - 1/print.ps, 1 / print.ps), results[,1]$v, type = "l", xaxt = "n", xlab = "n", ylab = "v", lwd = linewidth, col = 2) # Plot v
# plot(results[,1]$u, results[,1]$v, type = "l", ylab = "v", xlab = "u", lwd = linewidth, col = 1, xaxt = "n") # Plot u vs v
plot(results[,1]$u, results[,1]$v, type = "l", ylab = "v", xlab = "u", lwd = linewidth, col = 1, xlim = c(-2.1,0.5), ylim = c(-1.5,1.1), xaxt = "n") # Plot u vs v
legend("topleft", c(expression(epsilon == 0.1), expression(epsilon == 0.3), expression(epsilon == 1.0)), col = c(1,2,3), inset = c(0.01, 0.05), lty = 1, lwd = linewidth, bg="white", cex = 1.3, box.lwd = 0)
plot(results.2[,1]$u, results.2[,1]$v, type = "l", ylab = "v", xlab = "u", lwd = linewidth, xlim = c(-2.1,0.5), ylim = c(-1.5,1.1), col = 2, xaxt = "n") # Plot u vs v
plot(results.3[,1]$u, results.3[,1]$v, type = "l", ylab = "v", xlab = "u", lwd = linewidth, xlim = c(-2.1,0.5), ylim = c(-1.5,1.1), col = 3) # Plot u vs v
mtext("u", side=1, line=3,las=1, col="black", cex = .75)




