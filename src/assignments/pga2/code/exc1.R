setwd("~/compbio/src/assignments/pga2/code/")
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))
library(scales)

# Parameters / data
N = 100
timepoints = 0:6
data = data.frame(time = 0:6, n = c(28, 39, 72, 83, 98, 96, 99) / N)
x = data[, "time"]
xs = seq(0, 6, 0.1)

# Fraction development over time
q = function(q.init, sigma, t) {
  q.init * exp(sigma * t) / (1 - q.init + q.init * exp(sigma * t))
}

# Calculate the log-likelihood
nas = data[,"n"]*N
prefact = factorial(N)/(sapply(nas, factorial)*sapply(N-nas, factorial))
log.lik = function(q.init, sigma){
  qas = q(q.init, sigma, timepoints)
  # nas = qas*N
  ln.p = log(prefact*qas**(nas)*((1-qas)**(N-nas))) 
  -sum(ln.p)   
}

fit = mle(minuslogl = log.lik, start = list(q.init = .5, sigma = .5), method="L-BFGS-B", lower = 0, upper=c(1, Inf)) 
coef = coef(summary(fit))

# Calculate upper and lower bound of estimate
yhigh = q(coef[1, 1] + coef[1, 2], coef[2, 1] + coef[2, 2], xs)
ylow = q(coef[1, 1] - coef[1, 2], coef[2, 1] - coef[2, 2], xs)

# Plot it all
plot(1, type="n", axes=T, xlab="Time [days]", ylab="Allele frequency", ylim=c(0,100/N), xlim = c(0,6), bty = "n")
polygon(c(xs, rev(xs)), c(yhigh, rev(ylow)),
          col = alpha(2, .33), border = F, lty = 2, lwd = 2)
lines(x, data[,"n"], type = "b", pch = 17, cex = 1.5, lwd = 3, col = 1)
lines(xs, q(coef[1, 1], coef[2, 1], xs), col = 2, lwd = 3, lty = 2)

# Plot borders
lines(xs, ylow, col = alpha(2, 1), lwd = 1)
lines(xs, yhigh, col = alpha(2, 1), lwd = 1)

#################################################################
### Image plot
#################################################################

qvals = seq(0, 1, 0.001)
svals = seq(0, 7.5, 0.05)

likelihoods = sapply(qvals, function(x)
  sapply(svals, function(y)
    log.lik(x, y)))

filled.contour((likelihoods),   xlab = expression(q[0]),  ylab = expression(sigma),
               col = viridis(17),
               axes = F
)
axis(1, at = seq(0, .75, 0.15), labels = round(seq(.1, 1, 1/.75 * .25/2), 2))
axis(2, at = seq(0, 1, 0.2), labels = seq(0, 7.5, 7.5 * 0.2))
axis(4, at = seq(0, .9, 0.15), labels = seq(0, 3000, 500), las =2, line= -2)


################################################
### SIMULATED ANNEALING PART
##################################################
# Parameters
# generations = 300000
# temp        = 1e1
# temp.rate   = 1 - 1e-4
# 
# change.state = function(q.init, sigma){
#   # if(runif(1) > .5){
#     q.init = q.init + rnorm(1, sd = 0.025)
#     q.init = max(q.init, 0)
#     q.init = min(q.init, 1)
#   # }else{
#     sigma = sigma + rnorm(1, sd = 0.025)
#     sigma = max(sigma, 0)
#     sigma = min(sigma, 1)
#   # }
#   list(q=q.init, s=sigma, ll=log.lik(q.init, sigma))
# }
# 
# state = change.state(.5, .5)
# for(iter in 1:generations){
#   # Suggest update
#   new.state = change.state(state[[1]], state[[2]])
#   
#   # Is it better?
#   if (new.state[[3]] <= state[[3]]) {
#     state = new.state
#   }
#   # Should we still go there?
#   else if (runif(1) < exp((state[[3]] - new.state[[3]]) / temp)) {
#     state = new.state
#   }
#   # Change temperature
#   temp = temp * temp.rate
#   
#   # Print current result
#   if (iter %% 10000 == 0){
#     cat(iter, "\t", temp, "\t", state[[1]], "\t", state[[2]], "\t",  state[[3]], "\n")
#   }
# }
