#!/usr/bin/env Rscript
setwd("~/compbio/src/assignments/cnsa2/code/")
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 3
##########################
# Actual code
##########################

dvdt.nc = function(v, Ie=-1, ...){
  v*(1-v**2) + Ie
}
dudt.nc = function(v, ...){
  2*v
}

# Simulate
x = seq(-1,1, 0.001)
par(mfrow=c(1,1))
plot(x, dvdt.nc(x, 0), type = "l", col = "darkblue", ylim = c(-1.5,1), bty= "n", lwd = lw.s, xaxt = "n", yaxt="n", xlab = "v", ylab = "u")
axis(1, pos=0)
axis(2, pos=0)
lines(x, dvdt.nc(x, -1), type = "l", col = 2, lwd = lw.s)
lines(x, dudt.nc(x), type = "l", col = 1, lwd = lw.s)
legend("topleft", c("dvdt = 0, Ie = 0", "dvdt = 0, Ie = -1", "dudt = 0"), col = c("darkblue",2,1), inset = c(0.01, 0.05), lty = 1, lwd = lw.s, bg="white", cex = 1.1)

# Do newton to find roots
tolerance = 1e-10
reps = 200
newton = function(f, x0, args) {
  h = 0.001
  i = 2
  x1 = x0
  p = numeric(reps)
  p[1] = x0
  for (ii in 2:reps) {
    dfdx = (f(x0 + h, args) - f(x0, args)) / h
    x1 = (x0 - (f(x0, args) / dfdx))
    p[i] = x1
    i = i + 1
    if (abs(x1 - x0) < tolerance)
      break
    x0 = x1
  }
  return(p[1:(i - 1)])
}
print(newton(function(x, ...) x ** 3 + x + 1, 0, 0))

