sigma = 0#0.0
mu = 0.0
xi = 0
h = 0.1
N = 10
vol = 0
start.t = 0
end.t = 500
init.q = 0.5
no.iterations = (end.t - start.t) / h
no.replicates = 1000

derivs = function(t, q) {
  sigma * q * (1 - q) + mu * (1 - 2 * q)
}

euler = function(t, q, h) {
  q + h * derivs(t, q)
}

milstein = function(t, q, h) {
  vol = sqrt(abs(q * (1 - q) / N))
  q.pred = q + derivs(t, q) * h + 1.0 * sqrt(h)
  q + derivs(t, q) * h + 1.0 * rnorm(1) * vol * sqrt(h)
}

run.func = function(n, m) {
  qs = vector(mode = "numeric", end.t)
  qs[1] = init.q
  q = qs[1]
  sigma <<- runif(1, min = n, max = m)
  last.t = 0
  for (ii in 2:no.iterations) {
    q = euler(t, q, h)
    q = min(1, q)
    q = max(0, q)
    if (ii * h > last.t) {
      last.t = last.t + 1
      qs[last.t] = q
    }
  }
  return(qs)
}

qs1 = replicate(no.replicates, run.func(0, 1))
qs2 = replicate(no.replicates, run.func(-1, 0))

par(mfrow=c(1,1), mai = c(.5,.8,.5,.5))
plot(
  qs1[, 1],
  type = "l",
  col = "orange",
  bty = "n",
  ylab = "q",
  xlab = "", xaxt = "n", ylim = c(0,1)
)
lines(1:end.t, rep(init.q, end.t - start.t), col = "black", lwd = 2, lty = 2 )
out = apply(qs1[, -1], 2, function(x)
  lines(x, col = "orange"))
out = apply(qs2[, -1], 2, function(x)
  lines(x, col = "forestgreen"))


# print(count / no.replicates)
