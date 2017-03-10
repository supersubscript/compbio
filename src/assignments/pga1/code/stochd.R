sigma = 0
mu = 0.0
h = 0.1
N = 10000
vol = 0
start.t = 0
end.t = 500
init.q = 0.2
no.iterations = (end.t - start.t) / h
no.replicates = 500

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

run.func = function() {
  qs = vector(mode = "numeric", end.t - start.t)
  qs[1] = init.q
  q = qs[1]
  last.t = 0
  for (ii in 2:no.iterations) {
    q = milstein(t, q, h)
    q = min(1, q)
    q = max(0, q)
    if (ii * h > last.t) {
      last.t = last.t + 1
      qs[last.t] = q
    }
  }
  return(qs)
}

qs1 = replicate(no.replicates, run.func())
qs2 = replicate(no.replicates, run.func())
qs3 = replicate(no.replicates, run.func())

par(mfrow=c(1,1), mai = c(.5,.8,.5,.5))
plot(
  qs1[, 1],
  type = "l",
  col = sample(rainbow(10)),
  bty = "n",
  ylab = "q",
  xlab = "", xaxt = "n", ylim = c(0,1)
)
# lines(1:end.t, rep(init.q, end.t - start.t), col = "black", lwd = 2, lty = 2 )
out = apply(qs1[, -1], 2, function(x)
  lines(x, col = sample(rainbow(10))))
