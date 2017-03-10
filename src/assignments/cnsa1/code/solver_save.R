#!/usr/bin/env Rscript
# Set stuff
inhibitory = F
method = "milstein"
EL = -70.0 #mV
ES = ifelse(inhibitory, -80.0, 0.0) #mV
Vth = -54.0 #mV
Vreset = -80.0 # mV
tm = 20.0 # ms
ts = 10.0 # ms
rmgs = 0.15
RmIe = 18.0 # mV
Pmax = 0.5

# Init
start.p = c(0, 0)
start.v = runif(2, min = Vreset, max = Vth)
start.z = c(0, 0)
start.t = 0.0
end.t = 500.0
h = 0.01
no.iterations = (end.t - start.t) / h
vol = c(3, 999999999, 9999999999)

# Init
p = start.p
z = start.z
t = start.t
v = start.v
p.outs = data.frame(x = rep(NA, end.t), y = rep(NA, end.t))
z.outs = data.frame(x = rep(NA, end.t), y = rep(NA, end.t))
v.outs = data.frame(x = rep(NA, end.t), y = rep(NA, end.t))
prev = -1
firing.coord = c()
phase.shift = c()
firing.rate = c()

# Define derivatives
z.derivs = function(t, z) {
  -z / ts
}

p.derivs = function(t, p, z) {
  (exp(1) * Pmax * z - p) / ts
}

v.derivs = function(t, v, p) {
  (EL - v - rmgs * p * (v - ES) + RmIe) / tm
}

euler = function(t, v, p, z, h, vol) {
  z <<- z + z.derivs(t, z) * h
  p <<- p + p.derivs(t, p, z) * h
  v <<- v + v.derivs(t, v, p) * h
  t <<- t + h
}

milstein = function(t, v, p, z, h, vol) {
  for (ii in 1:length(v)) {
    z.pred.first = ifelse(z.derivs(t, z[ii]) > 0, sqrt(z.derivs(t, z[ii])),-sqrt(-z.derivs(t, z[ii])))
    z.pred = z[ii] + z.derivs(t, z[ii]) * h + z.pred.first * sqrt(h)
    z.first = ifelse(z.derivs(t, z[ii]) > 0, sqrt(z.derivs(t, z[ii])),-sqrt(-z.derivs(t, z[ii])))
    z.second = ifelse(z.derivs(t, z.pred) > 0,
                      sqrt(z.derivs(t, z.pred)), -sqrt(-z.derivs(t, z.pred)))
    z[ii] <<-
      z[ii] + z.derivs(t, z[ii]) * h +        z.first * rnorm(1) / vol[1] * sqrt(h) + 1.0 / (2.0 * sqrt(h)) *
      (z.second - z.first) * (rnorm(1) ** 2 / (vol[1] ** 2)) * h
    
    p.pred.first = ifelse(p.derivs(t, p[ii], z[ii]) > 0,
                          sqrt(p.derivs(t, p[ii], z[ii])), -sqrt(-p.derivs(t, p[ii], z[ii])))
    p.pred = p[ii] + p.derivs(t, p[ii], z[ii]) * h + p.pred.first * sqrt(h)
    p.first = ifelse(p.derivs(t, p[ii], z[ii]) > 0,
                     sqrt(p.derivs(t, p[ii], z[ii])), -sqrt(-p.derivs(t, p[ii], z[ii])))
    p.second = ifelse(p.derivs(t, p.pred, z[ii]) > 0,
                      sqrt(p.derivs(t, p.pred, z[ii])),-sqrt(-p.derivs(t, p.pred, z[ii])))
    p[ii] <<-
      p[ii] + p.derivs(t, p[ii], z[ii]) * h + p.first * rnorm(1) / vol[2] * sqrt(h) + 1.0 / (2.0 * sqrt(h)) *
      (p.second - p.first) * (rnorm(1) ** 2 / (vol[2] ** 2)) * h
    
    v.pred.first = ifelse(v.derivs(t, v[ii], p[ii]) > 0,
                          sqrt(v.derivs(t, v[ii], p[ii])), -sqrt(-v.derivs(t, v[ii], p[ii])))
    v.pred = v[ii] + v.derivs(t, v[ii], p[ii]) * h + v.pred.first * sqrt(h)
    v.first = ifelse(v.derivs(t, v[ii], p[ii]) > 0,
                     sqrt(v.derivs(t, v[ii], p[ii])), -sqrt(-v.derivs(t, v[ii], p[ii])))
    v.second = ifelse(v.derivs(t, v.pred, p[ii]) > 0,
                      sqrt(v.derivs(t, v.pred, p[ii])),-sqrt(-v.derivs(t, v.pred, p[ii])))
    v[ii] <<-
      v[ii] + v.derivs(t, v[ii], p[ii]) * h + v.first * rnorm(1) / vol[3] * sqrt(h) + 1.0 / (2.0 * sqrt(h)) *
      (v.second - v.first) * (rnorm(1) ** 2 / (vol[3] ** 2)) * h
  }
  t <<- t + h
}

rk4 = function(t, v, p, z, h, vol) {
  z.k1 = z.derivs(t, z)
  z.k2 = z.derivs(t + h / 2.0, z + h / 2.0 * z.k1)
  z.k3 = z.derivs(t + h / 2.0, z + h / 2.0 * z.k2)
  z.k4 = z.derivs(t + h, z + h * z.k3)
  z <<- z + h / 6.0 * (z.k1 + 2 * z.k2 + 2 * z.k3 + z.k4)
  
  p.k1 = p.derivs(t, p, z)
  p.k2 = p.derivs(t + h / 2.0, p + h / 2.0 * p.k1, z + h / 2.0 * z.k1)
  p.k3 = p.derivs(t + h / 2.0, p + h / 2.0 * p.k2, z + h / 2.0 * z.k2)
  p.k4 = p.derivs(t + h, p + h * p.k3, z + h * z.k3)
  p <<- p + h / 6.0 * (p.k1 + 2 * p.k2 + 2 * p.k3 + p.k4)
  
  v.k1 = v.derivs(t, v, p)
  v.k2 = v.derivs(t + h / 2.0, v + h / 2.0 * v.k1, p + h / 2.0 * p.k1)
  v.k3 = v.derivs(t + h / 2.0, v + h / 2.0 * v.k2, p + h / 2.0 * p.k2)
  v.k4 = v.derivs(t + h, v + h * v.k3, p + h * p.k3)
  v <<- v + h / 6.0 * (v.k1 + 2 * v.k2 + 2 * v.k3 + v.k4)
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

# RUN!
for (x in 1:no.iterations) {
  converging.limit = 10 ** -5
  phase.shift = c(9999999, 999999999)
  # while (abs(phase.shift[1] - phase.shift[2]) > converging.limit) {
  fct(t, v, p, z, h, vol)
  firing = which(v > Vth)
  for (ii in firing) {
    v[ii] = Vreset
    z[ii %% length(z) + 1] = 1.0
    firing.coord = append(firing.coord, t)
    if (ii == 2) {
      phase.shift = append(((firing.coord[length(firing.coord)] - firing.coord[length(firing.coord) - 1]) / (firing.coord[length(firing.coord)] - firing.coord[length(firing.coord) - 2])
      ), phase.shift)
      if (length(firing.coord) > 2)
        firing.rate = append(firing.rate, firing.coord[length(firing.coord)] - firing.coord[length(firing.coord) -
                                                                                              2])
    }
    
    if (phase.shift[1] == 0 & phase.shift[2] > .5)
      phase.shift[1] = 1
    if (phase.shift[1] == 1 & phase.shift[2] < .5)
      phase.shift[1] = 0
  }
  
  #cat(c(t, v[1], v[2], p[1], p[2], z[1], z[2], "\n"), sep = "\t")
  if (floor(t) > prev) {
    prev = floor(t)
    p.outs[t, ] = p
    z.outs[t, ] = z
    v.outs[t, ] = v
  }
}