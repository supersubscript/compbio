#!/usr/bin/env Rscript
inhibitory = F
EL = -70.0 #mV
ES = ifelse(inhibitory, 0.0, -80.0) #mV
Vth = -54.0 #mV
Vreset = -80.0 # mV
tm = 20.0 # ms
rmgs = 0.15
RmIe = 18.0 # mV
ts = 10.0 # ms
Pmax = 0.5

# Preset
start.p = runif(1)
start.v = -80#runif(1, min = Vreset, max = Vth)
start.z = runif(1)
start.t = 0.0
end.t = 500.0
h = 0.01
no.iterations = 1000.0#(start.t - end.t)/h
vol = c(1,1,999999999)

# Init
p = start.p
z = start.z
t = start.t
v = start.v
h = (end.t - start.t) / no.iterations

z.derivs = function(t, z) {
  -z / ts
}

p.derivs = function(t, p, z) {
  (exp(1) * Pmax * z - p) / ts
}

v.derivs = function(t, v, p) {
  (EL - v - rmgs * p * (v - ES) + RmIe) / tm
}

# Stochastic derivatives
z.derivs = function(t, z) {
  -z / ts
}

p.derivs = function(t, p, z) {
  (exp(1) * Pmax * z - p) / ts
}

v.derivs = function(t, v, p) {
  (EL - v - rmgs * p * (v - ES) + RmIe) / tm
}



milstein = function(t, v, p, z, h, vol) {
  z.pred.first = ifelse(z.derivs(t, z) > 0, sqrt(z.derivs(t, z)), -sqrt(-z.derivs(t, z)))
  z.pred = z + z.derivs(t, z) * h + z.pred.first * sqrt(h)
  z.first = ifelse(z.derivs(t, z) > 0, sqrt(z.derivs(t, z)), -sqrt(-z.derivs(t, z)))
  z.second = ifelse(z.derivs(t, z.pred) > 0, sqrt(z.derivs(t, z.pred)),-sqrt(-z.derivs(t, z.pred)))
  
  p.pred.first = ifelse(p.derivs(t, p, z) > 0, sqrt(p.derivs(t, p, z)),-sqrt(-p.derivs(t, p, z)))
  p.pred = p + p.derivs(t, p, z) * h + p.pred.first * sqrt(h)
  p.first = ifelse(p.derivs(t, p, z) > 0, sqrt(p.derivs(t, p, z)),-sqrt(-p.derivs(t, p, z)))
  p.second = ifelse(p.derivs(t, p, z.pred) > 0, sqrt(p.derivs(t, p, z.pred)),-sqrt(-p.derivs(t, p, z.pred)))
  
  v.vred.first = ifelse(v.derivs(t, v, p) > 0, sqrt(v.derivs(t, v, p)),-sqrt(-v.derivs(t, v, p)))
  v.vred = v + v.derivs(t, v, p) * h + v.vred.first * sqrt(h)
  v.first = ifelse(v.derivs(t, v, p) > 0, sqrt(v.derivs(t, v, p)),-sqrt(-v.derivs(t, v, p)))
  v.second = ifelse(v.derivs(t, v, p.pred) > 0, sqrt(v.derivs(t, v, p.pred)),-sqrt(-v.derivs(t, v, p.pred)))
  
  z <<- 
    z + z.derivs(t, z) * h + z.first * rnorm(1) / vol[1] * sqrt(h) + 1.0 / (2.0 * sqrt(h)) *
    (z.second - z.first) * (rnorm(1) ** 2 / (vol[1]**2) * h - h)
  p <<-
    p + p.derivs(t, p, z) * h + p.first * rnorm(1) / vol[2] * sqrt(h) + 1.0 / (2.0 * sqrt(h)) *
    (p.second - p.first) * (rnorm(1) ** 2 / (vol[2]**2) * h - h)
  v <<-
    v + v.derivs(t, v, p) * h + v.first * rnorm(1) / vol[3] * sqrt(h) + 1.0 / (2.0 * sqrt(h)) *
    (v.second - v.first) * (rnorm(1) ** 2 / (vol[3]**2) * h - h)
  t <<- t + h
  
}

rk4 = function(t, v, p, z, h) {
  z.k1 = z.derivs(t, z)
  p.k1 = p.derivs(t, p, z)
  z.k2 = z.derivs(t + h / 2.0, z + h / 2.0 * z.k1)
  p.k2 = p.derivs(t + h / 2.0, p + h / 2.0 * p.k1, z + h / 2.0 * z.k1)
  z.k3 = z.derivs(t + h / 2.0, z + h / 2.0 * z.k2)
  p.k3 = p.derivs(t + h / 2.0, p + h / 2.0 * p.k2, z + h / 2.0 * z.k2)
  z.k4 = z.derivs(t + h, z + h * z.k3)
  p.k4 = p.derivs(t + h, p + h * p.k3, z + h * z.k3)
  
  z <<- z + h / 6.0 * (z.k1 + 2 * z.k2 + 2 * z.k3 + z.k4)
  p <<- p + h / 6.0 * (p.k1 + 2 * p.k2 + 2 * p.k3 + p.k4)
  t <<- t + h
  for (ii in range(1, length(v))){
    v.k1 = v.derivs(t, v[ii], p)
    v.k2 = v.derivs(t + h / 2.0, v[ii] + h / 2.0 * v.k1, p + h / 2.0 * p.k1)
    v.k3 = v.derivs(t + h / 2.0, v[ii] + h / 2.0 * v.k2, p + h / 2.0 * p.k2)
    v.k4 = v.derivs(t + h, v[ii] + h * v.k3, p + h * p.k3)
    v[(ii %% length(v)) + 1] <<- v[(ii %% length(v)) + 1] + h / 6.0 * (v.k1 + 2 * v.k2 + 2 * v.k3 + v.k4)
  }
}

p.outs = c()
z.outs = c()
v.outs = c()
for (x in 1:no.iterations) {
  #rk4(t, v, p, z, h)
  milstein(t,v,p,z,h, vol)
  if(any(v > Vth)){
    v[which(v > Vth)] = Vreset
    z = 1.0
  }    
  
  #cat(c(t, v, p, z, "\n"), sep = "\t")
  p.outs = append(p.outs, p)
  z.outs = append(z.outs, z)
  v.outs = rbind(v.outs, v)
}

#par(mfrow=c(3,1), mar = c(0,0,0,0), mai = c(.0,.6,.5,.5), oma=c(5,2,2,2))
par(mfrow=c(4,1), oma = c(6,5,4,2) + 0.0,
    mai = c(.0,.5,.0,.5))
plot(z.outs, col = "red", xaxt ="n", xlab = "", type="l", ylab="z")
plot(p.outs, type = "l", xaxt="n", xlab ="", ylab="P")
plot(v.outs[,1], col = "blue", type= "l", ylab="V", xaxt="n")
#plot(v.outs[,2], col = "blue", type= "l", ylab="V")

title(xlab = "Time [ms]",
      ylab = "",
      outer = TRUE, line = 3)