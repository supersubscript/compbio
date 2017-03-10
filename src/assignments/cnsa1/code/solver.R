#!/usr/bin/env Rscript
inhibitory = F
method = "milstein"
EL = -70.0 #mV
ESS =  seq(-100, 100, 20) #ifelse(inhibitory, -80.0, 0.0) #mV
ES = -80
Vth = -54.0 #mV
Vreset = -80.0 # mV
tm = 20.0 # ms
ts = 10.0 # ms
rmgs = 0.15
RmIe = 18.0 # mV
Pmax = 0.5

# run.thing = function(aaa) {
#   cat("Now simulating with ES value: ", aaa, "\n")
#   ES = aaa
  # Preset
  start.p = c(0, 0) #runif(2, min = 0, max = Pmax)
  start.v = c(-70, -70)#runif(2, min = Vreset, max = Vth) #c(-80, -70) #runif(2, min = Vreset, max = Vth) #c(-80,-80)
  start.z = c(0, 0) #runif(2)
  start.t = 0.0
  end.t = 500.0
  h = 0.01
  no.iterations = (end.t - start.t) / h
  vol = c(3, 9999999, 99999999)
  
  # Init
  p = start.p
  z = start.z
  t = start.t
  v = start.v
  
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
      z.pred.first = ifelse(z.derivs(t, z[ii]) > 0, sqrt(z.derivs(t, z[ii])), -sqrt(-z.derivs(t, z[ii])))
      z.pred = z[ii] + z.derivs(t, z[ii]) * h + z.pred.first * sqrt(h)
      z.first = ifelse(z.derivs(t, z[ii]) > 0, sqrt(z.derivs(t, z[ii])), -sqrt(-z.derivs(t, z[ii])))
      z.second = ifelse(z.derivs(t, z.pred) > 0,
                        sqrt(z.derivs(t, z.pred)),-sqrt(-z.derivs(t, z.pred)))
      z[ii] <<-
        z[ii] + z.derivs(t, z[ii]) * h +        z.first * rnorm(1) / vol[1] * sqrt(h) + 1.0 / (2.0 * sqrt(h)) *
        (z.second - z.first) * (rnorm(1) ** 2 / (vol[1] ** 2)) * h
      
      p.pred.first = ifelse(p.derivs(t, p[ii], z[ii]) > 0,
                            sqrt(p.derivs(t, p[ii], z[ii])),-sqrt(-p.derivs(t, p[ii], z[ii])))
      p.pred = p[ii] + p.derivs(t, p[ii], z[ii]) * h + p.pred.first * sqrt(h)
      p.first = ifelse(p.derivs(t, p[ii], z[ii]) > 0,
                       sqrt(p.derivs(t, p[ii], z[ii])),-sqrt(-p.derivs(t, p[ii], z[ii])))
      p.second = ifelse(p.derivs(t, p.pred, z[ii]) > 0,
                        sqrt(p.derivs(t, p.pred, z[ii])),
                        -sqrt(-p.derivs(t, p.pred, z[ii])))
      p[ii] <<-
        p[ii] + p.derivs(t, p[ii], z[ii]) * h + p.first * rnorm(1) / vol[2] * sqrt(h) + 1.0 / (2.0 * sqrt(h)) *
        (p.second - p.first) * (rnorm(1) ** 2 / (vol[2] ** 2)) * h
      
      v.pred.first = ifelse(v.derivs(t, v[ii], p[ii]) > 0,
                            sqrt(v.derivs(t, v[ii], p[ii])),-sqrt(-v.derivs(t, v[ii], p[ii])))
      v.pred = v[ii] + v.derivs(t, v[ii], p[ii]) * h + v.pred.first * sqrt(h)
      v.first = ifelse(v.derivs(t, v[ii], p[ii]) > 0,
                       sqrt(v.derivs(t, v[ii], p[ii])),-sqrt(-v.derivs(t, v[ii], p[ii])))
      v.second = ifelse(v.derivs(t, v.pred, p[ii]) > 0,
                        sqrt(v.derivs(t, v.pred, p[ii])),
                        -sqrt(-v.derivs(t, v.pred, p[ii])))
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
  
  p.outs = data.frame(x = rep(NA, end.t), y = rep(NA, end.t))
  z.outs = data.frame(x = rep(NA, end.t), y = rep(NA, end.t))
  v.outs = data.frame(x = rep(NA, end.t), y = rep(NA, end.t))
  prev = -1
  firing.coord = c()
  phase.shift = c()
  firing.rate = c()
  
  if (tolower(method) == "milstein") {
    fct = milstein
  } else if (tolower(method) == "euler") {
    fct = euler
  } else if (tolower(method) == "rk4") {
    fct = rk4
  } else {
    stop("Error: Method not found.")
  }
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
        if (ii == 2){
          phase.shift = append(((firing.coord[length(firing.coord)] - firing.coord[length(firing.coord) - 1]) / (firing.coord[length(firing.coord)] - firing.coord[length(firing.coord) - 2])
          ), phase.shift)
          if(length(firing.coord) > 2)
            firing.rate = append(firing.rate, firing.coord[length(firing.coord)] - firing.coord[length(firing.coord)-2])
        }
        
        if (phase.shift[1] == 0 & phase.shift[2] > .5)
          phase.shift[1] = 1
        if (phase.shift[1] == 1 & phase.shift[2] < .5)
          phase.shift[1] = 0
      }
      
      #cat(c(t, v[1], v[2], p[1], p[2], z[1], z[2], "\n"), sep = "\t")
      if (floor(t) > prev) {
        prev = floor(t)
        p.outs[t,] = p
        z.outs[t,] = z
        v.outs[t,] = v
      }
    }
    
    p.outs = p.outs[!is.na(p.outs[, 1]),]
    z.outs = z.outs[!is.na(z.outs[, 1]),]
    v.outs = v.outs[!is.na(v.outs[, 1]),]
    phase.shift = rev(phase.shift)[-(1:3)]
#    return(list(phase.shift, firing.rate))
# }
# 
# things = lapply(ESS, function(x) run.thing(x))
# firing.rates = sapply(things, function(x) x[[2]][[length(x[[2]])]])
# phase.shifts = sapply(things, function(x) x[[1]])#[[length(x[[1]])]])
# par(mfrow=c(1,1))
# plot(1:length(firing.rates), firing.rates,  ylab = "", ylim=c(0,100), pch = 19, bty = "l", xaxt="n", yaxt = "n")
# lines(1:length(firing.rates), firing.rates, col = "black")


# plot(1, xlab = "", ylab = "", ylim=c(0,1), xlim=c(0,100))
# for(ii in 1:length(phase.shifts)){
#   lines(phase.shifts[[ii]], col = ii)
# }


# firing.rates = lapply(things, function(x) x[[2]][[length(x[[2]])]])
# plot(1:7, firing.rates,  ylab = "", ylim=c(0,30), pch = 19, bty = "l", xaxt="n", yaxt = "n")
# lines(1:7, firing.rates, col = "black")
# axis(1, at=1:7,labels=c(-3, -2, -1, 0, 1, 2, 3), col.axis="black", las=1, cex.axis=1, tck=.03, cex = 2, cex.axis=1.5)
# axis(2, at=c(0, 10,20,30), col.axis="black", las=1, tck=.03, cex = 2, cex.axis=1.5)
# mtext(expression(tau[s]/tau[m]), side=1, line=3,las=1, col="black", cex = 1.5)
# mtext("Period time [ms]", side=2, line=3,las=3, col="black", cex = 1.5)


  # things = replicate(300, run.thing())
  # save(things, file = "inhib.Rdata")
  # dev.off()
  # par(mfrow = c(1, 1))
  # #plot(1,type='n',xlim=c(1,10), ylim=c(0,max_y),xlab='ID', ylab='Frequency')
  # plot(
  #   things[[1]],
  #   xlim = c(1, 100),
  #   ylim = c(.0, 1),
  #   type = "l",
  #   xlab = "Iterations",
  #   ylab = "Firing events"
  # )
  # for (ii in things[-1]) {
  #   lines(ii, col = sample(rainbow(10)))
  # }
  
  par(mfrow=c(3,1), oma = c(6,5,4,2) + 0.0,
      mai = c(.0,.6,.0,.5))
  plot(z.outs[,1], col = "red", xaxt ="n", xlab = "", type="l", ylab="z", lwd= 2, cex.lab = 2)
  lines(z.outs[,2], col = "red", lty=2, lwd= 2)
  plot(p.outs[,1], type = "l", xaxt="n", xlab="", ylab="P", lwd= 2, ylim=c(0,.8), cex.lab = 2)
  lines(p.outs[,2], col = "black", lty=2, lwd = 2)
  plot(v.outs[,1], col = "blue", type= "l", ylab="V", xaxt="n", lwd=2, cex.lab = 2)
  lines(v.outs[,2], col = "blue", lty=2, lwd=2)

  # plot(v.outs[,2], col = "blue", type= "l", ylab="V")
  # par(mfrow=c(1,1))
  # plot(v.outs[,1], v.outs[,2], type = "l")
  #
  # title(xlab = "Time [ms]",
  #       ylab = "",
  #       outer = TRUE, line = 3)