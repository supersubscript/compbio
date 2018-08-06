#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/modelling/")
library(tidyverse)
library(gridExtra)
source("organism_wrapper.R")

simulator.path    = "~/Organism/bin/simulator"
model.path        = "models/new_rad1.model"
init.path         = "inits/new_rad1.init"
init.ss.path      = "inits/new_rad1_ss.init"
det.solver.path   = "solvers/rk4.solver"
# stoch.solver.path = "solvers/ms.solver"
# stoch.solver.path = "solvers/rk4.solver"
stoch.solver.path = "solvers/gillespie.solver"
model = "new_rad1"

t.end         = 100
no.prints     = 100
step.size     = 1e-2
output.format = 2
error         = 1e-2
noise         = .6

###############################################################################
### Define reaction parameters
###############################################################################
wp = 1
wD = .15
wd = 0.2
ws = wD

Cv = 150
Ck = 4
Cn = 2
Cd = 1

compute.contacts = function(x, y, z, r) {
  nCells = length(x)
  contacts = matrix(0, nrow = nCells, ncol = nCells)
  for (ii in 1:nCells)
    for (jj in 1:nCells) {
      if (sqrt((x[ii] - x[jj]) ** 2 + (y[ii] - y[jj]) ** 2 + (z[ii] - z[jj]) **
               2) <= r[ii] + r[jj]) {
        contacts[ii, jj] = 1
        contacts[jj, ii] = 1
      }
    }
  contacts
}

###########
# computes the stable state of column undergoing production, diffusion, degradation and sink degradation (= diffusion)
# rates: p: production, g: degradation, d: diffusion
# cp: productor vector, s: sink vector
pdds.abstract = function(params, cp, s, contacts) {
  # params is on p, D, d format
  rSum = rowSums(contacts)
  mat  = contacts * params[2]
  vect = -cp * params[1]
  diag(mat) = -params[2] * (s + rSum) - params[2]
  solve(mat, vect)
}

# example:
# cp = rep(0, length())
# params = c(1,2,.3)
# cp = init[,"aWUS"]
# s = init[,"Sink"]
# contacts = compute.contacts(init[, 1], init[, 2], init[, 3], init[, 4])

###############################################################################
### Setup init file
###############################################################################
init.solver.gillespie("solvers/gillespie.solver", t.start = 0, t.end = t.end, 
                      no.prints = no.prints, step.size = step.size, 
                      output.format = output.format)
init.solver.milsteinStrat("solvers/ms.solver", t.start = 0, t.end = t.end, 
                          no.prints = no.prints, step.size = step.size, 
                          output.format = output.format, noise = noise, vol.flag = 0)

###############################################################################
### Setup model file
###############################################################################
source("models.R")
model.obj = init.model(model)
addTopo(model.obj)    = topo.3d()
addSpecies(model.obj) = init.species("Sink")
addSpecies(model.obj) = init.species("aWUS")

l1                   = init.species("l1")
addReaction(l1)     = creationOne(index = model.obj@species$L1@index, rate = l1p)
addReaction(l1)     = diffusionSimple(rate = l1D)
addReaction(l1)     = degradationOne(rate = l1d)
addReaction(l1)     = degradationTwo(index = model.obj@species$Sink@index, rate = l1s)
addSpecies(model.obj) = l1

l2                    = init.species("l2")
addReaction(l2)       = creationOne(index = model.obj@species$aWUS@index, rate = l2p)
addReaction(l2)       = diffusionSimple(rate = l2D)
addReaction(l2)       = degradationOne(rate = l2d)
addReaction(l2)       = degradationTwo(index = model.obj@species$Sink@index, rate = l2s)
addSpecies(model.obj) = l2


clv3                  = init.species("CLV3")
addReaction(clv3)     = hill(act.ind = model.obj@species$pWUS@index, rep.ind = list(), vars = c(Cv, Ck, Cn))
addReaction(clv3)     = degradationOne(rate = Cd)
addSpecies(model.obj) = clv3

addNeighRule(model.obj) = neighborhoodDistanceSphere(index = 3) 
write.model(model.obj, file.path = model.path)

###############################################################################
### Set up init file
###############################################################################
# Set parameters
no.cells    = 8
no.species  = model.obj@no.species
no.topcells = 4
radius      = 2.5

# Randomize initial values
init.molecules = matrix(runif(no.species * no.cells), no.cells, no.species)
init           = cbind(seq(from = 0, length.out = no.cells, by = 4), 
                       rep(0, no.cells), rep(0, no.cells), rep(radius, no.cells), 
                       init.molecules)
colnames(init) = c("x", "y", "z", "r", names(model.obj@species))

# Set initial values of sink and anchor
init[, "aWUS"] = c(1, .5, 0, 0, rep(0, no.cells - 4))*10
init[, "pWUS"] = rep(0, no.cells)
init[, "CLV3"] = rep(0, no.cells)
init[, "Sink"] = c(rep(0, no.cells - 1), c(1))

# Write to file
write.init(outfile = init.path, init)

###############################################################################
### Simulate system
###############################################################################
result.ss = simulate(
  simulator.path,
  model.path,
  init.path,
  det.solver.path,
  toR         = FALSE,
  output.mode = 2,
  init.out    = init.ss.path
)
result = simulate(
  simulator.path,
  model.path,
  init.ss.path,
  stoch.solver.path,
  toR         = TRUE,
  output.mode = 2
)

colnames(result)[1:(model.obj@no.species + 8)] = c("iter", "t", "cell", "nNeigh", "x", "y", "z", "r", names(model.obj@species))

# p1 = ggplot(result, aes(x = factor(x/4), y = aWUS)) + geom_line() + geom_point()
# p2 = ggplot(result, aes(x = factor(x/4), y = pWUS)) + geom_violin()
# # p3 = ggplot(result, aes(x = factor(x/4), y = CLV3)) + geom_violin()
# aWUS.aes = aes(x=factor(x/4), y = aWUS)
# pWUS.aes = aes(x=factor(x/4), y = pWUS)
# CLV3.aes = aes(x=factor(x/4), y = CLV3)
# p1 = ggplot(result %>% filter(cell < 7)) + geom_line(aWUS.aes) + geom_point(aWUS.aes) + labs(x="", y = "WUS") + theme_bw()
# p2 = ggplot(result %>% filter(cell < 7)) + geom_jitter(pWUS.aes, alpha = .1, colour="blue") + geom_violin(pWUS.aes, bw = .7) + labs(x="", y = "wus") + theme_bw()
# p3 = ggplot(result %>% filter(cell < 7)) + geom_jitter(CLV3.aes, alpha = .1, colour="red") + geom_violin(CLV3.aes, bw = 6) + 
#   labs(x="Distance to apex", y = "CLV3") + theme_bw()
grid.arrange(p1, p2, p3)
# p2 = ggplot(result, aes(x = factor(x), y = pWUS)) + geom_boxplot()
# p3 = ggplot(result, aes(x = factor(x), y = CLV3)) + geom_boxplot()
