#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/modelling/")
library(plyr)
library(tidyverse)
library(gridExtra)
source("organism_wrapper.R")

simulator.path    = "~/oldOrganism/bin/simulator"
model.path        = "models/new_rad1.model"
init.path         = "inits/new_rad1.init"
init.ss.path      = "inits/new_rad1_ss.init"
det.solver.path   = "solvers/rk4.solver"
# stoch.solver.path = "solvers/ms.solver"
# stoch.solver.path = "solvers/rk4.solver"
stoch.solver.path = "solvers/gillespie.solver"
# stoch.solver.path = "solvers/heunito.solver"
model = "spherical1.R"

t.end         = 100
no.prints     = 100
step.size     = 1e-2
output.format = 2
error         = 1e-2
noise         = 1

###############################################################################
### Define reaction parameters
###############################################################################
wp = 1
wD = .9
wd = 0.2
ws = wD

Cv = 150
Ck = 2
Cn = 2
Cl1k = 4
Cl1n = 2 
Cd = 1

l1p = 1000
l1D = 10
l1d = 1
l1s = l1D

###############################################################################
### Setup init file
###############################################################################
init.solver.gillespie("solvers/gillespie.solver", t.start = 0, t.end = t.end, 
                      no.prints = no.prints, step.size = step.size, 
                      output.format = output.format)
init.solver.milsteinStrat("solvers/ms.solver", t.start = 0, t.end = t.end, 
                          no.prints = no.prints, step.size = step.size, 
                          output.format = output.format, noise = noise, vol.flag = 0)
init.solver.HeunItoDiffusive("solvers/heunito.solver", t.start = 0, t.end = t.end, 
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
addSpecies(model.obj) = init.species("L1")

pwus                  = init.species("pWUS")
addReaction(pwus)     = creationOne(index = model.obj@species$aWUS@index, rate = wp)
addReaction(pwus)     = diffusionSimple(rate = wD)
addReaction(pwus)     = degradationOne(rate = wd)
addReaction(pwus)     = degradationTwo(index = model.obj@species$Sink@index, rate = ws)
addSpecies(model.obj) = pwus

clv3                  = init.species("CLV3")
addReaction(clv3)     = hill(act.ind = list(model.obj@species$pWUS@index, 
                                            model.obj@species$pWUS@index + 2L), 
                             rep.ind = list(), vars = c(Cv, Ck, Cn, Cl1k, Cl1n))
addReaction(clv3)     = degradationOne(rate = Cd)
addSpecies(model.obj) = clv3

pl1                  = init.species("pl1")
addReaction(pl1)     = creationOne(index = model.obj@species$L1@index, rate = l1p)
addReaction(pl1)     = diffusionSimple(rate = l1D)
addReaction(pl1)     = degradationOne(rate = l1d)
addReaction(pl1)     = degradationTwo(index = model.obj@species$Sink@index, rate = l1s)
addSpecies(model.obj) = pl1

addNeighRule(model.obj) = neighborhoodDistanceSphere(index = 3) 
write.model(model.obj, file.path = model.path)

###############################################################################
### Set up init file
###############################################################################
# Set parameters
no.cells    = 7
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
init[, "aWUS"] = c(1, .5, 0, 0, rep(0, no.cells - 4)) * 10
init[, "L1"]   = c(rep(0, no.cells - 1), 10)
init[, "pl1"]   = c(rep(0, no.cells - 1), 1)
init[, "pWUS"] = rep(1, no.cells)
init[, "CLV3"] = rep(1, no.cells)
init[, "Sink"] = c(rep(1, no.cells - 1), c(1))

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
# p3 = ggplot(result, aes(x = factor(x/4), y = CLV3)) + geom_violin()
aWUS.aes = aes(x=factor(x/4), y = aWUS)
pWUS.aes = aes(x=factor(x/4), y = pWUS)
CLV3.aes = aes(x=factor(x/4), y = CLV3)
a.aes = aes(x=factor(x/4),    y = pl1)
a2.aes = aes(x=factor(x/4),    y = L1)
p1 = ggplot(result %>% filter(cell < 7)) + geom_line(aWUS.aes) + geom_point(aWUS.aes) + labs(x="", y = "WUS") + theme_bw()
p2 = ggplot(result %>% filter(cell < 7)) + geom_jitter(pWUS.aes, alpha = .1, colour="blue") + geom_violin(pWUS.aes, bw = .7) + labs(x="", y = "wus") + theme_bw()
p3 = ggplot(result %>% filter(cell < 7)) + geom_jitter(CLV3.aes, alpha = .1, colour="red") + geom_violin(CLV3.aes, bw = 6) + 
  labs(x="Distance to apex", y = "CLV3") + theme_bw()
p4 = ggplot(result %>% filter(cell < 7)) + geom_jitter(a.aes, alpha = .1, colour="red") + geom_violin(a.aes, bw = 6) + 
  labs(x="Distance to apex", y = "a") + theme_bw()
p5 = ggplot(result %>% filter(cell < 7)) + geom_jitter(a2.aes, alpha = .1, colour="red") + geom_violin(a2.aes, bw = 6) + 
  labs(x="Distance to apex", y = "a") + theme_bw()

grid.arrange(p1, p2, p3,p4, p5)
