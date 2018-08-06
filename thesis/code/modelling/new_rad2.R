#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/modelling/")
library(tidyverse)
source("organism_wrapper.R")

simulator.path    = "~/Organism/bin/simulator"
model.path        = "models/new_rad2.model"
init.path         = "inits/new_rad2.init"
init.ss.path      = "inits/new_rad2_ss.init"
det.solver.path   = "solvers/rk4.solver"
# stoch.solver.path = "solvers/ms.solver"
# stoch.solver.path = "solvers/rk4.solver"
stoch.solver.path = "solvers/gillespie.solver"
model = "new_rad2"

t.end         = 100
no.prints     = 100
step.size     = 1e-2
output.format = 2
error         = 1e-2
noise         = .05

###############################################################################
### Define reaction parameters
###############################################################################

Wv = 1200
WWk = 4
WWn = 2
Wck = 2
Wcn = 2
Wd = 1

wp = 6
wD = 3
wd = 1
ws = wD

Cv = 150
Cwk = 2
Cwn = 2
Cck = 2
Ccn = 2
Cd = 1
epsilon = 1

cp = 10
cD = 9
cs = cD
cd = 1


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

WUS                  = init.species("WUS")
addReaction(WUS)     = hill(act.ind = list(model.obj@species$aWUS@index), 
                            rep.ind = list(model.obj@species$aWUS@index + 4L), vars = c(Wv, WWk, WWn, Wck, Wcn))
addReaction(WUS)     = degradationOne(rate = Wd)
addReaction(WUS)     = creationZero(rate = epsilon)
addSpecies(model.obj) = WUS


pwus                  = init.species("pWUS")
addReaction(pwus)     = creationOne(index = model.obj@species$WUS@index, rate = wp)
addReaction(pwus)     = diffusionSimple(rate = wD)
addReaction(pwus)     = degradationOne(rate = wd)
addReaction(pwus)     = degradationTwo(index = model.obj@species$Sink@index, rate = ws)
addSpecies(model.obj) = pwus

clv3                  = init.species("CLV3")
addReaction(clv3)     = hill(act.ind = list(model.obj@species$pWUS@index, model.obj@species$pWUS@index + 2L), 
                             rep.ind = list(), vars = c(Cv, Cwk, Cwn, Cck, Ccn))
addReaction(clv3)     = degradationOne(rate = Cd)
addReaction(clv3)     = creationZero(rate = epsilon)
addSpecies(model.obj) = clv3

pclv3                  = init.species("pCLV3")
addReaction(pclv3)     = creationOne(index    = model.obj@species$CLV3@index, rate = cp)
addReaction(pclv3)     = diffusionSimple(rate = cD)
addReaction(pclv3)     = degradationOne(rate  = cd)
addReaction(pclv3)     = degradationTwo(index = model.obj@species$Sink@index, rate = cs)
addSpecies(model.obj)  = pclv3

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
init[, "aWUS"] = c(1, .75, .5, .25, rep(0, no.cells - 4))*10
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

p1 = ggplot(result, aes(x = factor(x/4), y = WUS)) + geom_line() + geom_point() + labs(x="")
p2 = ggplot(result, aes(x = factor(x/4), y = pWUS)) + geom_jitter() + labs(x="")
p3 = ggplot(result, aes(x = factor(x/4), y = CLV3)) + geom_jitter() + labs(x="")
p4 = ggplot(result, aes(x = factor(x/4), y = pCLV3)) + geom_jitter() + labs(x="Radius")
# p2 = ggplot(result, aes(x = factor(x), y = pWUS)) + geom_boxplot()
# p3 = ggplot(result, aes(x = factor(x), y = CLV3)) + geom_boxplot()
grid.arrange(p1, p2, p3, p4)
