#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/modelling/")
library(tidyverse)
source("organism_wrapper.R")

simulator.path    = "~/Organism/bin/simulator"
model.path        = "models/simpleRadial4.model"
init.path         = "inits/simpleRadial4.init"
init.ss.path      = "inits/simpleRadial4_ss.init"
det.solver.path   = "solvers/rk4.solver"
# stoch.solver.path = "solvers/ms.solver"
stoch.solver.path = "solvers/gillespie.solver"
model.name        = "simple4"

t.end         = 100
no.prints     = 100
step.size     = 1e-3
output.format = 2
error         = 1e-2
noise         = 1

###############################################################################
### Define reaction parameters
###############################################################################
wp = 1.0
wD = 0.3
wd = 0.2
ws = wD

ap = .3
aD = .2
ad = .1

Cv  = 20
Cwk = 19
Cwn = 2.0
Cak = .5
Can = 2.0
Cd  = 0.2

cp = 2
cD = 2.0 * wD
cd = 0.1

###############################################################################
### Setup init file
###############################################################################
init.solver.gillespie("solvers/gillespie.solver", t.start = 0, t.end = t.end, 
                      no.prints = no.prints, step.size = step.size, 
                      output.format = output.format)

###############################################################################
### Setup model file
###############################################################################
source("models.R")
model.obj             = init.model(model.name)
addTopo(model.obj)    = topo.3d()
addSpecies(model.obj) = init.species("Sink")
addSpecies(model.obj) = init.species("aWUS")
addSpecies(model.obj) = init.species("aA")

pwus                  = init.species("pWUS")
addReaction(pwus)     = creationOne(index = model.obj@species$aWUS@index, rate = wp)
addReaction(pwus)     = diffusionSimple(rate = wD)
addReaction(pwus)     = degradationOne(rate = wd)
addReaction(pwus)     = degradationTwo(index = model.obj@species$Sink@index, rate = ws)
addSpecies(model.obj) = pwus

pA                    = init.species("pA")
addReaction(pA)       = creationOne(index = model.obj@species$aA@index, rate = ap)
addReaction(pA)       = diffusionSimple(rate = aD)
addReaction(pA)       = degradationOne(rate = ad)
addReaction(pA)       = degradationTwo(index = model.obj@species$Sink@index, rate = ws)
addSpecies(model.obj) = pA

clv3                  = init.species("CLV3")
addReaction(clv3)     = hill(act.ind = c(model.obj@species$pWUS@index, model.obj@species$pA@index), 
                             rep.ind = list(), vars = c(Cv, Cwk, Cwn, Cak, Can))
addReaction(clv3)     = degradationOne(rate = Cd)
addSpecies(model.obj) = clv3

pclv3                 = init.species("pCLV3")
addReaction(pclv3)    = creationOne(index = model.obj@species$CLV3@index, rate = cp)
addReaction(pclv3)    = diffusionSimple(rate = cD)
addReaction(pclv3)    = degradationOne(rate = cd)
addReaction(pclv3)    = degradationTwo(index = model.obj@species$Sink@index, rate = ws)
addSpecies(model.obj) = pclv3


addNeighRule(model.obj) = neighborhoodDistanceSphere(index = 3) 
write.model(model.obj, file.path = model.path)

###############################################################################
### Set up init file
###############################################################################
# Set parameters
no.cells    = 12
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
init[, "aWUS"]  = c(1, .75, .5, .25, rep(0, no.cells - 4))*10
init[, "aA"]    = c(1, 0, 0, 0, rep(0, no.cells - 4))*10
init[, "pWUS"]  = rep(0, no.cells)
init[, "CLV3"]  = rep(0, no.cells)
init[, "pCLV3"] = rep(1, no.cells)
init[, "pA"]    = rep(1, no.cells)
init[, "Sink"]  = c(rep(0, no.cells - 1), c(1))

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
) %>% as.tibble

result[result < 0] = 0

colnames(result)[1:(model.obj@no.species + 8)] = c("iter", "t", "cell", "nNeigh", "x", "y", "z", "r", names(model.obj@species))

p1 = ggplot(result, aes(x = factor(x), y = aWUS))  + geom_line() + geom_point()
p2 = ggplot(result, aes(x = factor(x), y = pWUS))  + geom_point()
p3 = ggplot(result, aes(x = factor(x), y = CLV3))  + geom_point()
p4 = ggplot(result, aes(x = factor(x), y = pCLV3)) + geom_point()
p5 = ggplot(result, aes(x = factor(x), y = pA))    + geom_point()
grid.arrange(p1, p2, p3, p4, p5)
