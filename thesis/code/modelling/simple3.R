#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/modelling/")
library(tidyverse)
source("organism_wrapper.R")

simulator.path    = "~/Organism/bin/simulator"
model.path        = "models/simpleRadial3.model"
init.path         = "inits/simpleRadial3.init"
init.ss.path      = "inits/simpleRadial3_ss.init"
det.solver.path   = "solvers/rk4.solver"
# stoch.solver.path = "solvers/ms.solver"
stoch.solver.path = "solvers/gillespie.solver"
model.name        = "simple3"

t.end         = 100
no.prints     = 100
step.size     = 1e-3
output.format = 2
error         = 1e-2
noise         = 1

###############################################################################
### Define reaction parameters
###############################################################################
wp = 0.5
wD = 1.2
wd = 0.2
ws = wD

xp = 1
xD = .1
xd = .4
xs = xD

Cv = 20
Cwk = 5
Cwn = 2
Cxk = 10000000000
Cxn = 8
Cd = 0.8

cp = 1.9
cD = 2 * wD
cd = 0.1
cs = cD


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
addSpecies(model.obj) = init.species("aX")

pwus                  = init.species("pWUS")
addReaction(pwus)     = creationOne(index = model.obj@species$aWUS@index, rate = wp)
addReaction(pwus)     = diffusionSimple(rate = wD)
addReaction(pwus)     = degradationOne(rate = wd)
addReaction(pwus)     = degradationTwo(index = model.obj@species$Sink@index, rate = ws)
addSpecies(model.obj) = pwus

pX                    = init.species("pX")
addReaction(pX)       = creationOne(index = model.obj@species$aX@index, rate = xp)
addReaction(pX)       = diffusionSimple(rate = xD)
addReaction(pX)       = degradationOne(rate = xd)
addReaction(pX)       = degradationTwo(index = model.obj@species$Sink@index, rate = xs)
addSpecies(model.obj) = pX

clv3                  = init.species("CLV3")
addReaction(clv3)     = hill(act.ind = c(model.obj@species$pWUS@index),
                             rep.ind = c(model.obj@species$pX@index), vars = c(Cv, Cwk, Cwn, Cxk, Cxn))
addReaction(clv3)     = degradationOne(rate = Cd)
addSpecies(model.obj) = clv3

pclv3                 = init.species("pCLV3")
addReaction(pclv3)    = creationOne(index = model.obj@species$CLV3@index, rate = cp)
addReaction(pclv3)    = diffusionSimple(rate = cD)
addReaction(pclv3)    = degradationOne(rate = cd)
addReaction(pclv3)    = degradationTwo(index = model.obj@species$Sink@index, rate = cs)
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
radius      = 3

# Randomize initial values
init.molecules = matrix(runif(no.species * no.cells), no.cells, no.species)
init           = cbind(seq(from = 0, length.out = no.cells, by = 4), 
                       rep(0, no.cells), rep(0, no.cells), rep(radius, no.cells), 
                       init.molecules)
colnames(init) = c("x", "y", "z", "r", names(model.obj@species))

# Set initial values of sink and anchor
init[, "aWUS"]  = c(10, 8, 6, 4, 2, rep(0, no.cells - 5))
init[, "aX"]    = c(1, 0, 0, 0, rep(0, no.cells - 4))*10
init[, "pWUS"]  = rep(0, no.cells)
init[, "CLV3"]  = rep(0, no.cells)
init[, "pCLV3"] = rep(1, no.cells)
init[, "pX"]    = rep(1, no.cells)
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
  toR         = TRUE,
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

p1 = ggplot(result, aes(x = x, y = aWUS)) + geom_line() + geom_point() + labs(x = "Distance from top")
p2 = ggplot(result, aes(x = factor(x), y = pWUS))  + geom_point() + labs(x = "Distance from top")
p3 = ggplot(result, aes(x = factor(x), y = CLV3))  + geom_point() + labs(x = "Distance from top")
p4 = ggplot(result, aes(x = factor(x), y = pCLV3)) + geom_point() + labs(x = "Distance from top")
# p5 = ggplot(result, aes(x = factor(x), y = pX))    + geom_point()
grid.arrange(p1, p2, p3, p4, p5)
