#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/modelling/")
library(tidyverse)
source("organism_wrapper.R")

simulator.path    = "~/Organism/bin/simulator"
model.path        = "models/simpleRadial1.model"
init.path         = "inits/simpleRadial1.init"
init.ss.path      = "inits/simpleRadial1_ss.init"
det.solver.path   = "solvers/rk4.solver"
# stoch.solver.path = "solvers/ms.solver"
stoch.solver.path = "solvers/gillespie.solver"

t.end         = 100
no.prints     = 100
step.size     = 1e-3
output.format = 2
error         = 1e-2
noise         = 1

###############################################################################
### Define system targets
###############################################################################
load("../modelling/binned_d2t_12.RData")
binned.data = binned.data[, 1:3]
binned.data = binned.data[, -1]
bin.targets = unname(as.matrix(binned.data))
sum.targets = unname(colSums(binned.data))

###############################################################################
### Define reaction parameters
###############################################################################
wp = 2
wD = 0.4
wd = 0.2
ws = wD

Cv = 2
Ck = 1 
Cn = 2
Cd = 0.9

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
model.obj = init.model("simpleRadial1")
addTopo(model.obj)    = topo.3d()
addSpecies(model.obj) = init.species("Sink")
addSpecies(model.obj) = init.species("aWUS")

pwus                  = init.species("pWUS")
addReaction(pwus)     = creationOne(index = model.obj@species$aWUS@index, rate = wp)
addReaction(pwus)     = diffusionSimple(rate = wD)
addReaction(pwus)     = degradationOne(rate = wd)
addReaction(pwus)     = degradationTwo(index = model.obj@species$Sink@index, rate = ws)
addSpecies(model.obj) = pwus

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

p1 = ggplot(result, aes(x = x, y = aWUS)) + geom_line() + geom_point()
p2 = ggplot(result, aes(x = factor(x), y = pWUS)) + geom_point()
p3 = ggplot(result, aes(x = factor(x), y = CLV3)) + geom_point()
# p2 = ggplot(result, aes(x = factor(x), y = pWUS)) + geom_boxplot()
# p3 = ggplot(result, aes(x = factor(x), y = CLV3)) + geom_boxplot()
grid.arrange(p1, p2, p3)
