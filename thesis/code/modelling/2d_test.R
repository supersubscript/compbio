


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
model = "2d_test"

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

###
init = read.init("../modelling/inits/2d_template.init")
colnames(init)[5:ncol(init)] = c("L1", "Sink", "WUS", "CLV3", "KAN")
init = init[,-ncol(init)]

init.solver.gillespie("solvers/gillespie.solver", t.start = 0, t.end = t.end, 
                      no.prints = no.prints, step.size = step.size, 
                      output.format = output.format)
init.solver.milsteinStrat("solvers/ms.solver", t.start = 0, t.end = t.end, 
                          no.prints = no.prints, step.size = step.size, 
                          output.format = output.format, noise = noise, vol.flag = 0)
###
source("models.R")
model.obj = init.model(model)
addTopo(model.obj)    = topo.3d()
addSpecies(model.obj) = init.species("Sink")
addSpecies(model.obj) = init.species("WUS")

l1                   = init.species("L1")
addReaction(l1)     = creationOne(index = model.obj@species$L1@index, rate = l1p)
addReaction(l1)     = diffusionSimple(rate = l1D)
addReaction(l1)     = degradationOne(rate = l1d)
addReaction(l1)     = degradationTwo(index = model.obj@species$Sink@index, rate = l1s)
addSpecies(model.obj) = l1

clv3                  = init.species("aCLV3")
addReaction(clv3)     = hill(act.ind = model.obj@species$CLV3@index, rep.ind = list(), vars = c(Cv, Ck, Cn))
addReaction(clv3)     = degradationOne(rate = Cd)
addSpecies(model.obj) = clv3

addNeighRule(model.obj) = neighborhoodDistanceSphere(index = 3) 
write.model(model.obj, file.path = model.path)

