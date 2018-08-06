#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/modelling/")
source("organism_wrapper.R")

simulator.path = "~/compbio/thesis/software/organism/bin/simulator"
model.path     = "models/simpleRadial1.model"
init.path      = "inits/simpleRadial1.init"
init.ss.path   = "inits/simpleRadial1_ss.init"
solver.path    = "solvers/ms.solver"

###############################################################################
### Define reaction parameters
###############################################################################
wp = .5
wD = .3
wd = .1
ws = wD # Sink decay same as diffusion
Cp = .4
Cd = .2 


###############################################################################
### Setup init file
###############################################################################
init.solver.milsteinStrat("solvers/ms.solver", t.start = 0, t.end = 100, no.prints = 100, step.size = 1e-3, output.format = 2, noise = 10, vol.flag = 1)
init.solver.rk5adaptive("solvers/rk5.solver",  t.start = 0, t.end = 100, no.prints = 100, step.size = 1e-3, output.format = 2, error = 1e-2)
init.solver.rk4("solvers/rk4.solver",          t.start = 0, t.end = 100, no.prints = 100, step.size = 1e-2, output.format = 2)

###############################################################################
### Setup model file
###############################################################################
source("models.R")
model.obj = init.model("simpleRadial1")

### Add topology
addTopo(model.obj) = topo.3d()

### Add species
addSpecies(model.obj) = init.species("Sink")
addSpecies(model.obj) = init.species("aWUS")

# wus                   = init.species("WUS")
# addReaction(wus)      = creationOne(index = model.obj@species$aWUS@index, rate = Wp)
# addReaction(wus)      = degradationOne(rate = Wd)
# addSpecies(model.obj) = wus

pwus                  = init.species("pWUS")
addReaction(pwus)     = creationOne(index = model.obj@species$aWUS@index, rate = wp)
addReaction(pwus)     = diffusionSimple(rate = wD)
addReaction(pwus)     = degradationOne(rate = wd)
addReaction(pwus)     = degradationTwo(index = model.obj@species$Sink@index, rate = ws)
addSpecies(model.obj) = pwus

clv3                  = init.species("CLV3")
addReaction(clv3)     = creationOne(index = model.obj@species$pWUS@index, rate = Cp)
addReaction(clv3)     = degradationOne(rate = Cd)
addSpecies(model.obj) = clv3

### Add reactions
# No reaction

### Add neighborhoodRule
addNeighRule(model.obj) = neighborhoodDistanceSphere(index = 3) 

# Write model to file

write.model(model.obj, file.path = model.path)

###############################################################################
### Set up init file
###############################################################################
# Set parameters
no.cells    = 20
no.species  = model.obj@no.species
no.topcells = 1
radius      = 2.5

# Randomize initial values
init.molecules = matrix(runif(no.species * no.cells), no.cells, no.species)
init           = cbind(seq(from = 0, length.out = no.cells, by = 4), 
                       rep(0, no.cells), rep(0, no.cells), rep(radius, no.cells), 
                       init.molecules)
colnames(init) = c("x", "y", "z", "r", names(model.obj@species))

# Set initial values of sink and anchor
init[, "aWUS"] = c(rep(1, no.topcells), rep(0, no.cells - no.topcells))
init[, "pWUS"] = rep(0, no.cells)
init[, "CLV3"] = rep(0, no.cells)
init[, "Sink"] = c(rep(0, no.cells - 1), 1)

# Write to file
write.init(outfile = init.path, init)

###############################################################################
### Simulate system
###############################################################################
result.ss = simulate(simulator.path, model.path, init.path, "solvers/rk4.solver", toR = TRUE, output.mode = 2, init.out = "sim")
result = simulate(simulator.path, model.path, init.path, solver.path, toR = TRUE, output.mode = 2)

colnames(result)[1:(model.obj@no.species + 8)] = c("iter", "t", "cell", "nNeigh", "x", "y", "z", "r", names(model.obj@species))
head(result)
par(mfrow = c(4, 1))
plot(result[result[, "iter"] == 99, "cell"], result[result[, "iter"] == 99, "pWUS"], type = "l", main =
       "pWUS")
plot(result[, "t"], result[, "pWUS"], type = "p", main =
       "pWUS / t")
print(result[result[, "iter"] == 99, "pWUS"])
plot(result[result[, "iter"] == 99, "cell"], result[result[, "iter"] == 99, "CLV3"], type = "l", main = "CLV3")
plot(result[, "t"], result[, "CLV3"], main = "CLV3")
