#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/modelling/")
source("organism_wrapper.R")

simulator.path    = "/home/henrik/Organism/bin/simulator"
model.path        = "models/simpleRadial2_center.model"
init.path         = "inits/simpleRadial2_center.init"
init.ss.path      = "inits/simpleRadial2_center_ss.init"
det.solver.path   = "solvers/rk4.solver"
stoch.solver.path = "solvers/ms.solver"

t.end         = 100
no.prints     = 100
step.size     = 1e-2
output.format = 2
error         = 1e-2

###############################################################################
### Define reaction parameters
###############################################################################
noise = 20

wp = 0.2
wD = 2
wd = 0.1

Cp = 0.1
Cv = 0.0
Ck = 0.1
Cn = 2.0
Cd = 0.1

cp = 0.5
cD = 2 * wD
cd = 0.2

###############################################################################
### Setup init file
###############################################################################
init.solver.milsteinStrat("solvers/ms.solver", t.start = 0, t.end = t.end, 
                          no.prints = no.prints, step.size = step.size, 
                          output.format = output.format, noise = noise, 
                          vol.flag = 1)
init.solver.rk4("solvers/rk4.solver", t.start = 0, t.end = t.end, 
                no.prints = no.prints, step.size = step.size, 
                output.format = output.format)

###############################################################################
### Setup model file
###############################################################################
source("models.R")
model.obj = init.model("simpleRadial2")

### Add topology
addTopo(model.obj) = topo.3d()

### Add species
addSpecies(model.obj) = init.species("Sink")
addSpecies(model.obj) = init.species("aWUS")
addSpecies(model.obj) = init.species("aCLV3")


pwus                  = init.species("pWUS")
addReaction(pwus)     = creationOne(index = model.obj@species$aWUS@index, rate = wp)
addReaction(pwus)     = diffusionSimple(rate = wD)
addReaction(pwus)     = degradationOne(rate = wd)
addReaction(pwus)     = degradationTwo(index = model.obj@species$Sink@index, rate = wD)
addSpecies(model.obj) = pwus

clv3                  = init.species("CLV3")
addReaction(clv3)     = creationOne(rate = Cp)
addReaction(clv3)     = hill(act.ind = model.obj@species$pWUS@index + 2L, rep.ind = list(), vars = c(Cv, Ck, Cn))
addReaction(clv3)     = degradationOne(rate = Cd)
addSpecies(model.obj) = clv3

pclv3                 = init.species("pCLV3")
addReaction(pclv3)    = creationOne(index = model.obj@species$CLV3@index, rate = cp)
addReaction(pclv3)    = diffusionSimple(rate = cD)
addReaction(pclv3)    = degradationOne(rate = cd)
addReaction(pwus)     = degradationTwo(index = model.obj@species$Sink@index, rate = cD)
addSpecies(model.obj) = pclv3

### Add reactions
# No reaction

### Add neighborhoodRule
addNeighRule(model.obj) = neighborhoodDistanceSphere(index = 3) 

# Write model to file
write.model(model.obj, file.path = model.path)
write.model(model.obj, file.path = stdout())

###############################################################################
### Set up init file
###############################################################################
# Set parameters
no.cells    = 21
no.species  = model.obj@no.species
no.topcells = 3
radius      = 2.5

# Randomize initial values
init.molecules = matrix(runif(no.species * no.cells), no.cells, no.species)
init           = cbind(seq(from = 0, length.out = no.cells, by = 4), 
                       rep(0, no.cells), rep(0, no.cells), 
                       rep(radius, no.cells), init.molecules)
colnames(init) = c("x", "y", "z", "r", names(model.obj@species))

# Set initial values of sink and anchor
init[, "aWUS"]  = rep(0, no.cells); init[11, "aWUS"] = 1
init[, "pWUS"]  = rep(0, no.cells)
init[, "CLV3"]  = rep(0, no.cells); init[21, "CLV3"] = 1
init[, "pCLV3"] = rep(0, no.cells)
init[, "Sink"]  = rep(0, no.cells); init[1, "Sink"] = 1

# Write to file
write.init(outfile = init.path, init)

###############################################################################
### Simulate system
###############################################################################
result.ss = simulate(simulator.path, model.path, init.path, det.solver.path,      
                     toR = TRUE, output.mode = 2, init.out = init.ss.path)
result    = simulate(simulator.path, model.path, init.ss.path,
                     stoch.solver.path, toR = TRUE,  output.mode = 2)
colnames(result)[1:(model.obj@no.species + 8)] = c("iter", "t", "cell",
                                                   "nNeigh", "x", "y", "z", "r",
                                                   names(model.obj@species))

plot.molecule = function(x) plot(result[, "x"], result[, x], main = x)

par(mfrow = c(2, 2))
plot.molecule("aWUS")
plot.molecule("pWUS")
plot.molecule("CLV3")
plot.molecule("pCLV3")

