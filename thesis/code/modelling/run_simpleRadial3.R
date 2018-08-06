#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/modelling/") # Set path to 
source("organism_wrapper.R")

### Set paths
simulator.path    = "/home/henrik/Organism/bin/simulator"
model.path        = "models/simpleRadial3.model"
init.path         = "inits/simpleRadial3.init"
init.ss.path      = "inits/simpleRadial3_ss.init"
det.solver.path   = "solvers/rk4.solver"
stoch.solver.path = "solvers/ms.solver"
model.name        = "simpleRadial3"

### Set simulation parameters
t.end         = 100
no.prints     = 100
step.size     = 1e-2
output.format = 2
error         = 1e-2

###############################################################################
### Define system parameters
###############################################################################
noise = 4

wp = 0.1
wD = 1.2
wd = 0.02

Cp = .2
Ck = .9
Cv = Cp
Cn = 2
Cd = 0.2

cp = 1.2
cD = 2 * wD
cd = 0.1

Xp = 0.1
Xv = 0.4
Xk = 1.
Xn = 4.0
Xd = 0.9

xp = 0.9
xD = 0.002
xd = 0.4

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
model.obj = init.model(model.name)

### Add topology
addTopo(model.obj) = topo.3d()

### Add species
# Anchors
addSpecies(model.obj) = init.species("Sink")
addSpecies(model.obj) = init.species("aWUS")

# With reaction
X     = init.species("X") 
pX    = init.species("pX")
pwus  = init.species("pWUS")
clv3  = init.species("CLV3")
pclv3 = init.species("pCLV3")

# WUS signal
addReaction(pwus) = creationOne(index = model.obj@species$aWUS@index, rate = wp)
addReaction(pwus) = diffusionSimple(rate = wD)
addReaction(pwus) = degradationOne(rate = wd)
addReaction(pwus) = degradationTwo(index = model.obj@species$Sink@index, 
                                       rate = wD)
addSpecies(model.obj) = pwus

# CLV3
addReaction(clv3) = creationOne(index = model.obj@species$pWUS@index, rate = Cp)
addReaction(clv3) = hill(act.ind = list(), 
                         rep.ind = model.obj@species$pWUS@index + 3L, # repression by X 
                         vars = c(Cv, Ck, Cn))
addReaction(clv3) = degradationOne(rate = Cd)
addSpecies(model.obj) = clv3

# CLV3 signal
addReaction(pclv3) = creationOne(index = model.obj@species$pWUS@index, rate = cp)
addReaction(pclv3) = diffusionSimple(rate = cD)
addReaction(pclv3) = degradationOne(rate = cd)
addReaction(pclv3) = degradationTwo(index = model.obj@species$Sink@index, 
                                       rate = cD)
addSpecies(model.obj) = pclv3

# X
addReaction(X) = creationZero(Xp)
addReaction(X) = hill(act.ind = list(), 
                      rep.ind = model.obj@species$pCLV3@index, 
                      vars = c(Xv, Xk, Xn))
addReaction(X) = degradationOne(rate = Xd)
addSpecies(model.obj) = X

# X signal
addReaction(pX) = creationOne(rate = Xp, index = model.obj@species$X@index)
addReaction(pX) = diffusionSimple(rate = xD)
addReaction(pX) = degradationOne(rate = xd)
addReaction(pX) = degradationTwo(index = model.obj@species$Sink@index, rate = xD)
addSpecies(model.obj) = pX

### Add global reactions
# No reactions

### Add neighborhoodRule
addNeighRule(model.obj) = neighborhoodDistanceSphere(index = 3) 

# Write model to file
write.model(model.obj, file.path = model.path)
write.model(model.obj, file.path = stdout())

###############################################################################
### Set up init file
###############################################################################
# Set parameters
no.cells    = 20
no.species  = model.obj@no.species
no.topcells = 3
radius      = 2.5

# Randomize initial values (not really necessary)
init.molecules = matrix(runif(no.species * no.cells), no.cells, no.species)
init           = cbind(seq(from = 0, length.out = no.cells, by = 4), 
                       rep(0, no.cells), rep(0, no.cells), 
                       rep(radius, no.cells), init.molecules)
colnames(init) = c("x", "y", "z", "r", names(model.obj@species))

# Set initial values of sink and anchor
init[, "Sink"]  = rep(0, no.cells); init[1, "Sink"] = 1 # Set sink into stem
init[, "aWUS"]  = rep(0, no.cells); init[no.cells:(no.cells - no.topcells + 1), "aWUS"] = 1 # Set wus domain
init[, "X"]     = rep(0, no.cells); init[1, "X"] = 1 # Set top
init[, "pX"]    = rep(0, no.cells)
init[, "pWUS"]  = rep(0, no.cells)
init[, "CLV3"]  = rep(0, no.cells)
init[, "pCLV3"] = rep(0, no.cells)

# Write to file
write.init(outfile = init.path, init)

###############################################################################
### Simulate system
###############################################################################
result.ss = simulate(simulator.path, model.path, init.path, det.solver.path,      
                     toR = TRUE, output.mode = 2, init.out = init.ss.path)
result    = simulate(simulator.path, model.path, init.ss.path,
                     stoch.solver.path, toR = TRUE,  output.mode = 2)
# 
colnames(result.ss)[1:(model.obj@no.species + 8)] = c("iter", "t", "cell",
                                                   "nNeigh", "x", "y", "z", "r",
                                                   names(model.obj@species))

plot.molecule = function(x) plot(result.ss[, "x"], result.ss[, x], main = x)

par(mfrow = c(3, 2))
plot.molecule("aWUS")
plot.molecule("pWUS")
plot.molecule("CLV3")
plot.molecule("pCLV3")
plot.molecule("X")
plot.molecule("pX")

colnames(result)[1:(model.obj@no.species + 8)] = c("iter", "t", "cell",
                                                   "nNeigh", "x", "y", "z", "r",
                                                   names(model.obj@species))

plot.molecule = function(x) plot(result[, "x"], result[, x], main = x)

par(mfrow = c(3, 2))
plot.molecule("aWUS")
plot.molecule("pWUS")
plot.molecule("CLV3")
plot.molecule("pCLV3")
plot.molecule("X")
plot.molecule("pX")

