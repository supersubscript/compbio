#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/modelling/")
library(tidyverse)
library(gridExtra)
source("organism_wrapper.R")

simulator.path    = "~/Organism/bin/simulator"
model.path        = "models/new_rad1.model"
init.path         = "inits/new_rad1.init"
init.ss.path      = "inits/new_rad1_ss.init"
# det.solver.path   = "solvers/rk4.solver"
stoch.solver.path = "solvers/ms.solver"
# stoch.solver.path = "solvers/rk4.solver"
# stoch.solver.path = "solvers/gillespie.solver"
model = "new_rad1"

t.end         = 100
no.prints     = 100
step.size     = 1e-2
output.format = 2
error         = 1e-2
noise         = .5


###############################################################################
### Define reaction parameters
###############################################################################
wp = 5
wD = 5
wd = 1
ws = wD

Cv = 150
Ck = 5
Cn = 4
Cd = 1


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
init[, "aWUS"] = c(1, 0, 0, 0, rep(0, no.cells - 4))*10
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
# p3 = ggplot(result, aes(x = factor(x/4), y = CLV3)) + geom_violin()
aWUS.aes = aes(x=factor(x/4), y = aWUS)
aWUS.aes2 = aes(x=as.integer(x/4) + 1, y = aWUS)

pWUS.aes = aes(x=factor(x/4), y = pWUS)
CLV3.aes = aes(x=factor(x/4), y = CLV3)
results  = result %>% filter(cell < 7)
# ggplot(results) +

d2ap = expression(paste("Distance to apex [", mu, "m]"))
p1 = ggplot(result %>% filter(cell < 7)) +
  geom_line(aWUS.aes, size = 1, alpha = .7) +
  geom_line(aWUS.aes2, size = 1, alpha = .7) +
  geom_point(aWUS.aes, size = 3, alpha = .7) +
  labs(x = "", y = "WUS") + theme_bw()
p2 = ggplot(result %>% filter(cell < 7)) +
  geom_jitter(pWUS.aes, alpha = .05, colour = 1) +
  geom_violin(pWUS.aes, bw = .7) + labs(x = "", y = "wus") +
  theme_bw() + labs(x = "")
p3 = ggplot(result %>% filter(cell < 7)) +
  geom_jitter(CLV3.aes, alpha = .05, colour = 2) +
  geom_violin(CLV3.aes, bw = 3) +
  labs(x = "Distance to apex [#neighbours]", y = "CLV3") + theme_bw()

grid.arrange(p1, p2, p3)

