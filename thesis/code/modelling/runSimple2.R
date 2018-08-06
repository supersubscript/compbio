#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/modelling/")
library(tidyverse)
source("organism_wrapper.R")

simulator.path    = "~/Organism/bin/simulator"
model.path        = "models/simpleRadial2.model"
init.path         = "inits/simpleRadial2.init"
init.ss.path      = "inits/simpleRadial2_ss.init"
det.solver.path   = "solvers/rk4.solver"
# stoch.solver.path = "solvers/ms.solver"
stoch.solver.path = "solvers/gillespie.solver"

t.end         = 100
no.prints     = 100
step.size     = 1e-2
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
bin.targets[,1] = bin.targets[,1] / max(bin.targets[,1])
bin.targets[,2] = bin.targets[,2] / max(bin.targets[,2])
sum.targets = unname(colSums(bin.targets))

###############################################################################
### Define reaction parameters
###############################################################################
parameters = c(
  wp = 0.25,
  wD = 1,
  wd = 0.025,
  Cp = 0.6,
  Cv = 0.9,
  Ck = 0.5,
  # Cn = 2.0,
  Cd = 0.1,
  cp = 0.5,
  cD = 2,
  cd = 0.2
)
###############################################################################
### Setup model file
###############################################################################
source("models.R")
write.it = function(parameters) {
  if(is.null(parameters) | length(parameters) == 1) {
    print("ERROR")
    stop()
  }
  wp <<- parameters[1] 
  wD <<- parameters[2]  
  wd <<- parameters[3]  
  Cp <<- parameters[4]  
  Cv <<- parameters[5] 
  Ck <<- parameters[6] 
  Cn <<- 2
  Cd <<- parameters[7] 
  cp <<- parameters[8] 
  cD <<- parameters[9] 
  cd <<- parameters[10] 
  parameters <<-  c(wp, wD, wd, Cp, Cv, Ck, Cn=Cn, Cd, cp, cD, cd)
  
  model.obj             <<- init.model("simpleRadial2")
  addTopo(model.obj)    <<- topo.3d()
  addSpecies(model.obj) <<- init.species("Sink")
  addSpecies(model.obj) <<- init.species("aWUS")
  pwus                  <<- init.species("pWUS")
  addReaction(pwus)     <<- creationOne(index = model.obj@species$aWUS@index, rate = wp)
  addReaction(pwus)     <<- diffusionSimple(rate = wD)
  addReaction(pwus)     <<- degradationOne(rate = wd)
  addReaction(pwus)     <<- degradationTwo(index = model.obj@species$Sink@index, rate = wD)
  addSpecies(model.obj) <<- pwus
  clv3                  <<- init.species("CLV3")
  addReaction(clv3)     <<- creationOne(index = model.obj@species$pWUS@index, rate = Cp)
  addReaction(clv3)     <<- hill(act.ind = model.obj@species$pWUS@index + 2L, rep.ind = list(), vars = c(Cv, Ck, Cn))
  addReaction(clv3)     <<- degradationOne(rate = Cd)
  addSpecies(model.obj) <<- clv3
  pclv3                 <<- init.species("pCLV3")
  addReaction(pclv3)    <<- creationOne(index = model.obj@species$CLV3@index, rate = cp)
  addReaction(pclv3)    <<- diffusionSimple(rate = cD)
  addReaction(pclv3)    <<- degradationOne(rate = cd)
  addReaction(pwus)     <<- degradationTwo(index = model.obj@species$Sink@index, rate = cD)
  addSpecies(model.obj) <<- pclv3
  addNeighRule(model.obj) <<- neighborhoodDistanceSphere(index = 3) 
  write.model(model.obj, file.path = model.path)
}
write.it(parameters)

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
init[, "aWUS"]  = c(10, 8, 6, 4, 2, rep(0, no.cells - 5)) * 5 # remove modifier
init[, "pWUS"]  = rep(0, no.cells)
init[, "CLV3"]  = rep(0, no.cells)
init[, "pCLV3"] = rep(0, no.cells)
init[, "Sink"]  = c(rep(0, no.cells - 1), c(1))

# Write to file
write.init(outfile = init.path, init)

###############################################################################
### Simulate system
###############################################################################
# parameters = parameters[-c("Cn")]
parameters = parameters[-7]
cost.fct = function(parameters) {
  wp = parameters[1]
  wD = parameters[2]
  wd = parameters[3]
  Cp = parameters[4]
  Cv = parameters[5]
  Ck = parameters[6]
  Cn = 2
  Cd = parameters[7]
  cp = parameters[8]
  cD = parameters[9]
  cd = parameters[10]
  parameters = c(wp, wD, wd, Cp, Cv, Ck, Cn, Cd, cp, cD, cd)
  
  write.it(parameters)
  print(parameters)
  # print(parameters)
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
  output = result %>% group_by(x) %>% summarize(mean = mean(CLV3), 
                                                sd   = sd(CLV3) / mean(CLV3)) %>% .[, -1]
  output[,1] = output[,1] / max(output[,1])
  res = colSums(abs(output - bin.targets)) / sum.targets
  res[is.nan(res)] = 999999
  res = unname(res[1] + 3*res[2])
  cat(res, "\n")
  -res
}

# library(cmaes)
# library(GenSA)
# outp = cma_es(
#   parameters,
#   cost.fct,
#   lower = 10e-2,
#   upper = 10e3,
#   control = list(stopfitness = 1.0, keep.best = TRUE)
# )
# outp = GenSA(
#   parameters,
#   fn=cost.fct,
#   lower = rep(10e-3, length(parameters)),
#   upper = rep(10e0, length(parameters)),
#   control = list(verbose=TRUE,maxit=3,max.call=3)
# )

cl = makeCluster(3)
clusterExport(cl, ls())
registerDoParallel(cl)

GA <- ga(type = "real-valued", 
         fitness =  cost.fct,
         min = rep(10e-3, length(parameters)), max = rep(10e0, length(parameters)), 
         popSize = 10, maxiter = 100)


# print(outp$par)
print(GA@solution)
outp = list()
outp$par = GA@solution
wp = outp$par[1]
wD = outp$par[2]
wd = outp$par[3]
Cp = outp$par[4]
Cv = outp$par[5]
Ck = outp$par[6]
Cn = 2
Cd = outp$par[7]
cp = outp$par[8]
cD = outp$par[9]
cd = outp$par[10]

parameters = c(wp, wD, wd, Cp, Cv, Ck, Cd, cp, cD, cd)

write.it(parameters)
# save(file="state.RData", list=ls())
if(!is.null(outp$par))
  result = simulate(
    simulator.path,
    model.path,
    init.ss.path,
    stoch.solver.path,
    toR         = TRUE,
    output.mode = 2
  )
colnames(result)[1:(model.obj@no.species + 8)] = c("iter", "t", "cell", "nNeigh", "x", "y", "z", "r", names(model.obj@species))

library(gridExtra)
p1 = ggplot(result, aes(x = x, y = aWUS)) + geom_line() + geom_point()
p2 = ggplot(result, aes(x = factor(x), y = pWUS)) + geom_point()
p3 = ggplot(result, aes(x = factor(x), y = CLV3)) + geom_point()

grid.arrange(p1, p2, p3)
