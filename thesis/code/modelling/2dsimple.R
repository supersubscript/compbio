#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/modelling/")
library(tidyverse)
library(gridExtra)
source("organism_wrapper.R")

simulator.path    = "~/Organism/bin/simulator"
model.path        = "models/new_rad1.model"
init.path         = "inits/new_rad1.init"
init.ss.path      = "inits/new_rad1_ss.init"
det.solver.path   = "solvers/rk4.solver"
stoch.solver.path = "solvers/ms.solver"
# stoch.solver.path = "solvers/heun.solver"
# stoch.solver.path = "solvers/rk4.solver"
# stoch.solver.path = "solvers/gillespie.solver"
model = "new_rad1"

t.end         = 100
no.prints     = 100
step.size     = 1e-3
output.format = 2
error         = 1e-2
noise         = 1


###############################################################################
### Define reaction parameters
###############################################################################
wp = 5
wD = 2
wd = 1
ws = wD

Cv   = 150
Ck   = 6
Cn   = 4
Cl1k = 2
Cl1n = 4
Cd   = 1


###############################################################################
### Setup init file
###############################################################################
init.solver.gillespie("solvers/gillespie.solver", t.start = 0, t.end = t.end, 
                      no.prints = no.prints, step.size = step.size, 
                      output.format = output.format)
init.solver.milsteinStrat("solvers/ms.solver", t.start = 0, t.end = t.end, 
                          no.prints = no.prints, step.size = step.size, 
                          output.format = output.format, noise = noise, vol.flag = 0)
init.solver.HeunItoDiffusive("solvers/heun.solver", t.start = 0, t.end = t.end, 
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
addSpecies(model.obj) = init.species("L1")

pwus                  = init.species("pWUS")
addReaction(pwus)     = creationOne(index = model.obj@species$aWUS@index, rate = wp)
addReaction(pwus)     = diffusionSimple(rate = wD)
addReaction(pwus)     = degradationOne(rate = wd)
addReaction(pwus)     = degradationTwo(index = model.obj@species$Sink@index, rate = ws)
addSpecies(model.obj) = pwus

clv3                  = init.species("CLV3")
addReaction(clv3)     = hill(act.ind = list(model.obj@species$pWUS@index, model.obj@species$L1@index), 
                             rep.ind = list(), 
                             vars = c(Cv, Ck, Cn, Cl1k, Cl1n))
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
# init.molecules = matrix(runif(no.species * no.cells*2), no.cells*2, no.species)
# init           = cbind(c(seq(from = 0, length.out = no.cells, by = 4), 
#                          seq(from = 0, length.out = no.cells, by = 4)),
#                        c(rep(0, no.cells), rep(4, no.cells)), rep(0, no.cells*2), rep(radius, no.cells*2), 
#                        init.molecules)
# colnames(init) = c("x", "y", "z", "r", names(model.obj@species))

init.molecules = matrix(runif(no.species * no.cells), no.cells, no.species)
init           = cbind(seq(from = 0, length.out = no.cells, by = 4), 
                       rep(0, no.cells), rep(0, no.cells), rep(radius, no.cells), 
                       init.molecules)
colnames(init) = c("x", "y", "z", "r", names(model.obj@species))

ww = 2**(1:5) #seq(1, 101, 10)
ll = 2**(1:5) #seq(1, 101, 10)
# ww = seq(1, 10, 1)
# ll = seq(1, 10, 1)

for (ii in 1:length(ww)) {
  print(ii)
  for (jj in 1:length(ll)) {
    init.molecules = matrix(runif(no.species * no.cells), no.cells, no.species)
    init           = cbind(seq(from = 0, length.out = no.cells, by = 4), 
                           rep(0, no.cells), rep(0, no.cells), rep(radius, no.cells), 
                           init.molecules)
    colnames(init) = c("x", "y", "z", "r", names(model.obj@species))
    

# Set initial values of sink and anchor
init[, "aWUS"] = c(c(1, .0, .0, .0, rep(0, no.cells - 4)) * ww[ii])
init[, "pWUS"] = rep(0, no.cells)
init[, "CLV3"] = rep(0, no.cells)
init[, "L1"]   = c(rep(1, no.cells)) * ll[jj]
init[, "Sink"] = c(c(rep(0, no.cells - 1), c(1)))
# init[, "aWUS"] = c(c(1, 0, 0, 0, rep(0, no.cells - 4))*10, c(2, 1, 0, 0, rep(0, no.cells - 4))*10)
# init[, "pWUS"] = rep(0, no.cells*2)
# init[, "CLV3"] = rep(0, no.cells*2)
# init[, "L1"]   = c(rep(1, no.cells), rep(0, no.cells))*100
# init[, "Sink"] = c(c(rep(0, no.cells - 1), c(1)), c(rep(0, no.cells - 1), c(1)))

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

aWUS.aes = aes(x=factor(x/4), y = aWUS)
aWUS.aes2 = aes(x=as.integer(x/4) + 1, y = aWUS)

pWUS.aes = aes(x=factor(x/4), y = pWUS)
CLV3.aes = aes(x=factor(x/4), y = CLV3)
results  = result %>% filter(cell < 7)
# ggplot(results) +
result = result %>% as.tibble()

d2ap = expression(paste("Distance to apex [", mu, "m]"))
# p1 = ggplot(result) +
#   geom_line(aWUS.aes, size = 1, alpha = .7) +
#   geom_line(aWUS.aes2, size = 1, alpha = .7) +
#   geom_point(aWUS.aes, size = 3, alpha = .7) +
#   labs(x = "", y = "WUS") + theme_bw()
# p2 = ggplot(result) +
#   geom_jitter(pWUS.aes, alpha = .05, colour = 1) +
#   geom_violin(pWUS.aes, bw = .7) + labs(x = "", y = "wus") +
#   theme_bw() + labs(x = "")
p3 = ggplot(result %>% filter(y == 0)) +
  geom_jitter(CLV3.aes, alpha = .05, colour = 2, na.rm = TRUE) +
  # geom_violin(CLV3.aes, bw = 3) +
  theme_bw() + ylim(0,150) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
# p4 = ggplot(result %>% filter(y == 4)) +
#   geom_jitter(CLV3.aes, alpha = .05, colour = 2) +
#   geom_violin(CLV3.aes, bw = 3) +
#   labs(x = "Distance to apex [#neighbours]", y = "CLV3") + 
#   theme_bw() + ylim(0,150)

# grid.arrange(p1, p2, p3)
# grid.arrange(p1, p2, p3)
assign(paste0("p", ii, ".", jj), p3)

  }
}

# grid.arrange(p1.1, p1.2, p1.3, p1.4, p1.5, p1.6, p1.7, p1.8, p1.9, p1.10,
#              p2.1, p2.2, p2.3, p2.4, p2.5, p2.6, p2.7, p2.8, p2.9, p2.10,
#              p3.1, p3.2, p3.3, p3.4, p3.5, p3.6, p3.7, p3.8, p3.9, p3.10,
#              p4.1, p4.2, p4.3, p4.4, p4.5, p4.6, p4.7, p4.8, p4.9, p4.10,
#              p5.1, p5.2, p5.3, p5.4, p5.5, p5.6, p5.7, p5.8, p5.9, p5.10,
#              p6.1, p6.2, p6.3, p6.4, p6.5, p6.6, p6.7, p6.8, p6.9, p6.10,
#              p7.1, p7.2, p7.3, p7.4, p7.5, p7.6, p7.7, p7.8, p7.9, p7.10,
#              p8.1, p8.2, p8.3, p8.4, p8.5, p8.6, p8.7, p8.8, p8.9, p8.10,
#              p9.1, p9.2, p9.3, p9.4, p9.5, p9.6, p9.7, p9.8, p9.9, p9.10,
#              p10.1, p10.2, p10.3, p10.4, p10.5, p10.6, p10.7, p10.8, p10.9, p10.10, ncol = 10)
             
             
grid.arrange(p1.1, p1.2, p1.3, p1.4, p1.5, 
             p2.1, p2.2, p2.3, p2.4, p2.5, 
           p3.1, p3.2, p3.3, p3.4, p3.5, 
             p4.1, p4.2, p4.3, p4.4, p4.5, 
             p5.1, p5.2, p5.3, p5.4, p5.5, ncol = 5, 
             top  = textGrob(expression(L1 %->% ""), gp = gpar(fontsize = 18, col = 'black'),  hjust = 4),
             left = textGrob(expression("" %<-% WUS), gp = gpar(fontsize = 18, col = 'black'), rot = 90, hjust = -1.1))
# pdf(file = "test.pdf", width = 6, height=6)
# graphics.off()
          