setwd("/home/henrik/compbio/thesis/code/modelling/")
source("models.R")

model = init.model("simpleRadial1.model")

### Add topology
addTopo(model) = topo.3d()

### Add species
addSpecies(model) = init.species("Sink")
addSpecies(model) = init.species("aWUS")

wus = init.species("WUS")
addReaction(wus)  = creationOne(index = model@species$aWUS@index)
addReaction(wus)  = degradationOne()
addSpecies(model) = wus

pwus = init.species("pWUS")
addReaction(pwus) = creationOne(index = model@species$WUS@index)
addReaction(pwus) = diffusionSimple()
addReaction(pwus) = degradationOne()
addReaction(pwus) = degradationTwo(index = model@species$Sink@index)
addSpecies(model) = pwus

clv3 = init.species("CLV3")
addReaction(clv3) = creationOne(index = model@species$pWUS@index)
addReaction(clv3) = degradationOne()
addSpecies(model) = clv3

### Add reactions
# No reaction

### Add neighborhoodRule
addNeighRule(model) = neighborhoodDistanceSphere(index = 3) 

# Write model to file
write.model(model, file.path = "models/simpleRadial_autotest.model")
