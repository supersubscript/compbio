source("http://www.bioconductor.org/biocLite.R")
require(methods)
args    = commandArgs(trailingOnly=TRUE)
HOMEDIR = "/home/henrik/spa3/"
TSPDIR  = "tsp_data/"
setwd(HOMEDIR)
OPTFILE = ifelse(length(args) == 0, "tsp_data/a280.opt.tour", args[1])
FILE = gsub("\\.opt\\.tour", "\\.tsp", OPTFILE)
trim = function (x) gsub("^\\s+|\\s+$", "", x) # Trim whitespaces
source("aux.R")

# Read in data
file      = paste0(HOMEDIR, TSPDIR, FILE)
cities    = read.coord(FILE)
cities    = as.matrix(cities)
class(cities) = "numeric"   
  
# Define organism class
org = function(tour, fit) list(tour=tour, fit=fit)

# Define error measure and fitness
# Note: We could've just precalculated this
tour.length = function(ind){
  sum(
    sqrt((cities[ind[1:(no.cities-1)], 2] - cities[ind[2:(no.cities)], 2])**2 +
         (cities[ind[1:(no.cities-1)], 3] - cities[ind[2:(no.cities)], 3])**2)) +
    sqrt((cities[ind[no.cities], 2]       - cities[ind[1], 2])**2 +
         (cities[ind[no.cities], 3]       - cities[ind[1], 3])**2)  
}

# fitness = function(ind){1 / tour.length(ind)}
fitness = function(ind){tour.length(ind)}
random  = function(){tour = sample(1:no.cities, no.cities); return(org(tour=tour, fit=fitness(tour)))}

# Setting
generations     = 5000000
population.size = 14
no.cities       = nrow(cities)
population      = replicate(population.size, random(), simplify = FALSE) # Initialise population
mutation.rate   = 1.0 / no.cities
crossover.fract = 0.2
mutation.fract  = 0.9
# selection.prob  = 0.9 # For binomial
best.fit        = min(sapply(population, function(x) x$fit))

####################################
### Selection
####################################
tournament.selection = function(n = population.size)
{
  # Keep the best one
  newpop      = vector(mode="list", length = n)
  best.index  = which(sapply(population, function(x) x$fit) == best.fit)[1]
  newpop[[1]] = population[[best.index]]
  population  = population[-best.index]

  if(n < 2)
    return(newpop)

  # Take random pair, pick best
  for(count in 2:n)
  {
    rands = sample(1:length(population), 2, replace = FALSE) # Pick random pair
    index = ifelse(population[[rands[1]]]$fit < population[[rands[2]]]$fit, rands[1], rands[2]) # Sort
    newpop[[count]] = population[[index]] # Take chosen individual and add to new population
    population      = population[-index]  # Remove from population
  }
  return(newpop)
}
####################################
### Crossover
####################################
pillar.crossover = function(ind1, ind2)
{
  rands   = sample(no.cities, 2)
  rands   = sort(rands)
  newind1 = vector(mode = "integer", length = no.cities)
  newind2 = vector(mode = "integer", length = no.cities)

  # Swap first chunk and add in remnants
  newind1[seq(rands[1], rands[2])] = ind1[seq(rands[1], rands[2])]
  newind2[!1:no.cities %in% seq(rands[1], rands[2])] = ind1[!1:no.cities %in% seq(rands[1], rands[2])]
  newind1[-seq(rands[1], rands[2])] = ind2[!ind2 %in% newind1[seq(rands[1],  rands[2])]] 
  newind2[seq(rands[1],  rands[2])] = ind2[!ind2 %in% newind2[-seq(rands[1], rands[2])]] 
  return(list(org(tour=newind1, fit = fitness(newind1)), org(tour=newind2, fit=fitness(newind2))))
}

####################################
### Mutations
####################################
swap.mutation = function(ind, base){
  rand = sample(1:no.cities, length(base))
  replace(ind, c(base, rand), ind[c(rand, base)])
}

mutate.randomise = function(ind){
  mutate.offers = which(runif(no.cities) < mutation.rate * 2) # Which bases are to be swapped?
  if(length(mutate.offers) < 2) return(org(tour=ind$tour, fit=ind$fit))
  ind.seq = ind$tour
  ind.seq[mutate.offers] = sample(ind.seq[mutate.offers], length(mutate.offers))
  return(org(tour=ind.seq, fit=fitness(ind.seq)))
}

mutate.chunk = function(ind)
{
  # Get random intervals
  interval = sample(no.cities, 2)
  interval = sort(interval)
  width    = interval[2] - interval[1]
  pos      = sample(1:(no.cities-width), 1)
  
  # Get regions for exchange
  first.seq       = seq(pos, pos + width)
  second.seq      = seq(interval[1], interval[2])
  diff            = intersect(first.seq, second.seq)
  first.seq       = first.seq[!first.seq   %in% diff]
  second.seq      = second.seq[!second.seq %in% diff]
  
  # Swap chunks
  newind             = ind$tour
  first.batch        = newind[first.seq]
  second.batch       = newind[second.seq]
  newind[first.seq]  = second.batch
  newind[second.seq] = first.batch
  
  return(org(tour=newind,fit=fitness(newind)))
}

# Randomly select two indices, and reverse the sequence in between
mutate.reverse = function(ind){
  newind  = ind$tour
  randseq = sample(no.cities, 2)
  randseq = randseq[1]:randseq[2]
  newind[randseq] = rev(newind[randseq])
  return(org(tour=newind,fit=fitness(newind)))
}

##########################################
## MAIN METHOD
##########################################
start.time =  proc.time()
for(iter in 1:generations){
#while((proc.time() - start.time)["elapsed"] < 1800){ # Run for 30 min
  # Mutate
  # population = append(population, lapply(population, mutate.randomise))
  # population = append(population, lapply(population, mutate.chunk))
  mutation.pop = tournament.selection(ceiling(mutation.fract * population.size))
  
  # Apply crossover
  cross.pop     = tournament.selection(2 * trunc(crossover.fract * population.size/2.0))
  cross.pop.ind = sample(length(cross.pop), length(cross.pop))
  cross.pop.ind = matrix(cross.pop.ind, 2)

  # Append to population
  population   = append(population, lapply(mutation.pop, mutate.reverse))
  population   = append(population, apply(cross.pop.ind, 2, 
                    function(x) pillar.crossover(cross.pop[[x[1]]]$tour, 
                                                 cross.pop[[x[2]]]$tour))[[1]])

  # Apply selection
  population = tournament.selection(population.size)
  
  # Update / print best fitness
  best.fit = min(sapply(population, function(x) x$fit))
   if (iter %% 1000 == 0){
     cat(iter, "\t", best.fit, "\n")
      best = population[which(sapply(population, function(x) x$fit) == best.fit)][[1]]
     solution = rbind(cities[best$tour, ], cities[best$tour[1], ])
      plot(cities[best$tour, 2], cities[best$tour, 3], type='o', pch = 16)
  }
}

# Print best solution
best = population[which(sapply(population, function(x) x$fit) == best.fit)][[1]]
write(c(gsub("\\.tsp", "", FILE), best$fit, tour.length(read.opt(OPTFILE))), ncolumns=3, sep = "\t", file = paste0(gsub("\\.tsp", "", FILE), "_ga_solution.txt"))
