#!/usr/bin/env Rscript
require(methods)
HOMEDIR = "/home/henrik/spa3/"
TSPDIR  = "tsp_data/"
setwd(HOMEDIR)
OPTFILE = ifelse(length(args) == 0, "tsp_data/a280.opt.tour", args[1])
FILE = gsub("\\.opt\\.tour", "\\.tsp", OPTFILE)
trim = function (x) gsub("^\\s+|\\s+$", "", x) # Trim whitespaces
# source("aux.R")

read.coord = function(f){
  lines = readLines(f)
  lines = trim(lines)
  lines = lines[which(lines != "EOF")]
  lines = lines[which(sapply(lines, length) > 0)]
  first = grep("^1\\s+|^0001\\s+", lines)[1]
  lines = lines[-(1:(first-1))]
  lines = sapply(lines, function(x) strsplit(x, "\\s+")[[1]])
  lines = unlist(lines)
  lines = lines[which(lines != "")]
  data  = matrix(lines, byrow=TRUE, ncol=3)
  # print(data)
  data  = data.frame(data)
  data
}

read.opt = function(f){
  #  print(f)
  lines = readLines(f)
  lines = trim(lines)
  lines = lines[which(lines != "EOF")]
  lines = lines[which(sapply(lines, length) > 0)]
  first = grep("^1", lines)[1]
  lines = lines[-(1:(first-1))]
  lines = sapply(lines, function(x) strsplit(x, "\\s+")[[1]])
  lines = unlist(lines)
  lines = lines[which(lines != "")]
  lines = as.numeric(lines)
  lines = lines[which(lines > 0)]
  lines
}


# Read in data
file          = paste0(HOMEDIR, TSPDIR, FILE)
cities        = read.coord(FILE)
cities        = as.matrix(cities)
no.cities     = nrow(cities)
class(cities) = "numeric"   

# Define organism class
org = function(tour, fit) list(tour=tour, fit=fit)

# Define error measure and fitness
tour.length = function(ind){
  sum(
    sqrt((cities[ind[1:(no.cities-1)], 2] - cities[ind[2:(no.cities)], 2])**2 +
           (cities[ind[1:(no.cities-1)], 3] - cities[ind[2:(no.cities)], 3])**2)) +
    sqrt((cities[ind[no.cities], 2]       - cities[ind[1], 2])**2 +
           (cities[ind[no.cities], 3]       - cities[ind[1], 3])**2)  
}
fitness = function(ind){tour.length(ind)}
random  = function(){tour = sample(1:no.cities, no.cities); return(org(tour=tour, fit=fitness(tour)))}

# Parameters
generations = 5e6
temp        = 1e3
temp.rate   = 1 - 1e-4
state       = random()

# Same as our good old mutate.reverse
state.neighbour = function(state){
  newtour = state$tour
  randseq = sample(no.cities, 2)
  randseq = randseq[1]:randseq[2]
  newtour[randseq] = rev(newtour[randseq])
  return(org(tour=newtour,fit=fitness(newtour)))
}

start.time =  proc.time()
for(iter in 1:generations){
  # Suggest update
  proposed.state = state.neighbour(state)

  # Is it better?
  if(proposed.state$fit <= state$fit)  
    state = proposed.state
  # Should we still go there?
  else if (runif(1) < exp((state$fit - proposed.state$fit)/temp)){
    state <<- proposed.state
  }

  # Change temperature
  temp = temp * temp.rate
  
  # Print current result
  if (iter %% 10000 == 0)
  {
    cat(iter, "\t", temp, "\t", state$fit, "\n")
    solution = cities[state$tour, ]
    plot(solution[, 2], solution[, 3], type='b')
  }
}

# Print bezt
best = population[which(sapply(population, function(x) x$fit) == best.fit)][[1]]
write(c(gsub("\\.tsp", "", FILE), best$fit, tour.length(read.opt(OPTFILE))), ncolumns=3, sep = "\t", file = paste0(gsub("\\.tsp", "", FILE), "_simann_solution.txt"))
