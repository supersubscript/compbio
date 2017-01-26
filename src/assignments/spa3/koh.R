require("stats")
grid.width  = 5
grid.height = grid.width
initial.rad = 3
radii       = rep(initial.rad, grid.width * grid.height)
iterations  = 10000
data        = matrix(runif(iterations*2),  iterations, 2)
input.dim   = ncol(data)
rad.deg     = 0.9999
alpha       = seq(0.05, 0.01, length.out = iterations)
no.images   = 3

# Init
# weights     = expand.grid(1:grid.width, 1:grid.height) / grid.width
weights   = matrix(runif(grid.width**2*input.dim,min=.4,max=.6),grid.width*grid.height, input.dim)
nodes     = expand.grid(1:grid.width, 1:grid.height)
distances = as.matrix(stats::dist(nodes, method = "maximum")) # Distances between nodes (grid-space)

# Plot grid
plot.grid = function(main="", sub=""){
  plot(weights,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="", xaxt = "n",yaxt="n", bty="n", pch = 16)
  for (ii in 1:grid.width) 
    for (jj in 1:grid.height) {
      if (ii < grid.width)
        lines(rbind(weights[which(nodes[, 1] == ii & nodes[, 2] == jj), ], weights[which(nodes[, 1] == (ii + 1) & nodes[, 2] == jj),]))
      if (jj < grid.height) 
        lines(rbind(weights[which(nodes[, 1] == ii & nodes[, 2] == jj), ], weights[which(nodes[, 1] == ii & nodes[, 2] == (jj + 1)),]))
    }
  title(main=main, line=-1)
}

# GO!
plot.grid(main=paste0("Iteration: ",0))
for (iter in 1:iterations){
  if(iter %% trunc(iterations / no.images) == 0) plot.grid(main=paste0("Iteration: ",iter))
  
  # Get input-weights distance
  dists     = apply(weights, 1, function(x) sum((data[iter,] - x)**2))
  minimum   = which(dists == min(dists))[1]
  
  # Find the neighbours and update them accordingly
  neighbours = which(distances[minimum,] < radii[minimum])
  datum      = matrix(rep(data[iter, ],length(neighbours)),ncol=ncol(weights),byrow=TRUE)#rep(train.data[iter, ], length(neighbours)), ncol=ncol(weights), byrow=TRUE)
  weights[neighbours, ] = weights[neighbours, ] + alpha[iter] * (datum - weights[neighbours, ])
  
  # Update radius for the winning node
  radii[minimum] = rad.deg * radii[minimum]
  
  # Print progress  
  if(iter %% 1000 == 0)
    print(iter)
}