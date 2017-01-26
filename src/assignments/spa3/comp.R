#!/usr/bin/env Rscript
setwd("~/spa3")
library('pixmap')
library('webp')
library('grid')
library('gridExtra')

# Take in data
lena.file     = "lena512color.tiff"
fabio.file    = "fab_fabio.webp"
lena.img.name = strsplit(lena.file, "\\.")[[1]]
lena.img.name = ifelse(length(lena.img.name) > 1, lena.img.name[-length(lena.img.name)], lena.img.name)
lena.img.name = paste0(lena.img.name, ".ppm")
aesthet       = system(paste0("convert ", lena.file, " ", lena.img.name))
lena.image    = read.pnm(lena.img.name, cellres=1)
fab.file      = read_webp("fab_fabio.webp")

# Retrieve matrix specific data
lena.red.matrix    = matrix(lena.image@red,   nrow = ncol(lena.image@red), 
                            ncol = nrow(lena.image@red))
lena.green.matrix  = matrix(lena.image@green, nrow = ncol(lena.image@red), 
                            ncol = nrow(lena.image@red))
lena.blue.matrix   = matrix(lena.image@blue,  nrow = ncol(lena.image@red), 
                            ncol = nrow(lena.image@red))
fabio.red.matrix   = matrix(fab.file[,,1], nrow = ncol(fab.file[,,1]), 
                            ncol = nrow(fab.file[,,1]))
fabio.green.matrix = matrix(fab.file[,,2], nrow = ncol(fab.file[,,2]), 
                            ncol = nrow(fab.file[,,1]))
fabio.blue.matrix  = matrix(fab.file[,,3], nrow = ncol(fab.file[,,3]), 
                            ncol = nrow(fab.file[,,1]))

# Define parameters
iterations  = 10000
grid.width  = 10
grid.height = 10
rad.start   = 2
rad.deg     = 0.9
block.side  = 4
no.nodes    = grid.width*grid.height
nodes       = expand.grid(1:grid.width, 1:grid.height)
nodes       = nodes[, c(2,1)] # Loop-intuitive order

# Define node distances
nodes.dist = function(x) as.matrix(stats::dist(x, method = "euclidian"))
distances  = nodes.dist(nodes)

# Retrieve Kohonen interpreted matrix
get.matrix = function(mat, method = "uniform"){
  cat("\nInitialising new Kohonen compression\n")
  # Initialise
  alpha   = seq(0.2, 0.01, length.out = iterations)
  radii   = rep(rad.start, grid.width * grid.height)
  weights = matrix(runif(grid.height*grid.width*block.side**2, max = .1), 
                   grid.width * grid.height, block.side**2)

  # Have a guess
  get.random.block = function(){
    ii = sample(ncol(mat) / block.size, 1)
    jj = sample(nrow(mat) / block.size, 1)
    x.ind = ((ii-1)*block.side + 1):((ii)*block.side)
    y.ind = ((jj-1)*block.side + 1):((jj)*block.side)
    as.vector(mat[x.ind, y.ind])
  }

  # Create all the input arrays
  data = list()
  for(ii in 1:(ncol(mat)/block.side)) {
    for(jj in 1:(nrow(mat)/block.side)) {
      x.ind = ((ii-1)*block.side + 1):((ii)*block.side)
      y.ind = ((jj-1)*block.side + 1):((jj)*block.side)
      if(method == "uniform")
        data = append(data, list(as.vector(mat[x.ind, y.ind]))) 
      else if (method == "random")
        data = append(data, list(get.random.block()))
    }
  }

  # GO!
  for (iter in 1:iterations){
    # Get input-weights distance
    dists     = apply(weights, 1, function(x) sum((data[[((iter - 1) %% (length(data)))+1]] - x)**2))
    minimum   = which(dists == min(dists))[1]

    # Find the neighbours and update them accordingly
    neighbours = which(distances[minimum,] < radii[minimum])
    datum      = matrix(rep(data[[((iter - 1) %% (length(data)))+1]], 
                            length(neighbours)), ncol=ncol(weights), byrow=TRUE)
    weights[neighbours, ] = weights[neighbours, ] + alpha[iter] * (datum - weights[neighbours, ])

    # Update radius for the winning node
    radii[minimum] = rad.deg * radii[minimum]

    # Print progress
    if(iter %% 1000 == 0)
      print(iter)
  }

  # Convert back to complete matrix
  # Note: In practice we store the (minimum) indices and the net.
  new.image = matrix(0, ncol(mat), nrow(mat))
  sapply(1:(ncol(mat) / block.side),
    function(ii) sapply(1:(nrow(mat) / block.side),
      function(jj) {
        x.ind = ((ii - 1) * block.side + 1):((ii) * block.side)
        y.ind = ((jj - 1) * block.side + 1):((jj) * block.side)
        input = as.vector(mat[x.ind, y.ind])
        dists = apply(weights, 1, function(x) sum((input - x) ** 2))
        minimum = which(dists == min(dists))[1]
        new.image[x.ind, y.ind] <<- weights[minimum, ]
      }))
  return(new.image)
}

# Take out the images
lena.image.red    = get.matrix(lena.red.matrix,    method = "uniform")
lena.image.blue   = get.matrix(lena.blue.matrix,   method = "uniform")
lena.image.green  = get.matrix(lena.green.matrix,  method = "uniform")
fabio.image.red   = get.matrix(fabio.red.matrix,   method = "uniform")
fabio.image.blue  = get.matrix(fabio.blue.matrix,  method = "uniform")
fabio.image.green = get.matrix(fabio.green.matrix, method = "uniform")

# Plot images
par(mfrow=c(2,3), mar = c(1,1,1,1))
image(t(apply(lena.red.matrix,2,rev)),xaxt="n",yaxt="n",main="Red")  
image(t(apply(lena.blue.matrix,2,rev)),xaxt="n",yaxt="n",main="Green")  
image(t(apply(lena.green.matrix,2,rev)),xaxt="n",yaxt="n",main="Blue")
image(t(apply(fabio.red.matrix,2,rev)),xaxt="n",yaxt="n")
image(t(apply(fabio.blue.matrix,2,rev)),xaxt="n",yaxt="n")
image(t(apply(fabio.green.matrix,2,rev)),xaxt="n",yaxt="n")
image(t(apply(lena.image.red,2,rev)),xaxt="n",yaxt="n",main="Red")
image(t(apply(lena.image.blue,2,rev)),xaxt="n",yaxt="n",main="Green")
image(t(apply(lena.image.green,2,rev)),xaxt="n",yaxt="n",main="Blue")
image(t(apply(fabio.image.red,2,rev)),xaxt="n",yaxt="n")
image(t(apply(fabio.image.blue,2,rev)),xaxt="n",yaxt="n")
image(t(apply(fabio.image.green,2,rev)),xaxt="n",yaxt="n")

# Merge channels
lena.all.three       = rgb(lena.image.red,  lena.image.green,  lena.image.blue)
fabio.all.three      = rgb(fabio.image.red, fabio.image.green, fabio.image.blue)
dim(lena.all.three)  = dim(lena.image.red)
dim(fabio.all.three) = dim(fabio.image.red)

# Prepare to plot
grob.size  = .95
ll         = rgb(lena.image@red, lena.image@green, lena.image@blue)
rr         = rgb(fab.file[,,1], fab.file[,,2], fab.file[,,3])
dim(ll)    = dim(lena.image.red)
dim(rr)    = dim(fabio.image.red)
l.grob     = rasterGrob(ll, interpolate = FALSE, width = grob.size)
f.grob     = rasterGrob(rr, interpolate = FALSE, width = grob.size)
lena.grob  = rasterGrob(lena.all.three,  interpolate=FALSE, width=grob.size)
fabio.grob = rasterGrob(fabio.all.three, interpolate=FALSE, width=grob.size)

#################################
# Compare to SVD
#################################
do.svd = function(mat, value){
  mat.svd = svd(mat)
  d = mat.svd$d; u = mat.svd$u; v = mat.svd$v
  mat.reconstruction = u #latex# diag(d) #latex# t(v)
  mat.compressed = u[,1:value] #latex# diag(d[1:value]) #latex# t(v[,1:value])
}

# Do SVD stuff, account for overshooting
lena.svd.red    = do.svd(lena.red.matrix,    30)
lena.svd.green  = do.svd(lena.green.matrix,  30)
lena.svd.blue   = do.svd(lena.blue.matrix,   30)
fabio.svd.red   = do.svd(fabio.red.matrix,   30)
fabio.svd.green = do.svd(fabio.green.matrix, 30)
fabio.svd.blue  = do.svd(fabio.blue.matrix,  30)
lena.svd.red[which(lena.svd.red       < 0)] = 0
lena.svd.green[which(lena.svd.green   < 0)] = 0
lena.svd.blue[which(lena.svd.blue     < 0)] = 0
lena.svd.red[which(lena.svd.red       > 1)] = 1
lena.svd.green[which(lena.svd.green   > 1)] = 1
lena.svd.blue[which(lena.svd.blue     > 1)] = 1
fabio.svd.red[which(fabio.svd.red     < 0)] = 0
fabio.svd.green[which(fabio.svd.green < 0)] = 0
fabio.svd.blue[which(fabio.svd.blue   < 0)] = 0
fabio.svd.red[which(fabio.svd.red     > 1)] = 1
fabio.svd.green[which(fabio.svd.green > 1)] = 1
fabio.svd.blue[which(fabio.svd.blue   > 1)] = 1
lena.all.three.svd       = rgb(lena.svd.red,  lena.svd.green,  lena.svd.blue)
fabio.all.three.svd      = rgb(fabio.svd.red, fabio.svd.green, fabio.svd.blue)
dim(lena.all.three.svd)  = dim(lena.svd.red)
dim(fabio.all.three.svd) = dim(fabio.svd.red)

# Plot stuff
par(mfrow=c(1,2)); grid.newpage();
lena.svd.grob  = rasterGrob(lena.all.three.svd,  interpolate=FALSE, width=grob.size)
fabio.svd.grob = rasterGrob(fabio.all.three.svd, interpolate=FALSE, width=grob.size)
grid.arrange(grobs=list(l.grob, lena.grob, lena.svd.grob, f.grob,  
                        fabio.grob,  fabio.svd.grob), ncol = 3, nrow=2)
