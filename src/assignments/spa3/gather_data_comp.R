require(methods)
args    = commandArgs(trailingOnly=TRUE)
HOMEDIR = "/home/henrik/spa3/"
TSPDIR  = "tsp_data/"
setwd(HOMEDIR)

svd.img   = system("ls compress/*svd*", intern=TRUE)
svd.sizes = system("du compress/*svd*", intern=TRUE)
svd.sizes = sapply(svd.sizes, function(x) strsplit(x, "\t")[[1]][1])
svd.vals  = sapply(svd.img, function(x) strsplit(x, "_")[[1]][3])
svd.vals  = sapply(svd.vals, function(x) gsub("\\.png", "", x))

svd.vals  = as.numeric(svd.vals)
svd.sizes = as.numeric(svd.sizes)

########################

koh.img   = system("ls compress/*koh*", intern=TRUE)
koh.sizes = system("du compress/*koh*", intern=TRUE)
koh.sizes = sapply(koh.sizes, function(x) strsplit(x, "\t")[[1]][1])
koh.vals  = sapply(koh.img,   function(x) strsplit(x, "_")[[1]][3])
koh.vals  = sapply(koh.vals,  function(x) gsub("\\.png", "", x))

koh.vals  = as.numeric(koh.vals)
koh.vals  = koh.vals**2
koh.sizes = as.numeric(koh.sizes)

koh.vals  = append(koh.vals[-1],  koh.vals[1])
koh.sizes = append(koh.sizes[-1], koh.sizes[1])

par(mfrow=c(1,2), mar = c(5,4,4,1))
plot(koh.vals, koh.sizes, pch = 15, col = "blue", cex = 2, ylab = "Filesize [kB]", xlab = "Number of nodes")
plot(svd.vals, svd.sizes, pch = 16, col = "red", cex = 2, ylab = "Filesize [kB]", xlab = "Decomposition value")
