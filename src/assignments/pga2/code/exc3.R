setwd("~/compbio/src/assignments/pga2/code/")
library(RColorBrewer)
library(scales) #imports alpha
library(stats)
palette(brewer.pal(n = 8, name = "Set1"))
line.size = 2

data = read.csv("../data/four_population_sequencing_data.csv")
data.firstind = data[, which(grepl("0", colnames(data)))] # get first individual

get.D = function(test.data, frac = FALSE) {
  if (!frac) {
    n.abba = sum(apply(test.data, 1, function(x) x[1] == 0 & 
                         x[1] != x[2] & 
                         x[2] == x[3] & 
                         x[3] != x[4])) # 0 = ancestral
    n.baba = sum(apply(test.data, 1, function(x) x[1] == 1 & 
                         x[1] != x[2] & 
                         x[2] != x[3] & 
                         x[3] != x[4]))
    D = (n.abba - n.baba) / (n.abba + n.baba)
  } else{
    As = test.data[, which(grepl("A_", colnames(test.data)))]
    Bs = test.data[, which(grepl("B_", colnames(test.data)))]
    Cs = test.data[, which(grepl("C_", colnames(test.data)))]
    Ds = test.data[, which(grepl("D_", colnames(test.data)))]
    As = apply(As, 1, mean); Bs = apply(Bs, 1, mean); Cs = apply(Cs, 1, mean); Ds = Ds
    ddd = cbind(As, Bs, Cs, Ds)
    n1 = sum(apply(ddd, 1, function(x) (1 - x[1]) * x[2] * x[3] * (1 - x[4])))
    n2 = sum(apply(ddd, 1, function(x) (x[1]) * (1 - x[2]) * x[3] * (1 - x[4])))
    D = (n1 - n2) / (n1 + n2)
  }
  list(D=D, n.abba=n.abba, n.baba=n.baba)
}
first.data = get.D(data.firstind)
D = first.data[[1]]
n.abba = first.data[[2]]
n.baba = first.data[[3]]

abba.baba.test = function(data, n.folds = 10, frac = FALSE, Dval) {
  folds.i = sample(rep(1:n.folds, length.out = nrow(data)))
  Ds = vector(mode = "numeric", length = n.folds)
  for (ii in 1:n.folds) {
    test.i = which(folds.i == ii)
    train = data[-test.i,]
    Ds[ii] = get.D(train, frac)[[1]]
  }
  Z = Dval / sqrt(n.folds * var(Ds))
  p = 1 - pnorm(abs(Z))
  p
}

# Binom
p.binom = binom.test(n.abba, n.abba + n.baba, 0.5) 

# Using fractions instead
D.frac = get.D(data, frac = TRUE)[[1]]
p = abba.baba.test(data.firstind, n.folds = nrow(data.firstind), frac = FALSE, Dval = D)
p.frac = abba.baba.test(data, n.folds = nrow(data), frac = TRUE, Dval = D.frac)

# Histogram
ll = c()
for(ii in 2:(nrow(data)-1)){
  ll = c(ll, min(data[ii,1] - data[ii-1, 1], data[ii+1,1] - data[ii,1]))
}
ll = c(ll, data[2,1] - data[1,1])  
ll = c(ll, data[nrow(data),1] - data[nrow(data)-1,1])    
hist(ll, breaks = 7, col = 1, density = 100, xlab = "Minimal distance", main = "")

# Perform jack-blockknife
jack.blockknife = 30000
data[,1] = data[,1] / jack.blockknife
blocks = as.integer(data[,1]) + 1
frac = TRUE
n.folds = length(unique(blocks))#max(blocks)

Ds = vector(mode = "numeric", length = n.folds)
count = 1
for (ii in sort(unique(blocks))) {
  test.i = which(blocks == ii)
  train = data[-test.i, ]
  Ds[count] = get.D(train, frac)[[1]]
  count = count + 1 
}
Z.jack = D.frac / sqrt(n.folds * var(Ds, na.rm = T))
p.jack = 1 - pnorm(abs(Z.jack))
p.jack















