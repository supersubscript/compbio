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
    n.abba = sum(apply(test.data, 1, function(x)
      x[1] == 0 &
        x[1] != x[2] & x[2] == x[3] & x[3] != x[4])) # 0 = ancestral
    n.baba = sum(apply(test.data, 1, function(x)
      x[1] == 1 & x[1] != x[2] & x[2] != x[3] & x[3] != x[4]))
    D = (n.abba - n.baba) / (n.abba + n.baba)
  } else{
    As = test.data[, which(grepl("A_", colnames(test.data)))]
    Bs = test.data[, which(grepl("B_", colnames(test.data)))]
    Cs = test.data[, which(grepl("C_", colnames(test.data)))]
    Ds = test.data[, which(grepl("D_", colnames(test.data)))]
    As = apply(As, 1, function(x)
      sum(x) / length(x))
    Bs = apply(Bs, 1, function(x)
      sum(x) / length(x))
    Cs = apply(Cs, 1, function(x)
      sum(x) / length(x))
    Ds = Ds
    ddd = cbind(As, Bs, Cs, Ds)
    n1 = sum(apply(ddd, 1, function(x)
      (1 - x[1]) * x[2] * x[3] * (1 - x[4])))
    n2 = sum(apply(ddd, 1, function(x)
      (x[1]) * (1 - x[2]) * x[3] * (1 - x[4])))
    D = (n1 - n2) / (n1 + n2)
  }
  D
}
D = get.D(data.firstind)

cat("n.abba = ", n.abba, "\n")
cat("n.baba = ", n.baba, "\n")
cat("D = ", D, "\n")


# See some introgression between populations A and C 
# Wouldn't expect introgression into either in preference since they are evolutionarily on the same level hierarchically?
##############################
### Assuming independent
##############################
# Ask others? 
p.binom = binom.test(n.abba + n.baba, nrow(data.firstind), 0.5) 
cat("p, binom = ", p.binom$p.value, "\n")

##############################
### Assuming non-independent
##############################
abba.baba.test = function(data, n.folds = 10, frac = FALSE) {
  folds.i = sample(rep(1:n.folds, length.out = nrow(data)))
  Ds = vector(mode = "numeric", length = n.folds)
  for (ii in 1:n.folds) {
    test.i = which(folds.i == ii)
    train = data[-test.i,]
    Ds[ii] = get.D(train, frac)
  }
  Z = D / sqrt(n.folds * var(Ds))
  p = 1 - pnorm(Z)
  p
}
p = abba.baba.test(data.firstind, nrow(data.firstind))
cat("p = ", p, "\n")

##############################
### D
##############################
D.frac = get.D(data, frac = TRUE)
cat("D, frac = ", D.frac, "\n")
p.frac = abba.baba.test(data, nrow(data), frac = TRUE)
cat("p, frac = ", p.frac, "\n")

##############################
### G
##############################
















