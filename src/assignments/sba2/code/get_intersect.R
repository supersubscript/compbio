setwd("/home/henrik/compbio/src/assignments/sba2/code")
mol1.data.mine = as.character(read.table("../data/3uon_within_5.dat.unique")[, 1])
mol1.data.key = as.character(read.table("../data/3uon_key_binding_sites.txt", header = F)[, 1])
mol1.data.key = mol1.data.key[-c((length(mol1.data.key) - 1):length(mol1.data.key))]

mol2.data.mine = as.character(read.table("../data/4mqs_within_5.dat.unique")[, 1])
mol2.data.key = as.character(read.table("../data/4mqs_key_binding_sites.txt", header = F)[, 1])

mol3.2cu.data.mine = as.character(read.table("../data/4mqt_within_5_2CU.dat.unique")[, 1])
mol3.2cu.data.key = as.character(read.table("../data/4mqt_2CU_key_binding_sites.txt", header = F)[, 1])

mol3.ixo.data.mine = as.character(read.table("../data/4mqt_within_5_IXO.dat.unique")[, 1])
mol3.ixo.data.key = as.character(read.table("../data/4mqt_IXO_key_binding_sites.txt", header = F)[, 1])
mol1.paper.important = c("ASP103",
                         "ASN404",
                         "PHE181")

mol2.paper.important = c("ASP103",
                         "ASN404")

mol3.cu.paper.important = c("ASN410",
                            "TRP422",
                            "TYR426",
                            "TYR80",
                            "ASN419",
                            "TYR177",
                            "GLU172")

mol3.ixo.paper.important = c("ASP103",
                             "ASN404")



all.cases = list(
  mol1.data.key,
  mol1.data.mine,
  mol1.paper.important,
  mol2.data.key,
  mol2.data.mine,
  mol2.paper.important,
  mol3.2cu.data.key,
  mol3.2cu.data.mine,
  mol3.cu.paper.important,
  mol3.ixo.data.key,
  mol3.ixo.data.mine,
  mol3.ixo.paper.important
)

all.aa = sort(unique(
  c(
    mol1.data.key,
    mol1.data.mine,
    mol2.data.key,
    mol2.data.mine,
    mol3.2cu.data.key,
    mol3.2cu.data.mine,
    mol3.ixo.data.key,
    mol3.ixo.data.mine,
    mol1.paper.important,
    mol3.ixo.paper.important,
    mol3.cu.paper.important,
    mol2.paper.important
  )
))

out.data = data.frame(aa = all.aa)
for (case in all.cases) {
  out = vector(mode = "double", length = nrow(out.data))
  for (ii in 1:nrow(out.data))
    if (out.data[ii, 1] %in% case)
      out[ii] = 1
    out.data = cbind(out.data, out)
}
rownames(out.data) = out.data[, 1]
out.data = out.data[, -1]
colnames(out.data) = c(
  "3UON PDB",
  "3UON OWN",
  "3UON PAPER",
  "4MQS PDB",
  "4MQS OWN",
  "4MQS PAPER",
  "4MQT.2CU PDB",
  "4MQT.2CU OWN",
  "4MQT.2CU PAPER",
  "4MQT.IXO PDB",
  "4MQT.IXO OWN",
  "4MQT.IXO PAPER"
)

library(RColorBrewer)
library(scales) #imports alpha
library(stats)
palette(brewer.pal(n = 8, name = "Set1"))
line.size = 2

par(mar=c(8,10,2,4))
x = data.matrix(out.data)
image(x[, ncol(x):1], xaxt = "n", yaxt = "n", col = 0:5)
axis(2, at = seq(0, 1, length.out = ncol(x)), labels = rev(colnames(x)), las = 2)
axis(1, at = seq(0, 1, length.out = nrow(x)), labels = rownames(x), las = 2)
abline(h = mean(seq(0, 1, length.out = ncol(x))[9:10]), lwd = 2)
abline(h = mean(seq(0, 1, length.out = ncol(x))[6:7]), lwd = 2)
abline(h = mean(seq(0, 1, length.out = ncol(x))[3:4]), lwd = 2)
abline(h = mean(seq(0, 1, length.out = ncol(x))[2:3]), lwd = 1, lty = 2)
abline(h = mean(seq(0, 1, length.out = ncol(x))[1:2]), lwd = 1, lty = 2)
abline(h = mean(seq(0, 1, length.out = ncol(x))[4:5]), lwd = 1, lty = 2)
abline(h = mean(seq(0, 1, length.out = ncol(x))[5:6]), lwd = 1, lty = 2)
abline(h = mean(seq(0, 1, length.out = ncol(x))[7:8]), lwd = 1, lty = 2)
abline(h = mean(seq(0, 1, length.out = ncol(x))[8:9]), lwd = 1, lty = 2)
abline(h = mean(seq(0, 1, length.out = ncol(x))[10:11]), lwd = 1, lty = 2)
abline(h = mean(seq(0, 1, length.out = ncol(x))[11:12]), lwd = 1, lty = 2)



