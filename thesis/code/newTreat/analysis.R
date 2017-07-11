setwd("/home/henrik/compbio/thesis/code/newTreat/")
source("functions.R")
source("read_data.R")
plants = c(1, 2, 4, 13, 15, 18)

results = get.plant.data(18)
dd = results$division.data
sd = results$segm.data
qd = results$quant.data
tp = results$timepoints
last = length(tp)

nCells = sapply(sd["id", ], length)
nNucl  = sapply(qd, nrow); nNucl[nNucl == 0] = NA
nDivs  = as.numeric(table(dd$t))

# Q: How do the number of division events change with the number of cells?
# A: Not at all
# Q: How do the number of nuclei change with the number of cells?
# A: Not at all
par(mfcol = c(3, 2))
plot(tp,         nCells,                type = "b")
plot(tp,         nNucl,                 type = "b")
plot(tp[-last],  nDivs,                 type = "b")
plot(tp[-last],  nDivs / nCells[-last], type = "b")
plot(tp[-last],  nDivs / nNucl[-last],  type = "b")
plot(tp,         nNucl / nCells,        type = "b")
cor.test(nDivs,  nCells[-last]) 
cor.test(nDivs,  nNucl[-last])  
cor.test(nCells, nNucl)         

# How do the number of active CLV3 nuclei change with the number of cells found? 
par(mfrow = c(3, 1))
plot(tp, nCells, type = "b")
plot(tp, nNucl, type  = "b")
###############


