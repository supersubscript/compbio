BIOHOME = system(paste0("find ", path.expand("~"), "/", 
                        " -type d -name 'compbio'"), 
                        intern = TRUE)
setwd(paste0(BIOHOME, "/text/functional_genomics/assignment_1/"))
.libPaths(c(.libPaths(), "/local/data/public/hpa22/R/x86_64-pc-linux-gnu-library/3.1/"))
# source("https://bioconductor.org/biocLite.R")
# biocLite("affy")
# biocLite(affy)
# biocLite(affyPLM)
# biocLite(simpleaffy)
# biocLite(Harshlight)
# biocLite(limma)
# biocLite(lumi)


library(affy)
library(affyPLM)
library(simpleaffy)
library(Harshlight)
library(limma)
library(lumi)
library(panp)

Data  = ReadAffy(celfile.path = "./data/")
harsh = Harshlight(Data, na.sub = FALSE)
#eset  = mas5(harsh) ### Does log2 normalisation
#eset = expresso(harsh, bg.correct = bg.correct.rma(harsh))
eset <- expresso(harsh, normalize.method="quantiles.robust",
                 bgcorrect.method="rma", pmcorrect.method="pmonly",
                 summary.method="avgdiff")
# Expression set -- value for a single probe
### RMA does log2, mas5 not

### Set estrogen label
er.label = readLines("./data/GSE51450_series_matrix.txt", n = 37)[37]
er.label = strsplit(er.label, "\t")[[1]][-1]
er.label = gsub('\\\"', '', er.label)
er.label = ifelse(grepl("positive", er.label), "pos", "neg")
er.label = factor(er.label)
 
#### Histogramaflams
hist(Data, xlim=c(0, 20))
affy::plotDensity(log2(exprs(eset)), xlim = c(0, 20))

### Boxaplottaslots
boxplot(log2(exprs(Data)))
boxplot(log2(exprs(eset)))

### Quality controlitronix
Data.qc = qc(harsh, eset)
plot(Data.qc)

### 
affy::MAplot(harsh, plot.method = "smoothScatter")
affy::MAplot(eset,plot.method = "smoothScatter", log = TRUE)

# PCA
pca.unnorm = prcomp(exprs(harsh))[[2]][, 1:2]
pca.norm   = prcomp(exprs(eset))[[2]][, 1:2]
par(mfrow=c(1,2))
plot(pca.unnorm, cex = 3, pch = 17, col = as.integer(er.label), 
     xlim = c(-.38,-.22), ylim = c(-.6,.7))
plot(pca.norm, cex = 3, pch = 17, col = as.integer(er.label), 
     xlim = c(-.38,-.22), ylim = c(-.6,.7))

### FILTERING
### Remove if found in less than 6/12. 
PA = pa.calls(eset, looseCutoff = 0.1, tightCutoff = 0.1)
absent = rowSums(PA$Pcalls == 'A')
absent = which(absent >= ncol(PA$Pcalls) / 2) 
filtered.eset = eset[-absent, ]

pca.filtered.norm   = prcomp(exprs(filtered.eset))[[2]][, 1:2]
plot(pca.filtered.norm, cex = 3, pch = 17, col = as.integer(er.label))

affy::MAplot(filtered.eset, plot.method = "smoothScatter")
affy::plotDensity(log2(exprs(filtered.eset)))
boxplot(log2(exprs(filtered.eset)))

### Differential expressionz

par(mfrow=c(1,2))
conditions       = factor(c(rep(0,6), rep(1,6)), labels = c("pos","neg"))
design           = model.matrix(~ conditions - 1)
colnames(design) = levels(conditions)
fit              = lmFit(log2(exprs(eset)), design)  
contrast.matrix  = makeContrasts("pos-neg", levels = design)
cont.fit         = contrasts.fit(fit, contrast.matrix)
cont.fit         = eBayes(cont.fit)
volcanoplot(cont.fit,ylim=c(-5,-2))

fit              = lmFit(log2(exprs(filtered.eset)), design)  
contrast.matrix  = makeContrasts("pos-neg", levels = design)
cont.fit         = contrasts.fit(fit, contrast.matrix)
cont.fit         = eBayes(cont.fit)
volcanoplot(cont.fit,ylim=c(-5,-2))



### Any finds? 
topTable(cont.fit, n = 5)

### Formaline fixed, blalblbal -- expect high triangles
### GAPDH -- housekeeping expression, should be around 1 