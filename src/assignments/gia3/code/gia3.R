# Summary:	It has generally been assumed that most differences between males and females
# are due to developmental and hormonal differences between the sexes. Here we investigate
# the contribution of sex chromosomal complement to such sexual dimorphisms. These genome-wide
# transcription profiling showed that the expression of hundreds of autosomal genes was sensitive
# to sex chromosome complement, rather than gender. The existence of such differences between
# males and females holds important implications for understanding sexual dimorphisms in
# physiology and disease hitherto attributed solely to gender or hormonal effects.

# Overall design:	Thymic total RNA was isolated from 7-8 week old mice, with 3 biological
# replicates for each of four genotypes with different sex chromosome complements.

# Libs:
source("http://www.bioconductor.org/biocLite.R")
#biocLite("GEOquery")
#biocLite("affyPLM")
#biocLite("mouse4302.db")
#biocLite("limma")
#biocLite("edgeR")
#biocLite("affyQCreport")
# librar(edgeR)
# library(affyQCreport)
library(limma)
library(Vennerable)
library(GEOquery)
library(Biobase)
library(affyPLM)
library(mouse4302.db)
library(annotate)
library(simpleaffy)

# Read in preprocessed data.
# gset <- getGEO("GSE21822", GSEMatrix = TRUE, getGPL = FALSE)
# if (length(gset) > 1)
#   idx <- grep("GPL1261", attr(gset, "names"))
# else
#   idx <- 1
# gset <- gset[[idx]]
# 
# # Plot distribution of data.
# dev.new(width = 4 + dim(gset)[[2]] / 5, height = 6)
# par(mar = c(2 + round(max(
#   nchar(sampleNames(gset))
# ) / 2), 4, 2, 1))
# title = paste ("GSE21822", '/', annotation(gset), " selected samples", sep = '')
# boxplot(
#   exprs(gset),
#   boxwex = 0.7,
#   notch = TRUE,
#   main = title,
#   outline = FALSE,
#   las = 2
# )
# Raw data
#filePaths = getGEOSuppFiles("GSE21822")
setwd("~/gia3")
data = ReadAffy(celfile.path = paste0(getwd(), "/../GSE21822"))
pData(phenoData(data)) = cbind(
  pData(phenoData(data)),
  group = c(
    "XX",
    "XX",
    "XX",
    "X[Y-]sry",
    "X[Y-]sry",
    "X[Y-]sry",
    "XXsry",
    "XXsry",
    "XXsry",
    "X[Y-]",
    "X[Y-]",
    "X[Y-]"
  )
)
eset = rma(data) # Returns log2 values.
plotDensity(eset)

# MA-plot
# par(mfrow = c(1, 1))
# colnames(exprs(eset)) = cols
# MAplot(eset, plot.method = "smoothScatter")

# Female comparison
results.female = pairwise.comparison(
  eset, ## processed data
  "group", ## the factor in covdesc
  ## to use to group arrays
  c("XX", "X[Y-]"), ## groups to compare
  data, ## for PMA calls,
  method='unlogged',
  logged = TRUE
)
significant.female = pairwise.filter(
  results.female,
  min.present.no = 1,
  present.by.group = TRUE,
  fc = log2(1.2),
  tt = 0.05
)
# Male comparison
results.male = pairwise.comparison(
eset, ## processed data
"group",
c("XXsry", "X[Y-]sry"), ## groups to compare
data, ## for PMA calls,
method='unlogged',
logged = TRUE)
significant.male = pairwise.filter(
results.male,
min.present.no = 1,
present.by.group = TRUE,
fc = log2(1.2),
tt = 0.05
)

# Which probes are found?
probes.male   = means(significant.male)
probes.female = means(significant.female)
scs.genes     = intersect(rownames(probes.male), rownames(probes.female))
nrow(probes.male)
nrow(probes.female)
length(scs.genes) # Genes sensitive to X chromosome dosage in BOTH males and females. 

# # Plot 'em
# par(mfrow = c(1, 2))
# plot(results.male, showPMA = TRUE, main = "Male comparison")
# lines(0:40, 0:40 + .2, lw = 2)
# lines(0:40, 0:40 - .2, lw = 2)
# points(means(results.male)[which(rownames(means(results.male)) %in% scs.genes), ],
#   col = 'black',
#   pch = 16,
#   cex = 1)
# plot(results.female, showPMA = TRUE, main = "Female comparison")
# lines(0:40, 0:40 + .2, lw = 2)
# lines(0:40, 0:40 - .2, lw = 2)
# points(means(results.female)[which(rownames(means(results.female)) %in% scs.genes), ],
#   col = 'black',
#   pch = 16,
#   cex = 1)
# 
# plot(1, type="n", axes=F, xlab="", ylab="")
# venn.plot = venn.diagram(list(rownames(probes.male), rownames(probes.female)),NULL,fill = c("green", "red"),alpha = c(0.5, 0.5),cex = 3,cat.fontface = 1,category.names = c("Male", "Female"),cat.cex = 1.5)
# grid.draw(venn.plot)

#########
results.dimorphic = pairwise.comparison(eset, ## processed data
"group", ## the factor in covdesc
## to use to group arrays
c("XX", "X[Y-]sry"), ## groups to compare
data, ## for PMA calls,
method = 'unlogged',
logged = TRUE)

significant.dimorphic = pairwise.filter(
results.dimorphic,
min.present.no = 1,
present.by.group = TRUE,
fc = log2(1.2),
tt = 0.05
)
# plot(results.dimorphic)
# lines(0:40, 0:40 + .2, lw = 2)
# lines(0:40, 0:40 - .2, lw = 2)
# points(means(results.female)[which(rownames(means(results.female)) %in% scs.genes),],
#        col = 'black',
#        pch = 16,
#        cex = 1)

dimorphic.regulated.by.cc = intersect(scs.genes, rownames(means(significant.dimorphic)))
dimorphic.in.either.sex   = union(intersect(rownames(means(significant.dimorphic)), rownames(means(significant.male))), intersect(rownames(means(significant.dimorphic)), rownames(means(significant.female))))

length(dimorphic.regulated.by.cc)
length(dimorphic.in.either.sex)

getResults = function(x, y){
  pairwise.comparison(
    eset, ## processed data
    "group", ## the factor in covdesc
    c(x, y), ## groups to compare
    data, ## for PMA calls,
    method = 'unlogged',
    logged = TRUE)
}

results.xy      = getResults("X[Y-]", "X[Y-]sry")
results.xx      = getResults("XX", "XXsry")
results.xxxysry = getResults("XX", "XXsry")

plotPairwiseDiffExp = function(x)
{
  plot(x@means[,1], x@means[,2], cex = .4, col = 'gray', xlab = colnames(x@means)[1], ylab = colnames(x@means)[2], pch = 15)
  lines(0:40, 0:40 + .2, lw = 2, lt = 2)
  lines(0:40, 0:40 - .2, lw = 2, lt = 2)
  points(means(x)[which(rownames(means(x)) %in% scs.genes), ], col = 'black',pch = 16,cex = .7)
  points(means(x)[which(rownames(means(x)) %in% dimorphic.regulated.by.cc), ],col = 'green',pch = 17,cex = .7)
  points(means(x)[which(rownames(means(x)) %in% c(scs.genes, dimorphic.regulated.by.cc)), ], col = 'blue',pch = 16,cex = .7)
}

par(mfrow = c(3, 2), mar = c(4, 4, 1, 1))
plotPairwiseDiffExp(results.male)
text(3,12, "A", cex = 2)
plotPairwiseDiffExp(results.female)
text(3,12, "B", cex = 2)
plotPairwiseDiffExp(results.dimorphic)
text(3,12, "C", cex = 2)
#par(mfrow = c(2,2))
plotPairwiseDiffExp(results.xy)
text(3,12, "D", cex = 2)
plotPairwiseDiffExp(results.xx)
text(3,12, "E", cex = 2)
plotPairwiseDiffExp(results.xxxysry)
text(3,12, "F", cex = 2)

# Venn 1
males.length        = length(rownames(means(significant.male)))
females.length      = length(rownames(means(significant.female)))
mf.intersect.length = length(intersect(rownames(means(significant.male)), rownames(means(significant.female))))
plot(Venn(Sets = list(rownames(means(significant.male)), rownames(means(significant.female))), SetNames = c("Male", "Female")))

# Venn 2
scs.length          = length(scs.genes)
dimorphic.length    = length(rownames(means(significant.dimorphic)))
sd.intersect.length = length(intersect(scs.genes, rownames(means(significant.dimorphic))))
plot(Venn(Sets = list(scs.genes, rownames(means(significant.dimorphic))), SetNames = c("SCS", "Dimorphic")))

  # Actual gene names
# unique(select(
#   mouse4302.db,
#   rownames(probes.male),
#   c("SYMBOL", "ENTREZID", "GENENAME")
# )$ENTREZID)

# design = model.matrix( ~ group)
# fit    = lmFit(present.exprs, design)
# ebayes = eBayes(fit)
# genes.diffexp = which((ebayes$p.value)[,1:2] < .05)
# genes.diffexp
# significant.xy = pairwise.filter(
#   results.xy,
#   min.present.no = 1,
#   present.by.group = TRUE,
#   fc = log2(1.2),
#   tt = 0.05
# )
# significant.xxxysry = pairwise.filter(
#   results.xxxysry,
#   min.present.no = 1,
#   present.by.group = TRUE,
#   fc = log2(1.2),
#   tt = 0.05
# )
# significant.xx = pairwise.filter(
#   results.xx,
#   min.present.no = 1,
#   present.by.group = TRUE,
#   fc = log2(1.2),
#   tt = 0.05
# )


#### SNP stuffs
# http://www.ensembl.org/Homo_sapiens/Gene/Variation_Gene/Table?db=core;g=ENSG00000184895;r=Y:2786855-2787699;t=ENST00000383070

data  = read.table("../ensembl-export.csv", header = TRUE, sep = ",")
data  = data[which(data[, "Source"] == "SNP"),]
data  = data[grepl("Phenotype", data[,"Clin..Sig."]), ]
start = 2786855
stop  = 2787699

plot(seq(start, stop), rep(1, length(seq(start, stop))), type = "l", yaxt = "n", bty = "n", ylab = "", cex = 5, xaxt = "n", lwd = 3, xlab = "")
axis(1, line = -15, xlab = "Chromosome coordinate")
#lines(seq(2787150, 2787480), rep(1, length(seq(2787150, 2787480))), lwd = 10)
rect(2787150, .98, 2787480, 1.02, lwd = 3)
points(data[, "Chr"],
       rep(1, nrow(data)),
       type = "p",
       yaxt = "n", pch = 18, col = 'red', cex = 3)
points(data[, "Chr"],
       rep(1, nrow(data)),
       type = "p",
       yaxt = "n", pch = 5, col = 'black', cex = 2.5, lw = 2)
text(start + 75, 1.02, "SRY-001", cex = 1.5)
text(2787150 + 100, 1.04, "HMG-box", cex = 1.5)
