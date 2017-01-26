#!/usr/bin/Rscript
#source("https://bioconductor.org/biocLite.R")
#biocLite("DiffBind")
library(DiffBind)
options(stringsAsFactors = FALSE)
HOMEDIR = "/local/data/public/hpa22/assignments/fga3/"
setwd(HOMEDIR)

INPUT    = paste0(HOMEDIR, "names/sample_sheet_chip1.csv")
data     = dba(sampleSheet=INPUT)
onlyT47D = dba.mask(data, attribute = DBA_TISSUE, value = 'T47D')
onlyMCF7 = dba.mask(data, attribute = DBA_TISSUE, value = 'MCF7')

dba.plotHeatmap(data, ColAttributes = c(DBA_FACTOR, DBA_TREATMENT),
                mask = onlyT47D, bUsePval = FALSE, dendogram = 'none', 
                main = 'T-47D cells')

dba.plotPCA(data, DBA_TISSUE,label=DBA_CONDITION, main = 'T-47D cells')

# Cluster heatmap
plot(data)

# Count reads with peaks extending 5kb up- and downstream 
data = dba.count(data, summits=5000)

# Plot correlation heatmap for read counts (not included in paper, right?)
plot(data)

# Create contrast to establish which cell lines fall in which groups
data = dba.contrast(data, categories="DBA_CONDITION")

# Do the differential analysis
data = dba.analyze(data)

# Heatmap for differentially bound sites
plot(data, contrast=1)
save(a)
# Retrieve the differentially bound regions
# This should be exported and piped through bedtools in order
# to get count data. (To get Figure 2B.)
# $ bedtools multicov -bams aln1.bam aln2.bam aln3.bam -bed ivls-of-interest.bed
data.DB = dba.report(data)
df      = data.frame(seqnames=seqnames(data.DB),
                     starts=start(data.DB)-1,
                     ends=end(data.DB),
                     names=c(rep(".", length(data.DB))),
                     scores=c(rep(".", length(data.DB))),
                     strands=strand(data.DB))
write.table(df,file="thatFile.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# Make MA-plot (Figure 3C)
dba.plotMA(data)

# Do PCA (Extended Figure 2B; TODO: Both cell lines, all conditions -- don't get how this works)
dba.plotPCA(data, DBA_TISSUE,label=DBA_CONDITION)
