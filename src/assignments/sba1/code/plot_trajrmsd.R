setwd("~/compbio/src/assignments/sba1/code/")
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))


filenames = list.files("../data", pattern="*.dat", full.names=TRUE)
res = lapply(filenames, function(x) read.table(x, header = T))
names(res) = filenames #substr(filenames, 1, 30)
plot(1, type="n", axes=T, xlab="", ylab="", ylim=c(0,5), xlim = c(0,500))
for(ii in 1:length(res)){
  lines(res[[ii]], col = ii, lwd = 2)
}
conf.names = c("+ N CA C", "+ N CA C and helix", "+ N CA C and sheet", "+ N CA C and (sheet or helix)", "Protein")
# legend("bottomright", names(res), col = 1:length(res), inset = c(0.01, 0.05), lty = 1, lwd = 2)

legend("bottomright", conf.names, col = 1:length(res), inset = c(0.01, 0.05), lty = 1, lwd = 2)
