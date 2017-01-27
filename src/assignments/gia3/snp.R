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

