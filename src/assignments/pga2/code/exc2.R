setwd("~/compbio/src/assignments/pga2/code/")
.libPaths("/home/henrik/R/x86_64-pc-linux-gnu-library/3.3/")
library("RColorBrewer")
library(scales) #imports alpha
palette(brewer.pal(n = 8, name = "Set1"))
line.size = 3
lw.s = 3


# 1) ./analyse_first_multi <seed> <no mut> <no optimis> <input file>
# setwd("~/.local/optimist")
# no.mut = 1
# system(paste("./analyse_first_multi 2", no.mut, "10 Colour_freq.out"))
# system("rm Transfer_data.out")
# system(paste0("head -n ", 5 + no.mut, " Transfer_data_multi.out > Transfer_data.out"))
# system("./analyse_sys_multi 2 Colour_freq.out Transfer_data.out")
# system("cp Real_frequencies.out ~/compbio/src/assignments/pga2/data/Real_frequencies1.out")
# system("cp Model_frequencies.out ~/compbio/src/assignments/pga2/data/Model_frequencies1.out")

# no.mut = 2
# system(paste("./analyse_first_multi 2", no.mut, "10 Colour_freq.out"))
# system("rm Transfer_data.out")
# system(paste0("head -n ", 5 + no.mut, " Transfer_data_multi.out > Transfer_data.out"))
# system("./analyse_sys_multi 2 Colour_freq.out Transfer_data.out")
# system("cp Real_frequencies.out ~/compbio/src/assignments/pga2/data/Real_frequencies2.out")
# system("cp Model_frequencies.out ~/compbio/src/assignments/pga2/data/Model_frequencies2.out")
# 
# # 1) ./analyse_first_multi <seed> <no mut> <no optimis> <input file>
# no.mut = 3
# system(paste("./analyse_first_multi 2", no.mut, "10 Colour_freq.out"))
# system("rm Transfer_data.out")
# system(paste0("head -n ", 5 + no.mut, " Transfer_data_multi.out > Transfer_data.out"))
# system("./analyse_sys_multi 2 Colour_freq.out Transfer_data.out")
# system("cp Real_frequencies.out ~/compbio/src/assignments/pga2/data/Real_frequencies3.out")
# system("cp Model_frequencies.out ~/compbio/src/assignments/pga2/data/Model_frequencies3.out")
# setwd("~/compbio/src/assignments/pga2/code/")

# #################################### 
# ### PLOT
# ####################################
# 
d1 = read.table("../data/Model_frequencies1.out")
d2 = read.table("../data/Real_frequencies1.out")
d3 = read.table("../data/Model_frequencies2.out")
d4 = read.table("../data/Real_frequencies2.out")
d5 = read.table("../data/Model_frequencies3.out")
d6 = read.table("../data/Real_frequencies3.out")
#

# Log-likelihoods produced.
L1 = c( -92.5714,-92.5714,-92.5714,-92.5714,-92.5714,-92.5714,-92.5714,-92.5714,-92.5714,-92.5714)
L2 = c(-81.7473,-81.7689,-81.75,-81.7561,-81.8105,-81.7567,-81.7759,-81.7451,-81.7841,-81.8126 )
L3 = c(-81.1746,-81.7301,-81.0976,-81.4593,-81.3007,-81.3259,-81.7453,-81.356,-81.5501,-81.4329)

AIC = function(params, log.L){
  2*params - 2*log.L
}

BIC = function(n, k, log.L){
  log(n)*k - 2*log.L
}

Ls = rbind(L1, L2, L3)
means  = apply(Ls, 1, mean)
stds = apply(Ls, 1, function(x) sd(x)/sqrt(length(x)))

aics = AIC(1:3, means)
bics = BIC(10, 1:3, means)
aics.err.up = AIC(1:3, means + stds)
aics.err.do = AIC(1:3, means - stds)
bics.err.up = BIC(10, 1:3, means + stds)
bics.err.do = BIC(10, 1:3, means - stds)

# PLOT
par(mfrow=c(3,1))
plot(d1[,1], d1[,2], col = 1, type = "l", pch = 16, ylim=0:1, lwd = line.size, lty = 1)
points(d1[,1], d1[,3], col = 2, type = "l", pch = 17, lwd = line.size, lty = 1)
lines(d2[,1], d2[,2], lwd = line.size, type = "p", col = alpha(1, .5))
lines(d2[,1], d2[,3], lwd = line.size, type = "p", col = alpha(2, .5))

plot(d3[,1], d3[,2], col = 1, type = "l", pch = 16, ylim=0:1, lwd = line.size, lty = 1)
points(d3[,1], d3[,3], col = 2, type = "l", pch = 17, lwd = line.size, lty = 1)
lines(d4[,1], d4[,2], lwd = line.size, type = "p", col = alpha(1, .5))
lines(d4[,1], d4[,3], lwd = line.size, type = "p", col = alpha(2, .5))

plot(d5[,1], d5[,2], col = 1, type = "l", pch = 16, ylim=0:1, lwd = line.size, lty = 1)
points(d5[,1], d5[,3], col = 2, type = "l", pch = 17, lwd = line.size, lty = 1)
lines(d6[,1], d6[,2], lwd = line.size, type = "p", col = alpha(1, .5))
lines(d6[,1], d6[,3], lwd = line.size, type = "p", col = alpha(2, .5))



par(mfrow=c(1,1))
plot(1:3, aics, col = 1, type = "b", lwd = lw.s, xlab = "", ylab = "")
lines(1:3, bics, col = 2, type = "b", lwd = lw.s)
arrows(1:3, aics.err.up, 1:3, aics.err.do, length = 0.05, angle = 90, code = 3, col = "1", lwd = lw.s)
arrows(1:3, bics.err.up, 1:3, bics.err.do, length = 0.05, angle = 90, code = 3, col = "2", lwd = lw.s)
mtext(side=2, line = 3, las = 3, "Score")
mtext(side=1, line = 3, las = 1, "Number of mutations")
legend("topright", col =1:2, lty = 1, pch = 21, c("AIC", "BIC"), inset = c(0.02,0.02), lwd = lw.s)
# 2 beneficial mutations!
