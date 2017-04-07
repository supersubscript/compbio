setwd("~/compbio/src/assignments/pga2/code/")
library(RColorBrewer)
library(scales) #imports alpha
palette(brewer.pal(n = 8, name = "Set1"))
line.size = 2

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
# d1 = read.table("../data/Model_frequencies1.out")
# d2 = read.table("../data/Real_frequencies1.out")
# d3 = read.table("../data/Model_frequencies2.out")
# d4 = read.table("../data/Real_frequencies2.out")
# d5 = read.table("../data/Model_frequencies3.out")
# d6 = read.table("../data/Real_frequencies3.out")
# 
# 
# par(mfrow=c(3,1))
# 
# # Run0 Log L = -92.5714
# # Run1 Log L = -92.5714
# # Run2 Log L = -92.5714 
# # Run3 Log L = -92.5714
# # Run4 Log L = -92.5714
# # Run5 Log L = -92.5714
# # Run6 Log L = -92.5714
# # Run7 Log L = -92.5714
# # Run8 Log L = -92.5714
# # Run9 Log L = -92.5714
# plot(d1[,1], d1[,2], col = 1, type = "l", pch = 16, ylim=0:1, lwd = line.size, lty = 1)
# points(d1[,1], d1[,3], col = 2, type = "l", pch = 17, lwd = line.size, lty = 1)
# lines(d2[,1], d2[,2], lwd = line.size, type = "p", col = alpha(1, .5))
# lines(d2[,1], d2[,3], lwd = line.size, type = "p", col = alpha(2, .5))
# 
# # Run0 Log L = -81.7473
# # Run1 Log L = -81.7689
# # Run2 Log L = -81.75
# # Run3 Log L = -81.7561
# # Run4 Log L = -81.8105
# # Run5 Log L = -81.7567
# # Run6 Log L = -81.7759
# # Run7 Log L = -81.7451
# # Run8 Log L = -81.7841
# # Run9 Log L = -81.8126
# plot(d3[,1], d3[,2], col = 1, type = "l", pch = 16, ylim=0:1, lwd = line.size, lty = 1)
# points(d3[,1], d3[,3], col = 2, type = "l", pch = 17, lwd = line.size, lty = 1)
# lines(d4[,1], d4[,2], lwd = line.size, type = "p", col = alpha(1, .5))
# lines(d4[,1], d4[,3], lwd = line.size, type = "p", col = alpha(2, .5))
# 
# # Run0 Log L = -81.1746
# # Run1 Log L = -81.7301
# # Run2 Log L = -81.0976
# # Run3 Log L = -81.4593
# # Run4 Log L = -81.3007
# # Run5 Log L = -81.3259
# # Run6 Log L = -81.7453
# # Run7 Log L = -81.356
# # Run8 Log L = -81.5501
# # Run9 Log L = -81.4329
# plot(d5[,1], d5[,2], col = 1, type = "l", pch = 16, ylim=0:1, lwd = line.size, lty = 1)
# points(d5[,1], d5[,3], col = 2, type = "l", pch = 17, lwd = line.size, lty = 1)
# lines(d6[,1], d6[,2], lwd = line.size, type = "p", col = alpha(1, .5))
# lines(d6[,1], d6[,3], lwd = line.size, type = "p", col = alpha(2, .5))

# 2 beneficial mutations!
