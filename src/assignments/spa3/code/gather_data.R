require(methods)
args    = commandArgs(trailingOnly=TRUE)
HOMEDIR = "/home/henrik/spa3/"
TSPDIR  = "tsp_data/"
setwd(HOMEDIR)
# ga.files = system("cat tsp_data/*ga_solution*",     intern=TRUE)
# sa.files = system("cat tsp_data/*simann_solution*", intern=TRUE)
# ga.files = unname(sapply(ga.files, function(x) strsplit(x, "\t")[[1]]))
# sa.files = unname(sapply(sa.files, function(x) strsplit(x, "\t")[[1]]))
# ga.files = t(ga.files)
# sa.files = t(sa.files)
# sa.data = data.frame(sa.files)
# ga.data = data.frame(ga.files)
# sa.data[,1] = sapply(sa.data[,1], function(x) strsplit(as.character(x), "/")[[1]][2])
# ga.data[,1] = sapply(ga.data[,1], function(x) strsplit(as.character(x), "/")[[1]][2])
# sa.data = sa.data[order(sa.data[,1]), ]
# ga.data = ga.data[order(ga.data[,1]), ]
# all.data = data.frame(Tour=sa.data[,1], Genalg=ga.data[,2], Simann=sa.data[,2], Optimal=sa.data[,3])
# all.data[,2] = as.character.factor(all.data[,2])
# all.data[,3] = as.character.factor(all.data[,3])
# all.data[,4] = as.character.factor(all.data[,4])
# all.data[,2] = as.numeric(all.data[,2])
# all.data[,3] = as.numeric(all.data[,3])
# all.data[,4] = as.numeric(all.data[,4])
# all.data[,-1] = all.data[,-1] / all.data[,4]
# par(mar=c(10,7,7,6))
# barplot(t(as.matrix(all.data[,2:3])), beside=T, 
#         col=c("aquamarine3","coral"),
#         names.arg = all.data[,1], las=3,  cex.axis = 1.4, cex.names=1.4)
# legend("topleft", c("Genetic algorithm","Simulated annealing"), pch=15, 
#        col=c("aquamarine3","coral"), 
#        bty="n", cex=2)
# 
# file = system("find tsp_data/ -name *.opt.tour",intern=T)
# normal = vector(mode="logical")
# for(ii in gsub("\\.opt\\.tour", "\\.tsp", file)){
# a = any(grepl("EUC_2D", readLines(ii)))
# b = any(grepl("NODE_COORD_SECTION", readLines(ii))) 
# c = a & b
# normal = append(normal, c)
# # normal = append(normal, any(grepl("NODE_COORD_SECTION", readLines(ii))))
# }
# file = file[normal]
# write(file, "target_tours.txt", sep = "\n")

read.coord = function(f){
  lines = readLines(f)
  lines = trim(lines)
  lines = lines[which(lines != "EOF")]
  lines = lines[which(sapply(lines, length) > 0)]
  first = grep("^1\\s+|^0001\\s+", lines)[1]
  lines = lines[-(1:(first-1))]
  lines = sapply(lines, function(x) strsplit(x, "\\s+")[[1]])
  lines = unlist(lines)
  lines = lines[which(lines != "")]
  data  = matrix(lines, byrow=TRUE, ncol=3)
  # print(data)
  data  = data.frame(data)
  data
}

read.opt = function(f){
  #  print(f)
  lines = readLines(f)
  lines = trim(lines)
  lines = lines[which(lines != "EOF")]
  lines = lines[which(sapply(lines, length) > 0)]
  first = grep("^1", lines)[1]
  lines = lines[-(1:(first-1))]
  lines = sapply(lines, function(x) strsplit(x, "\\s+")[[1]])
  lines = unlist(lines)
  lines = lines[which(lines != "")]
  lines = as.numeric(lines)
  lines = lines[which(lines > 0)]
  lines
}

