#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/")
source("aux.R", local = TRUE)
source("read_data.R", local = TRUE)
source("process_data.R", local = TRUE)
source("main.R")

# plants    = c(1, 2, 4, 13, 15, 18)
plants = 18
map.it    = c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE)
# plants = 1
no.plants = length(plants)
no_cores  = detectCores() - 1
cl        = makeCluster(no_cores)
silent    = clusterExport(cl, ls())
results   = parLapply(cl, 1:length(plants), function(x) get.plant.data(plants[x], map.it[x]))
stopCluster(cl)

all.data = do.call(rbind, lapply(results, "[[", "all.data"))

### Cell divisions vs radius
no.breaks = 30
par(mfrow = c(1, 1))
hist((sqrt(all.data[, "x"] ** 2 + all.data[, "y"] ** 2)),    breaks = no.breaks, main = "rad")
hist(all.data[,"dist2top"], breaks = no.breaks, main = "rad")
hist(results[[2]][["all.data"]][,"dist2top"], breaks = no.breaks, main = "rad")
hist(results[[3]][["all.data"]][,"dist2top"], breaks = no.breaks, main = "rad")
hist(results[[4]][["all.data"]][,"dist2top"], breaks = no.breaks, main = "rad")
hist(results[[5]][["all.data"]][,"dist2top"], breaks = no.breaks, main = "rad")
# hist(results[[6]][["all.data"]][,"dist2top"], breaks = no.breaks, main = "rad")

# Plot all the quantification data
par(mfrow = c(4, 3))
no.breaks = 1000
plot.stat.hists(all.data, no.breaks=no.breaks, density = F)
for(ii in 6:ncol(all.data)) if(!all(is.na(all.data[,ii]))) plotDensity(data.frame(x=all.data[, ii]), main = colnames(all.data)[ii])
plant = 4
plot.stat.hists(results[[plant]]$all.data, no.breaks=no.breaks, density = F)
for(ii in 6:ncol(results[[plant]][["all.data"]])) if(!all(is.na(results[[plant]][["all.data"]][,ii]))) plotDensity(data.frame(x=results[[plant]][["all.data"]][, ii]), main = colnames(results[[plant]][["all.data"]])[ii])

plots = plot.all.stats(all.data)
plots[[1]]
plots[[2]]
plots[[3]]
plots[[4]]
plots[[5]]
plots[[6]]
plots[[7]]
plots[[8]]
plots[[9]]

plant = 4
plots = plot.all.stats(results[[plant]]$all.data)

##################################
# plot.12(cell.lines.data[order(-sapply(cell.lines.data, function(x) sum(nna(unlist(x)))))])

# Sort sublineages by mean dist2top
plant = 2
cell.lines.data = lapply(results[[plant]]$lineage.sublines.data, collapse.lineage) 
lineage.sublines.data = results[[plant]]$lineage.sublines.data[order(
  unlist(lapply(cell.lines.data, function(x) mean(x[,"dist2top"], na.rm=T))))] 
plot.12(unlist(lineage.sublines.data, recursive = F), batch = 1)

par(mfcol = c(2, 3), mar = c(2,2,2,2))
batch = 1; range  = (((batch-1)*3 + 1):(3*batch))
for(ii in range) plot.lineage(lineage.sublines.data[[ii]], var = "Mean.cell.intensity", ylim=c(100,160))
# for(ii in range) plot.lineage(lineage.sublines.data[[ii]], var = "dist2top")


# df = data.frame(ID = 1:12, do.call(rbind, lapply(subline.data[1:12], function(x) unlist(x[,"Mean.cell.intensity"]))))
# t  = data.frame(ID = 1:12, do.call(rbind, lapply(subline.data[1:12], function(x) unlist(x[,"Time"]))))


# library(traj)
# s1 = step1measures(df, t, ID = TRUE)
# s2 = step2factors(s1, discard = TRUE)
# s3 = step3clusters(s2, nclusters = 4)
# plant.data = results[[2]]$all.data
# t = plant.data[,"t"]; expr = plant.data[,"vol"]; there = nna(t) & nna(expr)
# p = ggplot(data.frame(t=t[there], expr=expr[there]), aes(factor(t), expr)); p + geom_violin()
# t = plant.data[,"t"]; expr = plant.data[,"expr"]; there = nna(t) & nna(expr)
# p = ggplot(data.frame(t=t[there], expr=expr[there]), aes(factor(t), expr)); p + geom_violin()
# t = plant.data[,"t"]; expr = plant.data[,"dist2top"]; there = nna(t) & nna(expr)
# p = ggplot(data.frame(t=t[there], expr=expr[there]), aes(factor(t), expr)); p + geom_violin()

