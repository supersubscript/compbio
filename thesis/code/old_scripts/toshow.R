setwd("/home/henrik/compbio/thesis/code/")
source("aux.R",          local = TRUE)
source("read_data.R",    local = TRUE)
source("process_data.R", local = TRUE)
source("plot.R")
source("main.R")

plants          = c(1, 2, 4, 13, 15) # 18
quant.missing   = list(c(-1), c(-1), c(24), c(-1), c(12))
mapping.missing = list(c(-1), c(-1), c(-1), c(-1), c(-1))
no.plants       = length(plants)

no_cores  = detectCores() - 1
cl        = makeCluster(no_cores)
silent    = clusterExport(cl, ls())
results   = parLapply(cl, 1:no.plants, function(x) 
  get.plant.data(plants[x], 
                 quant.missing = quant.missing[[x]],
                 map.quant     = TRUE,
                 only.l1       = FALSE,
                 only.l2       = FALSE))
silent    = stopCluster(cl)
all.data  = do.call(rbind, lapply(results, "[[", "all.data"))

### Cell divisions at certain distance to top. 
# Problem: Analysis obstructed by different distances (and noise) between plants
par(mfrow = c(3,2))
no.breaks = 30
hist(all.data[,"dist2top"], breaks = no.breaks, main = "All")
hist(results[[2]][["all.data"]][,"dist2top"], breaks = no.breaks, main = "2")
hist(results[[3]][["all.data"]][,"dist2top"], breaks = no.breaks, main = "4")
hist(results[[4]][["all.data"]][,"dist2top"], breaks = no.breaks, main = "13")
hist(results[[5]][["all.data"]][,"dist2top"], breaks = no.breaks, main = "15")

### Statistics distribution histograms @ DIVISION
#   1) AGE:   Initial age values are a defect. They come from cells dropping in / out of the tracking
#   2) TIME:  Is the bumpiness in time a consequence of the division time, or a signifier of division commitments? TODO: Correlate with no.nuclei.
#   3) VOL:   Bumpiness? Hard to tell. In that case, link to spatial structure? Volume decays with distance to top.
#   4) EXP:   Bumpiness. Link to spatial structure? Signifier of higher concentration compared to neighbours? 
#   5) D2T:   Idk. No clear patterns due to noise. Multiple peaks due to avg cell width?
#   6) Dvol:  Peak at around .3 of mother cell overall. Per plant. Bumpiness -- splits into fractions? 
#   7) Dexp:  Two peaks in interval 0-1? 
par(mfrow = c(3, 3))
plot.stat.hists(all.data, 30)
plot.stat.hists(results[[1]]$all.data, 30)
plot.stat.hists(results[[2]]$all.data, 30)
plot.stat.hists(results[[3]]$all.data, 30)
plot.stat.hists(results[[4]]$all.data, 30)
plot.stat.hists(results[[5]]$all.data, 30)

### Cell lines as whole, sorted by distance from the top
#   1) All plants show tendencies to oscillate -- not only at top, but clearer due to no decay
#   2) Thresholding effect: Decay seems to initiate at around ~2 from the top? See 1.3, 2, 4.2, 4.3 for plant 2
#   3) Cells alternate with who's on top
plant = 2; cell.lines.data = lapply(results[[plant]]$lineage.sublines.data, collapse.lineage) 
lineage.sublines.data = results[[plant]]$lineage.sublines.data[order(unlist(lapply(cell.lines.data, function(x) mean(x[,"dist2top"], na.rm=T))))] 
lineage.sublines.data = lineage.sublines.data[sapply(lineage.sublines.data, function(x) sum(nna(x))) > 6*8]
plot.12(unlist(lineage.sublines.data, recursive = F), batch = 1)
par(mfcol = c(2, 3)); for(ii in 10:12) plot.lineage(lineage.sublines.data[[ii]], ylim = c(0, 160))

### Naive attempt to classify derivate as a function of distance to top
all.results = do.call(rbind, lapply(lineage.sublines.data, get.roc))
par(mfrow = c(1, 1))
plot(all.results[, 1], all.results[, 2], ylim = c(-1.5, 1.5))
cor.test(all.results[, 1], all.results[, 2])

###############################################################################
### Population scale analysis
###############################################################################
### Number of CLV3 nuclei. Differences between plants.
# par(mfrow = c(3,2))
# for(ii in 1:nrow(no.clv3.nuclei)) plot(timepoints[ii,], no.clv3.nuclei[ii,], type = "b", main = paste0("no.clv3.nuclei ", plants[ii]), xlim=c(0,84))
# plot(seq(0,84,4), mean.no.clv3.nuclei, type = "b", main = "mean.no.clv3.nuclei all")

# Mean, a bit more grained out
par(mfrow = c(2,1))
plot(seq(0,84,4), mean.no.clv3.nuclei, type = "b", main = "mean.no.clv3.nuclei")
arrows(seq(0,84,4), mean.no.clv3.nuclei - sd.no.clv3.nuclei / sqrt(apply(no.clv3.nuclei, 2, function(x) sum(nna(x)))),
       seq(0,84,4), mean.no.clv3.nuclei + sd.no.clv3.nuclei / sqrt(apply(no.clv3.nuclei, 2, function(x) sum(nna(x)))), length = 0.05, angle = 90, code = 3)
plot(seq(0,84,4), sd.no.clv3.nuclei, type = "b", main = "sd.no.clv3.nuclei")
# plot(seq(0,84,4), sd.no.clv3.nuclei / mean.no.clv3.nuclei, type = "b", main = "sd.no.clv3.nuclei / mean.no.clv3.nuclei")

### Mean and SD varies significantly between plants
par(mfrow = c(3, 2)); lapply(1:length(clv3.mean), function(x) plot(clv3.mean[[x]], type = "b", main = paste0("clv3.mean ", plants[x])))
par(mfrow = c(3, 2)); lapply(1:length(clv3.sd),   function(x) plot(clv3.sd[[x]],   type = "b", main = paste0("clv3.sd ",   plants[x])))
par(mfrow = c(3, 2)); lapply(1:length(vol.mean),  function(x) plot(vol.mean[[x]],  type = "b", main = paste0("vol.mean ",  plants[x])))
par(mfrow = c(3, 2)); lapply(1:length(vol.sd),    function(x) plot(vol.sd[[x]],    type = "b", main = paste0("vol.sd ",    plants[x])))

