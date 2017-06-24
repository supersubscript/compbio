###############################################################################
### GENERAL IMPORT AND PLOT SETTINGS
###############################################################################
setwd("/home/henrik/compbio/thesis/code/")
library(parallel)
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 2
par(lwd = lw.s, cex = 1.5, ps = 15)

###############################################################################
### AUXILLIARY FUNCTIONS
###############################################################################
order.na = function(frame.list) {
  order(unlist(lapply(frame.list, function(x)
    sum(is.na(
      unlist(x)
    )))))
}

distance = function(v1, v2) {
  sqrt(sum((as.vector(v1) - as.vector(v2)) ** 2))
}

order.by.columns = function(data) data[do.call(order, lapply(1:NCOL(data), function(ii) data[, ii])), ]
order.by.rows    = function(data) data[do.call(order, lapply(1:NROW(data), function(ii) data[ii, ])), ]

init.counter = function() {
  x = 0
  function() {
    x <<- x + 1
    x
  }
}

get.lineages = function(sublines){
  subline.ids = rownames(sublines)
  lineage.ids = factor(sapply(subline.ids, function(x)
    strsplit(x, "\\.")[[1]][1]))
  return(split(as.data.frame(sublines), lineage.ids))
}

get.subline = function(sublines, line.id, all.lineages = NULL) {
  if (!is.null(all.lineages))
    return(all.lineages[line.id])
  else
    return(get.lineages(sublines)[line.id])
}

# Is not NA
nna = function(data){!is.na(unlist(data))}

### Extract numbers from string
extract.numbers = function(strs) {
  unname(sapply(strs, function(str) {
    matches = gregexpr('[0-9]+', str)
    sapply(regmatches(str, matches), as.numeric)
  }))
}


###############################################################################
### CUSTOM PLOT FUNCTIONS
###############################################################################
plot.lineage = function(lineage.sublines.data){
  for(ii in lineage.sublines.data) lines(ii[,"Time"], ii[,"Mean.cell.intensity"], type = "b", lwd = lw.s) 
}

###############################################################################
### HELPER PLOT FUNCTIONS
###############################################################################

plot.12 = function(data,
                   x1 = "Time",
                   y1 = "Mean.cell.intensity",
                   x2 = "Time",
                   y2 = "dist2top", batch = 1) {
  par(mfrow = c(4, 3), mar = c(2, 2, 2, 2))
  for (ii in ((batch-1)*12 + 1):((batch)*12)) {
    plot(data[[ii]][, x1],
         data[[ii]][, y1],
         ylim = c(0, 160),
         xlim = c(0, 84))
    par(new = TRUE)
    plot(
      data[[ii]][, x2],
      data[[ii]][, y2],
      col = 2,
      ylim = c(0, 10),
      xaxt = "n",
      yaxt = "n",
      xlim = c(0, 84)
    )
    axis(4)
  }
}
