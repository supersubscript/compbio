#Plot.R

###############################################################################
### CUSTOM PLOT FUNCTIONS
###############################################################################
plot.lineage = function(lineage.sublines.data,
                        var1 = "expr",
                        var2 = "d2t",
                        ylim = c(0, 160),
                        both = T) {
  plot(
    NA,
    xlim = c(0, 84),
    ylim = c(0, 160),
    xlab = "Time",
    ylab = "Expression",
    xaxt = "n",
    yaxt = "n"
  )
  for (ii in 1:length(lineage.sublines.data)) {
    lines(
      lineage.sublines.data[[ii]][, "t"]  %>% unlist(),
      lineage.sublines.data[[ii]][, var1] %>% unlist(),
      lwd = lw.s,
      col = length(lineage.sublines.data) + 1 - ii
    )
  }
  
  if (both) {
    plot(
      NA,
      xlim = c(0, 84),
      ylim = c(0, 7),
      xlab = "Time",
      ylab = "Distance to apex",
      xaxt = "n",
      yaxt = "n"
    )
    for (ii in 1:length(lineage.sublines.data))
      lines(
        lineage.sublines.data[[ii]][, "t"]  %>% unlist(),
        lineage.sublines.data[[ii]][, var2] %>% unlist(),
        lwd = lw.s,
        col = length(lineage.sublines.data) + 1 - ii
      )
  }
  print(lineage.sublines.data[[ii]]$layer, max.levels = 0)
}

###############################################################################
### HELPER PLOT FUNCTIONS
###############################################################################

plot.12 = function(data, x1 = "Time", y1 = "Mean.cell.intensity", x2 = "Time", y2 = "dist2top", batch = 1) {
  par(mfrow = c(4, 3), mar = c(2, 2, 2, 2))
  for (ii in ((batch-1)*12 + 1):((batch)*12)) {
    plot(data[[ii]][, x1], data[[ii]][, y1], ylim = c(100, 160), xlim = c(0, 84), main = names(data[[ii]]))
    par(new = TRUE)
    plot(data[[ii]][, x2], data[[ii]][, y2], col = 2, ylim = c(0, 10), 
         xaxt = "n", yaxt = "n", xlim = c(0, 84), main = names(data[[ii]]))
    axis(4)
  }
}

plot.stat.hists = function(all.data, no.breaks, density = FALSE){
  par(mfrow = c(4, 3), mar=c(2,2,2,2))
  if(density == FALSE)
    for(ii in 6:ncol(all.data)) if(!all(is.na(all.data[,ii]))) hist(all.data[, ii], breaks = no.breaks, main = colnames(all.data)[ii])
  # else
  # for(ii in 6:ncol(all.data)) if(!all(is.na(all.data[,ii]))) plotDensity(data.frame(x=all.data[, ii]), main = colnames(all.data)[ii])
}

plot.all.stats = function(all.data){
  graphics.off(); plot.new()
  par(mfrow = c(3, 3), mar = c(2, 2, 2, 2))
  plots       = vector(mode = "list", length = (ncol(all.data) - 5))
  interesting = (1:ncol(all.data))[-c(1:5, 8:10)] # Remove xyz and indices
  counter     = 1
  for (ii in interesting) {
    not.ii = interesting[interesting != ii]
    for (jj in not.ii) {
      if (length(all.data[!is.na(all.data[, ii]), ii]) != 0 &
          length(all.data[!is.na(all.data[, jj]), jj]) != 0) {
        x.lab       = colnames(all.data)[jj]
        y.lab       = colnames(all.data)[ii]
        is.not.na   = !is.na(all.data[, y.lab]) & !is.na(all.data[, x.lab])
        correlation = cor(all.data[is.not.na, x.lab] , all.data[is.not.na, y.lab])
        title       = paste0(y.lab, " vs ", x.lab, "\n cor: ", round(correlation, 2))
        plot(all.data[, jj], all.data[, ii], main = title, xlab = x.lab, ylab = y.lab)
      }
    }
    plots[[counter]] = recordPlot()
    plot.new()
    counter = counter + 1
  }
  return(plots)
}

plot.violin = function(data, var1, var2){
  t = data[,var1]; expr = data[, var2]; there = nna(t) & nna(expr)
  p = ggplot(data.frame(t=t[there], expr=expr[there]), aes(factor(t), expr)); p1 = p + geom_violin() + labs(y=var2, x=var1)
  return(p1)
}