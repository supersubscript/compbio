#!/usr/bin/env Rscript
source("aux.R")
source("read_data.R")
source("process_data.R")

plants = c(2)
plant.no = 15
data = read.data(plant.no, map.quant = TRUE)

mapping.data = data$mapping.data
quant.data   = data$quant.data
timepoints   = data$timepoints

# quant.data = lapply(quant.data, function(x) cbind(x, dist2top = rep(NA, nrow(x))))
quant.data[[1]] = quant.data[[1]][quant.data[[1]][, 1] %in% get.mom(mapping.data[[1]]), ]

for (ii in 2:length(mapping.data)){
  quant.data[[ii]] = quant.data[[ii]][
    quant.data[[ii]][,1] %in% union(get.mom(mapping.data[[ii]]), get.dau(mapping.data[[ii-1]])), ]
}

quant.data[[length(quant.data)]] =
  quant.data[[length(quant.data)]][quant.data[[length(quant.data)]][, 1] %in% 
                                     get.dau(mapping.data[[length(mapping.data)]]), ]


# return(list(timepoints=timepoints, quant.data=quant.data))
# }
# )

#############################source("process_data.R")

# timepoints = lapply(results, "[[", "timepoints")
# all.data =   lapply(results, "[[", "quant.data")
all.data = quant.data  
timepoints = unlist(timepoints)

no.clv3.nuclei = unlist(lapply(all.data, nrow))
no.clv3.nuclei.within = unlist(lapply(all.data, function(x) nrow(x[which(x$dist2top < 5), ])))

par(mfrow = c(2, 1))
plot(timepoints[-7], no.clv3.nuclei, type = "b")
plot(timepoints[-7], no.clv3.nuclei.within, type = "b")

par(mfrow=c(5,2))
for(ii in 1:10){
  no.clv3.nuclei.within = unlist(lapply(all.data, function(x) nrow(x[which(x$dist2top < ii), ])))
  plot(timepoints[-7], no.clv3.nuclei.within/no.clv3.nuclei, type = "b", main = ii, ylim=0:1)
}

############################################
top.ones = lapply(all.data, function(x) x[order(x$dist2top), ][1:5, ])
# all.data = lapply(all.data, function(x) x[which(x$dist2top < 2), ])
all.but.top = lapply(all.data, function(x) x[order(x$dist2top), ][-c(1:5), ])

plot.20 = function(plt.fct, x, y, xlim = NULL, ylim = NULL) {
  par(mfrow = c(5, 4))
  if(!missing(y))
  for(ii in 1:20) {
    rbPal = colorRampPalette(1:2)
    ccc = rbPal(10)[as.numeric(cut(all.data[[ii]][, "dist2top"], breaks = 10))]
    plt.fct(all.data[[ii]][, x], all.data[[ii]][, y], type = "p", main = ii, col = ccc, xlim = xlim, ylim = ylim)
  }
  else
    for(ii in 1:20) {
      rbPal = colorRampPalette(1:2)
      ccc = rbPal(10)[as.numeric(cut(all.data[[ii]][, "dist2top"], breaks = 10))]
      plt.fct(all.data[[ii]][, x], type = "p", main = ii, col = ccc, xlim = xlim, ylim = ylim)
    } 
}
plot.20(plot,"Mean.cell.intensity", "Boa.volume", xlim = c(0,165), ylim = 0:1)
plot.20(plot,"Mean.cell.intensity", "dist2top", xlim = c(0,165), ylim=c(0,7))
plot.20(plot,"dist2top", "Boa.volume", xlim = c(0,7), ylim = 0:1)


par(mfrow=c(5,4))
for(ii in 1:20){
  a = hist(all.but.top[[ii]]$Boa.volume, breaks = 30)
}
par(mfrow=c(5,4))
for(ii in 1:20)
  hist(all.but.top[[ii]]$Mean.cell.intensity, breaks = 30, xlim = c(0,160))

par(mfrow=c(5,4))
for(ii in 1:20)
  hist(all.but.top[[ii]]$dist2top, breaks = 30)

