

plants = c(2, 4, 13, 15, 18)
no.plants = length(plants)

get.coordinates = function(plant.no){
  quant.files = list.files(paste0("../data/msb_plants/plant", plant.no, "/results"), full.names = TRUE)
  no.quant.files = length(quant.files)
  
  # Put in right order
  order.quant = sapply(quant.files, function(x) strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])]) 
  order.quant = unname(sapply(order.quant, function(x) as.integer(gsub("t|\\.txt", "", x))))
  quant.files = quant.files[order(order.quant)]
  
  # Read in quantified data
  quant.data = lapply(quant.files, function(x) read.table(x, header=T, sep = "\t"))
  return(lapply(quant.data, "[", 2:4))
}

coordinate.data = lapply(plants, get.coordinates)
max.time = max(sapply(coordinate.data, length))
max.cells = max(unlist(sapply(coordinate.data, function(x) sapply(x, nrow))))

xs = lapply(coordinate.data, function(x) lapply(x, function(y) c(y[,1], rep(NA, max.cells - length(y[,1])))))
ys = lapply(coordinate.data, function(x) lapply(x, function(y) c(y[,2], rep(NA, max.cells - length(y[,2])))))
zs = lapply(coordinate.data, function(x) lapply(x, function(y) c(y[,3], rep(NA, max.cells - length(y[,3])))))

xs = lapply(xs, function(x) c(x, rep(list(rep(NA, max.cells)), max.time - length(x))))
ys = lapply(ys, function(x) c(x, rep(list(rep(NA, max.cells)), max.time - length(x))))
zs = lapply(zs, function(x) c(x, rep(list(rep(NA, max.cells)), max.time - length(x))))

xs = t(xs)
ys = t(ys)
zs = t(zs)

# all.xs = t(do.call(cbind, lapply(1:no.plants, function(x) do.call(rbind, xs[,x][[1]]))))
# all.ys = t(do.call(cbind, lapply(1:no.plants, function(x) do.call(rbind, ys[,x][[1]]))))
# all.zs = t(do.call(cbind, lapply(1:no.plants, function(x) do.call(rbind, zs[,x][[1]]))))
# 
# boxplot(all.xs)
# boxplot(all.ys)
# boxplot(all.zs)
# par(mfrow = c(1, 1))
# for (plant in plants) {
#   plot(NA, xlim = c(0, 84), ylim = c(0, 14))
#   for (ii in 1:max.time) {
#     points(rep((ii - 1) * 4, length(all.zs[, ii])), all.zs[, ii])
#   }
# }
# 
# new = normalize.loess(all.zs[,-c(21:22)])
# boxplot(new)
# 
# 
# # Put all on same mean
# means = apply(all.zs, 2, mean, na.rm = T)
# means.mean = mean(means)
# for(ii in 1:nrow(all.zs))
#   all.zs[ii, ] = all.zs[ii, ] - means 
# all.zs = all.zs + means.mean
# 
# plotDensity(all.zs)
# 
# par(mfrow = c(1, 1))
# for (plant in plants) {
#   plot(NA, xlim = c(0, 84), ylim = c(0, 14))
#   for (ii in 1:max.time) {
#     points(rep((ii - 1) * 4, length(all.zs[, ii])), all.zs[, ii])
#   }
# }



# df = data.frame(means=means, x = 1:22)
# fit = lm(means~x, data = df)
# out = predict(fit, df$x)
# predict(fit, interval = "prediction")


# no.plants = 5
# # Normalise every timepoint
# for(ii in 1:max.time){
#   xs.norm = scale(do.call(cbind, lapply(xs[1,], "[[", ii)))  
#   ys.norm = scale(do.call(cbind, lapply(ys[1,], "[[", ii))) 
#   zs.norm = scale(do.call(cbind, lapply(zs[1,], "[[", ii))) 
#   # xs.norm = normalize.quantiles(do.call(cbind, lapply(xs[1,], "[[", ii)))
#   # ys.norm = normalize.quantiles(do.call(cbind, lapply(ys[1,], "[[", ii)))
#   # zs.norm = normalize.quantiles(do.call(cbind, lapply(zs[1,], "[[", ii)))
#   for(jj in 1:no.plants){
#     xs[[jj]][[ii]] = xs.norm[,jj]
#     ys[[jj]][[ii]] = ys.norm[,jj]
#     zs[[jj]][[ii]] = zs.norm[,jj]
#   }
# }
  colnames(xs) = plants
colnames(ys) = plants
colnames(zs) = plants

# par(mfrow=c(1,1))
# for(ii in 1:max.time){
#   dist = do.call(cbind, lapply(zs[1,], "[[", ii))
#   plotDensity(dist)
#   Sys.sleep(.5)
# }


