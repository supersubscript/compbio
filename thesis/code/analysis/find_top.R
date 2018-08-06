#!/usr/bin/env Rscript
setwd("/home/henrik/compbio/thesis/code/newTreat/")
library(plyr)
library(tidyverse)
source("functions.R")
source("read_data.R")
source("aux.R")
source("../plot.R")

pp = c(1, 2, 4, 13, 15, 18)
p = c(2, 4, 13, 15)
# re1   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 1)
# re2   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 2)
# re4   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 4)
# re8   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 8)
# re16  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 16)
# re32  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 32)
# re64  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 64)
# reAll = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = "all")
# 
# rs1   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy", d2t.n = 1)
# rs2   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy", d2t.n = 2)
# rs4   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy", d2t.n = 4)
# rs8   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy", d2t.n = 8)
# rs16  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy", d2t.n = 16)
# rs32  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy", d2t.n = 32)
# rs64  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy", d2t.n = 64)
# rsAll = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy", d2t.n = "all")
# 
# rw1   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "w.expr", d2t.n = 1)
# rw2   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "w.expr", d2t.n = 2)
# rw4   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "w.expr", d2t.n = 4)
# rw8   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "w.expr", d2t.n = 8)
# rw16  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "w.expr", d2t.n = 16)
# rw32  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "w.expr", d2t.n = 32)
# rw64  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "w.expr", d2t.n = 64)
# rwAll = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "w.expr", d2t.n = "all")

# results = list(re1, re2, re4, re8, re16, re32, re64, reAll, 
               # rs1, rs2, rs4, rs8, rs16, rs32, rs64, rsAll, 
               # rw1, rw2, rw4, rw8, rw16, rw32, rw64, rwAll)
# save(file = "data/d2t_comparisons2.RData", results)
# load(file = "data/d2t_comparisons2.RData")

# qd  = do.call(rbind, lapply(results[[19]], "[[", "quant.data"))
# sd  = do.call(rbind, lapply(results[[19]], "[[", "segm.data"))
# dd  = do.call(rbind, lapply(results[[19]], "[[", "division.data"))
# cld = lapply(results[[19]], "[[", "lineage.sublines.data")
# save(qd, sd, dd, cld, file = "data/filtered_weightexpr4.RData")

# 


load(file="data/d2t_comparisons2.RData")
results = do.call(cbind, results)
results = results[-c(1,6), ]
rownames(results) = p
colnames(results) = c(rep(c(1, 2, 4, 8, 16, 32, 64, 128), 3))
# 
# # results = apply(results, 1:2, function(x) x["quant.data"])
results = apply(results, 1, function(x){
  x = lapply(x, "[[", "quant.data")
  ns  = colnames(results)
  types = c(c(rep("e", 8), rep("s", 8), rep("w", 8)))

  x = lapply(1:length(x), function(y) cbind(type = types[y], n = rep(ns[y], nrow(x[[y]])), x[[y]]) %>% as.tibble())
  x = do.call(rbind, x) %>% as.tibble()
  #
  data = x %>%
    filter(circ == 0) %>%
    select(plant, type, t, n, n.x, n.y, n.z, d2t, circ, layer) %>%
    group_by(plant, type, n, t) %>%
    summarize(x = mean(n.x), y = mean(n.y), z = max(n.z)) %>% arrange(type, n, t)
  data$n = as.integer(as.character(data$n))
  data = arrange(data, type, n)
  data$n = as.factor(data$n)
  data %>% ungroup()
})
results = do.call(bind_rows, results)
# write.table(results, file = "topcoords_tomax.dat", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# 
# # results = #lapply(results, "[[", "quant.data")
# # ns  = colnames(results)
# # types = c(c(rep("e", 6), rep("s", 6), rep("w", 6)), c("e", "e", "w", "w", "s", "s"))
# # 
# # 
# # results = lapply(1:length(results), function(x) cbind(type = types[x], n = rep(ns[x], nrow(results[[x]])), results[[x]]) %>% as.tibble())
# # results = do.call(rbind, results) %>% as.tibble()
# # # 
# data = results %>%
#   filter(circ == 0) %>%
#   select(type, t, n, n.x, n.y, n.z, d2t, circ, layer) %>%
#   group_by(type, n, t) %>%
#   summarize(x = mean(n.x), y = mean(n.y), z = max(n.z)) %>% arrange(type, n, t)
# data$n = as.integer(as.character(data$n))
# data = arrange(data, type, n)
# data$n = as.factor(data$n)
# data %>% ungroup()
# 
# 
data = results
p1 = ggplot(data %>% filter(plant == 2, t < 20), aes(x = n, y = x, colour = type, group = type)) +
  geom_point() +
  geom_line(aes(x = as.integer(n), y = x)) +
  facet_grid(~t, scales = "free") +
  theme_bw()
p2 = ggplot(data %>% filter(plant == 2, t < 20), aes(x = n, y = y, colour = type, group = type)) +
  geom_point() +
  geom_line(aes(x = as.integer(n), y = y)) +
  facet_grid(~t, scales = "free") +
  theme_bw()
grid.arrange(p1,p2, ncol = 1)











# 
# #!/usr/bin/env Rscript
# setwd("/home/henrik/compbio/thesis/code/newTreat/")
# library(plyr)
# library(tidyverse)
# source("functions.R")
# source("read_data.R")
# source("aux.R")
# source("../plot.R")
# 
# pp = c(1, 2, 4, 13, 15, 18)
# p = c(2, 4, 13, 15)
# re1   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 1)
# re2   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 2)
# re4   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 4)
# re8   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 8)
# re16  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 16)
# reAll = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = "all")
# 
# rs1   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy", d2t.n = 1)
# rs2   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy", d2t.n = 2)
# rs4   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy", d2t.n = 4)
# rs8   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy", d2t.n = 8)
# rs16  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy", d2t.n = 16)
# rsAll = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy", d2t.n = "all")
# 
# rw1   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "w.expr", d2t.n = 1)
# rw2   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "w.expr", d2t.n = 2)
# rw4   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "w.expr", d2t.n = 4)
# rw8   = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "w.expr", d2t.n = 8)
# rw16  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "w.expr", d2t.n = 16)
# rwAll = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "w.expr", d2t.n = "all")


# re32  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 32)
# rs32  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy",    d2t.n = 32)
# rw32  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "w.expr",    d2t.n = 32)
# re64  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 64)
# rs64  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy",    d2t.n = 64)
# rw64  = lapply(pp, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "w.expr",    d2t.n = 64)

# results = c(results, list(re32, re64, rw32, rw64, rs32, rs64))
# # results = list(re1, re2, re4, re8, re16, reAll, rs1, rs2, rs4 , rs8, rs16, rsAll, rw1, rw2, rw4, rw8, rw16, rwAll)
# save(file="data/d2t_comparisons.RData", results)
# load(file="data/d2t_comparisons.RData")
# results = do.call(cbind, results)
# results = results[-c(1,6), ]
# rownames(results) = p
# colnames(results) = c(rep(c(1, 2, 4, 8, 16, 128), 3), c(32, 64, 32, 64, 32, 64))
# 
# # results = apply(results, 1:2, function(x) x["quant.data"])
# apply(results, 1, function(x){
#   ns  = colnames(results)
#   types = c(c(rep("e", 6), rep("s", 6), rep("w", 6)), c("e", "e", "w", "w", "s", "s"))
#   
#   
#   results = lapply(1:length(results), function(x) cbind(type = types[x], n = rep(ns[x], nrow(results[[x]])), results[[x]]) %>% as.tibble())
#   results = do.call(rbind, results) %>% as.tibble()
#   # 
#   data = results %>%
#     filter(circ == 0) %>%
#     select(type, t, n, n.x, n.y, n.z, d2t, circ, layer) %>%
#     group_by(type, n, t) %>%
#     summarize(x = mean(n.x), y = mean(n.y), z = max(n.z)) %>% arrange(type, n, t)
#   data$n = as.integer(as.character(data$n))
#   data = arrange(data, type, n) 
#   data$n = as.factor(data$n)
#   data %>% ungroup()
#   
# })
# # results = #lapply(results, "[[", "quant.data")
# ns  = colnames(results)
# types = c(c(rep("e", 6), rep("s", 6), rep("w", 6)), c("e", "e", "w", "w", "s", "s"))
# 
# 
# results = lapply(1:length(results), function(x) cbind(type = types[x], n = rep(ns[x], nrow(results[[x]])), results[[x]]) %>% as.tibble())
# results = do.call(rbind, results) %>% as.tibble()
# # 
# data = results %>%
#   filter(circ == 0) %>%
#   select(type, t, n, n.x, n.y, n.z, d2t, circ, layer) %>%
#   group_by(type, n, t) %>%
#   summarize(x = mean(n.x), y = mean(n.y), z = max(n.z)) %>% arrange(type, n, t)
# data$n = as.integer(as.character(data$n))
# data = arrange(data, type, n) 
# data$n = as.factor(data$n)
# data %>% ungroup()
# 
# 
# p1 = ggplot(data %>% filter(t < 20), aes(x = n, y = x, colour = type, group = type)) +
#   geom_point() +
#   geom_line(aes(x = as.integer(n), y = x)) +
#   facet_grid(~t, scales = "free") + 
#   theme_bw()
# p2 = ggplot(data %>% filter(t < 20), aes(x = n, y = y, colour = type, group = type)) +
#   geom_point() +
#   geom_line(aes(x = as.integer(n), y = y)) +
#   facet_grid(~t, scales = "free") + 
#   theme_bw()
# grid.arrange(p1,p2, ncol = 1)

