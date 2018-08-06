setwd("/home/henrik/compbio/thesis/code/newTreat/")
library(plyr)
library(tidyverse)
source("functions.R")
source("read_data.R")
source("aux.R")
source("../plot.R")
plants = c(1, 2, 4, 13, 15, 18)

# results = lapply(plants, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 1)
# save(file="data/allbut18corresp_d2t_peakexpr_1.RData", results)
# results = lapply(plants, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 2)
# save(file="data/allbut18corresp_d2t_peakexpr_2.RData", results)
# results = lapply(plants, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 3)
# save(file="data/allbut18corresp_d2t_peakexpr_3.RData", results)
# results = lapply(plants, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 4)
# save(file="data/allbut18corresp_d2t_peakexpr_4.RData", results)
# results = lapply(plants, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy",    d2t.n = 1)
# save(file="data/allbut18corresp_d2t_topxy_1.RData", results)
# results = lapply(plants, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy",    d2t.n = 2)
# save(file="data/allbut18corresp_d2t_topxy_2.RData", results)
# results = lapply(plants, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy",    d2t.n = 3)
# save(file="data/allbut18corresp_d2t_topxy_3.RData", results)
# results = lapply(plants, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "top.xy",    d2t.n = 4)
# save(file="data/allbut18corresp_d2t_topxy_4.RData", results)

# I'm an idiot.
load("data/allbut18corresp_d2t_peakexpr_1.RData")
sd.e1 = do.call(bind_rows, lapply(results, "[[", "segm.data"))  %>% ungroup() %>% select(plant, t, m.x, m.y, m.z, d2t, circ)
load("data/allbut18corresp_d2t_peakexpr_2.RData") 
sd.e2 = do.call(bind_rows, lapply(results, "[[", "segm.data"))  %>% ungroup() %>% select(plant, t, m.x, m.y, m.z, d2t, circ)
load("data/allbut18corresp_d2t_peakexpr_3.RData") 
sd.e3 = do.call(bind_rows, lapply(results, "[[", "segm.data"))  %>% ungroup() %>% select(plant, t, m.x, m.y, m.z, d2t, circ)
load("data/allbut18corresp_d2t_peakexpr_4.RData") 
sd.e4 = do.call(bind_rows, lapply(results, "[[", "segm.data"))  %>% ungroup() %>% select(plant, t, m.x, m.y, m.z, d2t, circ)
load("data/allbut18corresp_d2t_topxy_1.RData")
sd.xy1 = do.call(bind_rows, lapply(results, "[[", "segm.data")) %>% ungroup() %>% select(plant, t, m.x, m.y, m.z, d2t, circ)
load("data/allbut18corresp_d2t_topxy_2.RData")
sd.xy2 = do.call(bind_rows, lapply(results, "[[", "segm.data")) %>% ungroup() %>% select(plant, t, m.x, m.y, m.z, d2t, circ)
load("data/allbut18corresp_d2t_topxy_3.RData")
sd.xy3 = do.call(bind_rows, lapply(results, "[[", "segm.data")) %>% ungroup() %>% select(plant, t, m.x, m.y, m.z, d2t, circ)
load("data/allbut18corresp_d2t_topxy_4.RData")
sd.xy4 = do.call(bind_rows, lapply(results, "[[", "segm.data")) %>% ungroup() %>% select(plant, t, m.x, m.y, m.z, d2t, circ)

#######################################
fcts = c("e1", "e2", "e3", "e4", "xy1", "xy2", "xy3", "xy4", "para")
sd.e1  = sd.e1  %>% mutate(var = factor("e1",  levels = c("e1", "e2", "e3", "e4", "xy1", "xy2", "xy3", "xy4", "para")))
sd.e2  = sd.e2  %>% mutate(var = factor("e2",  levels = c("e1", "e2", "e3", "e4", "xy1", "xy2", "xy3", "xy4", "para")))
sd.e3  = sd.e3  %>% mutate(var = factor("e3",  levels = c("e1", "e2", "e3", "e4", "xy1", "xy2", "xy3", "xy4", "para")))
sd.e4  = sd.e4  %>% mutate(var = factor("e4",  levels = c("e1", "e2", "e3", "e4", "xy1", "xy2", "xy3", "xy4", "para")))
sd.xy1 = sd.xy1 %>% mutate(var = factor("xy1", levels = c("e1", "e2", "e3", "e4", "xy1", "xy2", "xy3", "xy4", "para")))
sd.xy2 = sd.xy2 %>% mutate(var = factor("xy2", levels = c("e1", "e2", "e3", "e4", "xy1", "xy2", "xy3", "xy4", "para")))
sd.xy3 = sd.xy3 %>% mutate(var = factor("xy3", levels = c("e1", "e2", "e3", "e4", "xy1", "xy2", "xy3", "xy4", "para")))
sd.xy4 = sd.xy4 %>% mutate(var = factor("xy4", levels = c("e1", "e2", "e3", "e4", "xy1", "xy2", "xy3", "xy4", "para")))

data = bind_rows(sd.e1, sd.e2, sd.e3, sd.e4, sd.xy1, sd.xy2, sd.xy3, sd.xy4)
# data = data %>% 
  # group_by(plant, t, var) %>% 
  # mutate(m.x = m.x / max(m.x), m.y = m.y / max(m.y), m.z = m.z / max(m.z)) %>% ungroup() # Make plant a factor
data = data %>% mutate(plant = factor(plant))

top.coords = data %>% 
  filter(!is.na(d2t), circ == 0) %>% 
  group_by(plant, t, var)        %>% 
  summarize(x = mean(m.x),
            y = mean(m.y),
            z = mean(m.z)) %>% ungroup()

#######################################
### Parabloid top
#######################################
files = sapply(list.files("/home/henrik/plant_data", full.names = TRUE),
               list.files,
               full.names = TRUE) %>%
  unlist %>%
  unname
para.data = lapply(files, read.csv, sep = ",", header = TRUE)

# Tidy the para.data: Add plant and time info
time  = extract.numbers(files)[1, ]
plant = extract.numbers(files)[2, ]
para.data = lapply(1:length(para.data), function(x) 
  cbind(plant = rep(plant[x], nrow(para.data[[x]])),
        t     = rep(time[x], nrow(para.data[[x]])),
        para.data[[x]]))
para.data = do.call(rbind, para.data) %>% as.tibble()

para.data = para.data %>% filter(para_apex_x > 0, para_apex_y > 0, para_apex_z > 0) %>%
  arrange(plant, t) %>%
  select(plant, t, para_apex_x, para_apex_y, para_apex_z) %>%
  dplyr::rename(z = para_apex_y, y = para_apex_x, x = para_apex_z) %>% # Not sure about xyz mappings
  mutate(var = factor("para", levels = fcts)) %>% mutate(plant = factor(plant)) %>%
  select(plant, t, var, x, y, z) 

top.coords = bind_rows(top.coords, para.data) %>% arrange(plant, t, var)
top.coords[which(top.coords$var == "para"), c(4:6)] = 0.26 * top.coords[which(top.coords$var == "para"), c(4:6)]  # Set correct spacing

#######################################
### Parabloid top
#######################################

# top.coords[which(top.coords$var == "para"), 4:5] = top.coords[which(top.coords$var == "para"), 5:4]

ggplot(top.coords %>% filter(!plant %in% c(1,18), var == "e1"),
       aes(x = fct_inx, y = y, colour = var)) + 
  geom_point(data  = top.coords %>% filter(!plant %in% c(1,18),  var == "e4"), aes(x = x, y = y, colour = var)) +
  geom_path(data   = top.coords %>% filter(!plant %in% c(1,18),  var == "e4"),
            aes(x  = top.coords %>% filter(!plant %in% c(1,18),  var == "e4") %>% .$x, 
                y  = top.coords %>% filter(!plant %in% c(1,18),  var == "e4") %>% .$y, colour = var))  +
  # geom_point(data  = top.coords %>% filter(!plant %in% c(1,18),  var == "para"),
  #            aes(x = top.coords %>% filter(!plant %in% c(1,18),  var == "para") %>% .$x,
  #                y = top.coords %>% filter(!plant %in% c(1,18),  var == "para") %>% .$y, colour = var))  +
  # geom_path(data   = top.coords %>% filter(!plant %in% c(1,18),  var == "para"),
  #           aes(x  = top.coords %>% filter(!plant %in% c(1,18),  var == "para") %>% .$x,
  #               y  = top.coords %>% filter(!plant %in% c(1,18),  var == "para") %>% .$y, colour = var))  +
  geom_point(data  = top.coords %>% filter(!plant %in% c(1,18),  var == "e2"),
             aes(x  = top.coords %>% filter(!plant %in% c(1,18), var == "e2") %>% .$x,
                 y  = top.coords %>% filter(!plant %in% c(1,18), var == "e2") %>% .$y, colour = var))  +
  geom_path(data   = top.coords %>% filter(!plant %in% c(1,18),  var == "e2"),
            aes(x  = top.coords %>% filter(!plant %in% c(1,18),  var == "e2") %>% .$x,
                y  = top.coords %>% filter(!plant %in% c(1,18),  var == "e2") %>% .$y, colour = var))  +
  geom_point(data  = top.coords %>% filter(!plant %in% c(1,18),  var == "e1"),
            aes(x  = top.coords %>% filter(!plant %in% c(1,18),  var == "e1") %>% .$x,
                y  = top.coords %>% filter(!plant %in% c(1,18),  var == "e1") %>% .$y, colour = var))  +
  geom_path(data   = top.coords %>% filter(!plant %in% c(1,18),  var == "e1"),
            aes(x  = top.coords %>% filter(!plant %in% c(1,18),  var == "e1") %>% .$x,
                y  = top.coords %>% filter(!plant %in% c(1,18),  var == "e1") %>% .$y, colour = var))  +
  facet_wrap( ~ plant, scales = "free") + 
  scale_colour_discrete(name = "Measure",
                        labels = c("Top1-expr", "Paraboloid", "Top4-z"))


# top.coords %>% group_by(plant, t, var) %>% mutate(x = x / max(x), y = y / max(y), z = z / max(z))
p.x  = (top.coords %>% filter(plant == 2,  var == "para") %>% .$x)
e1.x = (top.coords %>% filter(plant == 2,  var == "e1") %>% .$x)
ggplot(data.frame(cbind(p.x, e1.x)), aes(x=p.x, y = e1.x)) + geom_path(size = 1)

p.y  = (top.coords %>% filter(plant == 15,  var == "para") %>% .$y)
e1.y = (top.coords %>% filter(plant == 15,  var == "e1") %>% .$y)
ggplot(data.frame(cbind(p.y, e1.y)), aes(x=p.y, y = e1.y)) + geom_path(size = 1)

