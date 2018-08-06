setwd("/home/henrik/compbio/thesis/code/newTreat/")
library(tidyverse)
# source("functions.R")
# source("read_data.R")
source("aux.R")
library(gridExtra)
plants = c(1, 2, 4, 13, 15) #18

load("allbut18.RData")

cld = lapply(results, "[[", "lineage.sublines.data")

qd = do.call(rbind, lapply(results, "[[", "quant.data") %>% unlist(recursive = F))
dd = do.call(rbind, lapply(results, "[[", "division.data"))
sd = do.call(cbind, lapply(results, "[[", "segm.data"))
topcells = qd %>% filter(d2t==0) %>% select(plant, t, n.x, n.y, n.z)

# sd = sd[c("plant", "id", "m.vol", "neighs", "layer"), ]
# dd = dd[, c("t", "age", "id", "plant", "n.vol", "expr", "d2t", "m.vol", "layer")]
qd = qd[, -(c(4:6))]

sd.n = colnames(sd)
sd = lapply(1:ncol(sd), function(x) as.tibble(sd[, x])) 
sd = lapply(1:length(sd), function(x) {
  add_column(sd[[x]], t = rep(as.integer(sd.n[x]), nrow(sd[[x]])), .before = 1)
}) %>% do.call(rbind, .)

d2ts = sd %>% merge(topcells, by=c("plant", "t")) %>% mutate(d2t = sqrt((m.x - n.x)**2 + (m.y - n.y)**2 + (m.z - n.z)**2)) %>% select(t, plant, id, d2t) %>% as.tibble
sd = sd %>% merge(d2ts, by=c("t","plant","id")) %>% as.tibble

# m.d2ts = dd %>% merge(sd, by=c("plant", "t")) %>% merge(topcells, by=c("plant", "t")) %>% mutate(d2t = sqrt((m.x - n.x)**2 + (m.y - n.y)**2 + (m.z - n.z)**2)) %>% select(t, plant, id, d2t) %>% as.tibble
# dd = dd %>% merge(m.d2ts, by=c("t","plant","id")) %>% as.tibble

############################
### Plots
############################
dd %>% ggplot(aes(x = m.vol)) + geom_histogram(bins = 200)
dd %>% ggplot(aes(x = m.vol, y = child1.m.vol, colour = d2t)) + geom_point() + xlim(0, 1000) + ylim(0, 1000)


dd %>% ggplot(aes(x = age)) + geom_histogram(bins = 50) + facet_wrap(~layer, ncol = 1) + labs(x = "Age", y = "Frequency")
dd %>% ggplot(aes(x = n.vol)) + geom_histogram(bins = 101)
dd %>% filter(plant != 1) %>% ggplot(aes(x = expr)) + geom_histogram(bins = 30) + facet_wrap(~plant)

p1 = qd %>% ggplot(aes(x = expr)) + geom_histogram(bins = 30) + facet_wrap(~plant, nrow=1) + labs(title="All data", x = " ")
p2 = dd %>% filter(plant!=1) %>% ggplot(aes(x = expr)) + geom_histogram(bins = 30) + facet_wrap(~plant, nrow=1) + labs(title="Division data", x = "Expression")
grid.arrange(p1,p2, nrow = 2)




### Histograms: Age v. layer
# dd %>% ggplot(aes(x = age)) + geom_histogram(bins = 50) + facet_wrap(~ layer)
# dd %>% ggplot(aes(x = age)) + geom_histogram(bins = 100) + facet_wrap(~ plant)

### Layer: Volume v. gene expression values
p1 = ggplot(qd %>% group_by(plant, t) %>% summarise(n = n()),
            aes(x = t, y = n)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ plant) +
  labs(title = "Nuclei per time")

p2 = ggplot(sd %>% filter(plant != 1) %>% group_by(t, plant) %>% summarize(n = n()),
            aes(x = t, y = n)) +
  geom_line() + geom_point() +
  facet_wrap( ~ plant) +
  labs(title = "Cells per time")

p3 = ggplot(dd %>% filter(plant != 1) %>% group_by(t, plant) %>% summarize(n = n()),
            aes(x = t, y = n)) +
  geom_line() + geom_point() +
  facet_wrap( ~ plant) +
  labs(title = "Divisions per time")

grid.arrange(p2, p3, ncol = 1)
grid.arrange(p1, p2, ncol = 1)

qd %>% ggplot(aes(x = n.vol, y = expr, colour = d2t)) + geom_point() + facet_wrap( ~ layer) 
qd %>% ggplot(aes(x = n.vol, y = expr, colour = layer)) + geom_point() + labs(x = "Nuclear volume", y = "Expression")

# # nDivs correlated w/ nNuclei
# cor.test(
#   qd %>% filter(plant != 1, t != 84) %>%
#     group_by(t, plant) %>% summarize(n = n()) %>% .$n,
#   dd %>% filter(plant != 1) %>%
#     group_by(t, plant) %>% summarize(n = n()) %>% .$n
# )
# 
#   group_by(t, plant) %>% summarize(n = n()) %>% .$n

# nNuclei correlated w/ nCells
# cor.test(
#   qd %>% filter(plant != 1) %>%
#     group_by(t, plant) %>% summarize(n = n()) %>% .$n,
#   sd %>% filter(plant != 1, t != 84 && (t!=76 & !plant %in% c(2,4))) %>%
#     group_by(t, plant) %>% summarize(n = n()) %>% .$n
# )
# 
# nDivs weakly correlated with number of nCells
# cor.test(
#   sd %>% filter(plant != 1, t != 84) %>%
#     group_by(t, plant) %>% summarize(n = n()) %>% .$n,
#   dd %>% filter(plant != 1, id %in% sd$id) %>%
#     group_by(t, plant) %>% summarize(n = n()) %>% .$n
# )

p1 = qd %>% ggplot(aes(x = d2t, y = expr, colour = layer)) + geom_point() + labs(title = "All")
p2 = qd %>% filter(layer == "L1") %>% ggplot(aes(x = d2t, y = expr, colour = t)) + stat_binhex() + labs(title = "L1")
grid.arrange(p1, p2, ncol = 2)

