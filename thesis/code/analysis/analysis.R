setwd("/home/henrik/compbio/thesis/code/newTreat/")
library(plyr)
library(tidyverse)
source("functions.R")
source("read_data.R")
source("aux.R")
source("../plot.R")
plants = c(1, 2, 4, 13, 15, 18)

# results = lapply(plants, get.plant.data, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 4)

# save(file="data/allbut18corresp.RData_d2t_peakexpr_4_expr_filtering.RData", results)
# results = re4
# save(file="data/allbut18corresp_d2t_peakexpr_4_expr_filtering.RData", results)
load(file="data/allbut18corresp_d2t_peakexpr_4_expr_filtering.RData")

# load("data/allbut18corresp.RData")
# load("data/allbut18corresp_d2t_peakexpr_4.RData")
# load("ages.RData")

qd = do.call(bind_rows, lapply(results, "[[", "quant.data")) 
sd = do.call(bind_rows, lapply(results, "[[", "segm.data"))
dd = do.call(bind_rows, lapply(results, "[[", "division.data"))
dd[which(dd$plant == 18), "d2t"] = NA
sd[which(sd$plant == 18), "d2t"] = NA
qd[which(qd$plant == 18), "d2t"] = NA

# qd = qd %>% select(-n.x, -n.y, -n.z)
# sd = sd %>% select(-n.x, -n.y, -n.z)

# data.tojose = bind_cols(sd, age=as.data.frame(ages)) %>% ungroup() %>% select(t, plant, id, layer, m.vol, n.vol, expr, m.x,m.y,m.z,n.x,n.y,n.z,d2t,circ,ages) %>% filter(!is.na(layer)) 
# colnames(data.tojose)[16] = "age"
# 
# write.table(data.tojose, row.names = FALSE, col.names = TRUE, file = "all_data_with_age_tojose.dat", quote = FALSE)
# for (ii in plants)
#   write.table(data.tojose %>% filter(plant == ii), row.names = FALSE, col.names = TRUE, file = paste0("all_data_with_age_tojose_plant", ii, ".dat"), quote = FALSE)

######################
cld = lapply(results, "[[", "lineage.sublines.data")
cld[[2]] = cld[[2]][sapply(cld[[2]], function(lineage) all(c(0, 76) %in% do.call(rbind, lineage)$t) & all(do.call(rbind, lineage)$layer == "L1"))]
cld[[3]] = cld[[3]][sapply(cld[[3]], function(lineage) all(c(0, 76) %in% do.call(rbind, lineage)$t) & all(do.call(rbind, lineage)$layer == "L1"))]
cld[[4]] = cld[[4]][sapply(cld[[4]], function(lineage) all(c(0, 84) %in% do.call(rbind, lineage)$t) & all(do.call(rbind, lineage)$layer == "L1"))]
cld[[5]] = cld[[5]][sapply(cld[[5]], function(lineage) all(c(0, 84) %in% do.call(rbind, lineage)$t) & all(do.call(rbind, lineage)$layer == "L1"))]
cld = cld[-c(1,6)]

var2 = "circ"
par(mfcol = c(6, 6), mar = c(1, 1, 1, 1)); for(ii in 1:18) {cat(ii, "\t"); plot.lineage(cld[[1]][[ii]], var2=var2, both = TRUE) }
par(mfcol = c(6, 6), mar = c(1, 1, 1, 1)); for(ii in 1:18) {cat(ii, "\t"); plot.lineage(cld[[2]][[ii]], var2=var2, both = TRUE) }
par(mfcol = c(6, 6), mar = c(1, 1, 1, 1)); for(ii in 1:18) {cat(ii, "\t"); plot.lineage(cld[[3]][[ii]], var2=var2, both = TRUE) }
par(mfcol = c(6, 6), mar = c(1, 1, 1, 1)); for(ii in 1:18) {cat(ii, "\t"); plot.lineage(cld[[4]][[ii]], var2=var2, both = TRUE) }

###############################################################################
### Plot data
###############################################################################
# Number of nuclei
# qd %>% 
  # group_by_(.dots = c("t", "plant")) %>% 
  # summarize(count = n(), mean.expr = mean(expr), median.expr = median(expr)) %>% gather(key, value) %>%
  # ggplot(aes(x = t, y = count, colour=plant)) + geom_line()
nNucl = qd %>%
group_by_(.dots = c("t", "plant")) %>% filter(d2t < 7.5) %>%
summarize(count = n()) %>%
  spread("t", "count") %>% apply(., 1, function(x)
  c(x[1], unlist(x[-1]))) %>% t %>% as.tibble
mExpr = qd %>%
  group_by_(.dots = c("t", "plant")) %>% filter(d2t < 7.5) %>%
  summarize(count = mean(expr, na.rm=T)) %>%
  spread("t", "count") %>% apply(., 1, function(x)
    c(x[1], unlist(x[-1]))) %>% t %>% as.tibble

dfm  = reshape2::melt(as.data.frame(nNucl), id.vars = c("plant"))
dfm2 = reshape2::melt(as.data.frame(mExpr), id.vars = c("plant"))
dfm %>% filter(plant !=18) %>% ggplot(aes(x = as.integer(variable) * 4, y = value, colour = "Number of nuclei")) +
  geom_line(data = dfm2 %>% filter(plant !=18), aes(x = as.integer(variable) * 4, y = value, colour = "Mean expression"), 
            linetype="dashed", size = 1) + 
  geom_line(size = 1) + 
  theme_bw() +
  facet_wrap( ~ factor(plant), ncol = 2, scales = "free") + 
  # scale_colour_brewer(type = "seq", palette = 1, direction = 1) +
  scale_colour_manual(name="",
                      values=c("Mean expression" = palette()[1], "Number of nuclei" = palette()[2])) +
  xlab('Time [h]') + 
  scale_y_continuous("Active CLV3 nuclei", 
                     sec.axis = sec_axis( ~ . * max(dfm2$value, na.rm = TRUE) /
                                                max(dfm$value,  na.rm = TRUE),
                                          name = "Mean CLV3 intensity")) +
  theme(
    legend.position = c(.68, .38),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "cm"),
    legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent")
  )

# e.dfm = reshape2::melt(as.data.frame(mExpr), id.vars = c("plant"))
# p2 = e.dfm %>% ggplot(aes(x=as.integer(variable), y = value)) + geom_line() + facet_wrap(~factor(plant), ncol = 2)
# grid.arrange(p1, p2)


# Distributions
plot.dist = function(qd, sd) {
  par(mfrow = c(3, 2))
  hist(qd$t,     breaks = 100, main = "t")
  hist(qd$expr,  breaks = 100, main = "expr")
  hist(qd$d2t,   breaks = 100, main = "d2t")
  hist(qd$n.vol, breaks = 100, main = "n.vol")
  hist(sd$m.vol, breaks = 5000, xlim = c(0,1000), main = "m.vol")
  hist(sd$neighs %>% sapply(length), breaks = 300, xlim = c(0,30), main = "nNeighs")
}
plot.dist(qd, sd)

# Cross-comparisons
plot.cc = function(qd) {
  par(mfrow = c(3, 2))
  boxplot(formula = qd$n.vol~qd$t)
  boxplot(formula = qd$expr~qd$t)
  boxplot(formula = qd$d2t~qd$t)
  plot(qd$n.vol, qd$expr)
  plot(qd$d2t,   qd$n.vol)
  plot(qd$d2t,   qd$expr)
}
plot.cc(qd)

# Division data
plot.dd = function(dd) {
  par(mfrow = c(3,2))
  hist(dd$t,     breaks = 30)
  hist(dd$age,   breaks = 100)
  hist(dd$n.vol, breaks = 50)
  hist(dd$m.vol, breaks = 50)
  hist(dd$expr,  breaks = 30)
  hist(dd$d2t,   breaks = 30)
}
plot.dd(dd)

###############################
# L1 analysis
###############################
qd.l1 = qd[qd$layer == "L1",  ]
qd.l2 = qd[qd$layer == "L2",  ]
qd.l3 = qd[qd$layer == "L3+", ]
sd.l1 = sd[sd$layer == "L1",  ]
sd.l2 = sd[sd$layer == "L2",  ]
sd.l3 = sd[sd$layer == "L3+", ]
dd.l1 = dd[dd$layer == "L1",  ]
dd.l2 = dd[dd$layer == "L2",  ]
dd.l3 = dd[dd$layer == "L3+", ]

plot.dist(qd.l1, sd.l1)
plot.dist(qd.l2, sd.l2)
plot.dist(qd.l3, sd.l3)

plot.cc(qd)
plot.cc(qd.l1)
plot.cc(qd.l2)
plot.cc(qd.l3)

plot.dd(dd.l1)
plot.dd(dd.l2)
plot.dd(dd.l3)

###############################
# CZ analysis
###############################
dist = 2
plot.dist(qd %>% filter(d2t < dist), sd)
plot.dist(qd.l1  %>% filter(d2t < dist), sd.l1)
plot.dist(qd.l2  %>% filter(d2t < dist), sd.l2)
plot.dist(qd.l3  %>% filter(d2t < dist), sd.l3)

plot.cc(qd)
plot.cc(qd    %>% filter(d2t < dist))
plot.cc(qd.l1 %>% filter(d2t < dist))
# plot.cc(qd.l2)
plot.cc(qd.l2 %>% filter(d2t < dist))
plot.cc(qd.l3 %>% filter(d2t < dist))

plot.dd(dd    %>% filter(d2t < dist))
plot.dd(dd.l1 %>% filter(d2t < dist))
plot.dd(dd.l2 %>% filter(d2t < dist))
plot.dd(dd.l3 %>% filter(d2t < dist))

####
# qd %>% group_by(t, plant) %>% summarise(n=n())
# 
p1 = ggplot(qd %>% filter(layer=="L1") %>% group_by(plant, t) %>% summarise(n = n()),
       aes(x = t, y = n)) +
  geom_line() + 
  geom_point() +
  facet_wrap( ~ plant) + 
  labs(title = "Nuclei per time")

p2 = ggplot(sd %>% filter(layer=="L1") %>%  filter(plant != 1) %>% group_by(t, plant) %>% summarize(n = n()),
            aes(x = t, y = n)) +
  geom_line() + geom_point() +
  facet_wrap(~ plant) +
  labs(title = "Cells per time")

p3 = ggplot(dd %>%  filter(layer=="L1") %>% filter(plant != 1) %>% group_by(t, plant) %>% summarize(n = n()),
       aes(x = t, y = n)) +
  geom_line() + geom_point() +
  facet_wrap(~ plant) +
  labs(title = "Divisions per time")

grid.arrange(p1, p2, p3, ncol = 1)

nCirc = sd %>% filter(plant != 1, plant != 18, layer == "L1", circ < 7.5) %>% 
  group_by(plant, circ) %>% 
  summarize(count = n()) %>% 
  spread("circ", "count") %>% as.tibble

nDivs  = dd %>% merge(sd, by = c("t", "plant", "id", "layer", "circ")) %>%
  select(t, plant, id, circ, layer) %>%
  filter(plant != 1, plant != 18, layer == "L1", circ < 7.5) %>%
  group_by(plant, circ) %>%
  summarize(count = n()) %>%
  spread("circ", "count") %>% as.tibble

mExpr  = dd %>% merge(sd, by = c("t", "plant", "id", "layer", "expr", "circ")) %>%
  select(t, plant, id, circ, layer, expr) %>%
  filter(plant != 1, plant != 18, layer == "L1", circ < 7.5) %>%
  group_by(plant, circ) %>%
  summarize(count = mean(expr, na.rm = TRUE)) %>% 
  spread("circ", "count") %>% as.tibble

mExpr = mExpr[,-1]
# nCirc = nCirc[, -((length(nCirc) - 1):length(nCirc))] # remove last 2
nDivs = nDivs / nCirc
nDivs = nDivs[,-1]

par(mfrow = c(1, 1))
mExpr.stat = mExpr %>% data.frame() %>% melt() %>% group_by(variable) %>% summarise(
  mean = mean(value, na.rm = T),
  sd = sd(value, na.rm = T),
  len = sum(!is.na(value))
) %>% group_by(variable) %>% mutate(sdm = sd / sqrt(len))


nDivs.stat = nDivs %>% data.frame() %>% melt() %>% group_by(variable) %>% summarise(
  mean = mean(value, na.rm = T),
  sd = sd(value, na.rm = T),
  len = sum(!is.na(value))
) %>% group_by(variable) %>% mutate(sdm = sd / sqrt(len))

p1 = ggplot(nDivs.stat,
       aes(
         x = as.integer(variable) - 1L,
         y = mean,
         ymin = mean - sdm,
         ymax = mean + sdm
       )) + 
  geom_errorbar(width = .2, size = 1.2) + geom_line(size = 1.2) + geom_point(size = 1.2) +
  labs(x = "", y = "Number of divisions (normalised)") + 
  theme_bw() + xlim(-0, 6) + theme(axis.title.x=element_blank(),
                                       axis.text.x=element_blank())
p2 = qd %>% 
  filter(layer=="L1") %>%
  ggplot(aes(x = as.factor(circ), y = expr)) + 
  geom_jitter(colour = "red", alpha = 0.1) +
  geom_violin() + 
  labs(x = "Distance to apex [#neighbours]", y = "Expression") + theme_bw()


p3 = sd %>% left_join(sd, by = c("t", "plant", "id", "layer", "m.vol", "circ")) %>%
  select(t, plant, id, circ, layer, m.vol, circ)  %>% filter(layer == "L1", circ < 7, plant != 1, plant != 18, !is.na(circ), !is.na(m.vol)) %>%
  select(t, plant, id, m.vol, circ) %>%
  group_by(circ) %>%
  summarize(mean = mean(m.vol), sd = sd(m.vol) / sqrt(n())) %>%
  ggplot(aes(x=circ, y = mean, ymin = mean - sd, ymax = mean + sd)) +
  geom_errorbar(width = .2, size = 1.2) +
  geom_line(size = 1.2) +
  theme_bw() +
  labs(x = "",
       y = "Mean membrane volume") + xlim(-0,6) + theme(axis.title.x = element_blank(),
                                                        axis.text.x  = element_blank())

# grid.arrange(p1, p3, p2)
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p3), ggplotGrob(p2), size = "last"))


nCells = sd %>% filter(plant != 1, layer=="L1") %>% group_by_(.dots = c("t", "plant")) %>% summarize(count = n()) %>% spread("t", "count") %>% apply(., 1, function(x) c(x[1], unlist(x[-1]))) %>% t %>% as.tibble
nNucl  = qd %>% filter(plant != 1, layer=="L1") %>% group_by_(.dots = c("t", "plant")) %>% summarize(count = n()) %>% spread("t", "count") %>% apply(., 1, function(x) c(x[1], unlist(x[-1]))) %>% t %>% as.tibble
nDivs  = dd %>% filter(plant != 1, layer=="L1") %>% group_by_(.dots = c("t", "plant")) %>% summarize(count = n()) %>% spread("t", "count") %>% apply(., 1, function(x) c(x[1], unlist(x[-1]))) %>% t %>% as.tibble


left_join(qd %>% filter(plant != 1, plant != 18) %>% group_by(t, plant, layer) %>% summarize(nNucl=n()),
dd %>% filter(plant != 1, plant != 18) %>% group_by(t, plant, layer) %>% summarize(nDivs=n()), by = c("t","plant", "layer")) %>% filter(layer == "L1") %>% 
  ggplot(aes(nNucl, nDivs, colour=factor(plant))) + geom_point(size = 2) + theme_bw() + labs(colour="Plant", x = "Number of divisions", y = "Number of CLV3 nuclei")


divs.vs.nucl = nDivs[-5,-1] %>% unlist %>% as.data.frame() %>% bind_cols(nNucl[,-c(1,ncol(nCells))]  %>% unlist %>% as.data.frame()) 
colnames(divs.vs.nucl) = c("nDivs", "nNucl")
divs.vs.nucl %>% ggplot(aes(x = nNucl, y = nDivs)) + geom_point(size = 2) + theme_bw()
cor.test(divs.vs.nucl[, 1], divs.vs.nucl[, 2])

# Correlations
cat("Indiv:\n")
for(ii in 1:nrow(nNucl)) print(cor.test(nCells[ii, -1] %>% unlist, 
                                         nNucl[ii, -1]  %>% unlist, na.action = na.omit())$p.value)
for(ii in 1:nrow(nNucl)) print(cor.test(nDivs[ii, -1] %>% unlist, 
                                        nNucl[ii, -c(1, ncol(nCells))]  %>% unlist, na.action = na.omit())$p.value)

print(cor.test(nCells[-5, -1] %>% unlist, nNucl[, -1] %>% unlist, na.action = na.omit())$p.value)
print(cor.test(nDivs[-5, -1] %>% unlist, nNucl[, -c(1, ncol(nCells))] %>% unlist, na.action = na.omit())$p.value)
# cor.test(nDivs[-5, -1] %>% unlist, nNucl[, -c(1, ncol(nCells))] %>% unlist, na.action = na.omit())

plot(nDivs[, -1] %>% unlist, nCells[, -c(1, ncol(nCells))] %>% unlist) # Not correlated!
plot(nDivs[-5, -1] %>% unlist, nNucl[, -c(1, ncol(nCells))] %>% unlist)  # Correlated!
plot(nNucl[, -1] %>% unlist, nCells[-5, -c(1)] %>% unlist)  # Correlated!



#######################################
# Nuclear vs. membrane volume
######################################3
p1 = sd %>% arrange(t, plant, layer) %>%  filter(layer == "L1") %>%
  ggplot(aes(x = m.vol, y = n.vol, colour = d2t)) + #xlim(0,350) +
  geom_point(na.rm = T)

# Per layer
p2 = sd %>% arrange(plant, t, layer) %>%
  ggplot(aes(x = m.vol, y = n.vol, colour = d2t)) + #xlim(0,350) +
  geom_point(na.rm = T) + facet_wrap( ~ layer) 
grid.arrange(p2, p1)

#######################################
# Violin plots
######################################3
dd %>% ggplot(aes(x=n.vol, y=expr, colour = d2t)) + geom_point(na.rm=T) + facet_wrap(~layer)
dd %>% filter(plant != 1) %>% ggplot(aes(x=factor(t), y=expr)) + geom_violin(na.rm=T) + facet_wrap(~plant) + labs(title="DIV:EXPR")
dd %>% filter(plant != 1) %>% ggplot(aes(x=factor(t), y=n.vol)) + geom_violin(na.rm=T) + facet_wrap(~plant) + labs(title="DIV:NVOL")
qd %>% filter(plant != 1) %>% ggplot(aes(x=factor(t), y=expr)) + geom_violin(na.rm=T) + facet_wrap(~plant) + labs(title="QD:EXPR")
qd %>% filter(plant != 1, layer=="L1") %>% ggplot(aes(x=factor(t), y=n.vol)) + geom_violin(na.rm=T) + facet_wrap(~plant) + labs(title="QD:NVOL")
qd %>% filter(layer =="L1") %>% ggplot(aes(x=factor(t), y=n.vol)) + geom_boxplot(na.rm=T)  + labs(title="QD:NVOL")
sd %>% filter(layer == "L1", m.vol < 300) %>% ggplot(aes(x=factor(t), y=m.vol)) + geom_boxplot(na.rm=T) + facet_wrap(~plant) + labs(title="SD:MVOL")

#######################################
# BIN EXPR-D2T IN L1
#######################################
qd.l1 = qd %>% filter(layer == "L1", plant != 1, plant != 18) %>% arrange(d2t)
qd.l1 = qd.l1 %>% filter(plant != 1, plant != 18, d2t < 8, !(d2t < 2 & expr < 50), plant ==2)
qd.l1 %>% ggplot(aes(x = d2t, y = expr, colour = d2t)) +
  geom_point(size = 2) + 
  theme_bw() + 
  labs(x = "Distance to center", y = "CLV3 expression") + facet_wrap(~plant)
nBins = 14
bin.l1 = add_column(as.data.frame(qd.l1), bin = 
                      cut(as.data.frame(qd.l1)$d2t, 
                          breaks = seq(min(as.data.frame(qd.l1)$d2t, na.rm = TRUE),
                                       max(as.data.frame(qd.l1)$d2t, na.rm = TRUE),
                                       length.out = nBins), 
                          include.lowest = TRUE, 
                          labels = 1:(nBins - 1))) %>% as.tibble
binned.data = bin.l1 %>% group_by(bin) %>% summarize(mean=mean(expr, na.rm=T), sd = sd(expr, na.rm=T) / (mean(expr, na.rm=T)), pos = mean(d2t, na.rm=T)) 
binned.data$bin = binned.data$bin %>% as.integer
ggplot(binned.data, aes(x=bin, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sd*mean, ymax=mean+sd*mean), colour="black", width=.1) +
  geom_line() +
  geom_point() + labs(x="Bin", y = "Mean expression")
# save(file="binned_d2t_12.RData", binned.data)
par(mar = c(4, 4, 4, 4))
plot(binned.data$sd)

# Violin plot of thingy
p1 = qd %>% 
  filter(layer=="L1") %>% 
  mutate(bin = cut(
  d2t,
  seq(0, max(d2t, na.rm=T), length.out = nBins),
  labels = 1:(nBins - 1),
  include.lowest = TRUE
)) %>% 
  group_by(bin) %>% 
    ggplot(aes(x = bin, y = expr)) +
    geom_jitter(colour="red", alpha=.1) +
    geom_violin() + 
    facet_wrap( ~ plant) + 
    labs(x = "Distance to apex (binned)", y = "Expression", title = "Indiv plants") 

p2 = qd %>% 
  filter(layer=="L1") %>%
  mutate(bin = cut(
  d2t,
  seq(0, max(d2t, na.rm=T), length.out = nBins),
  labels = 1:(nBins - 1),
  include.lowest = TRUE
)) %>%
  group_by(bin) %>% 
  ggplot(aes(x = bin, y = expr)) + 
  geom_jitter(colour = "red", alpha = 0.1) +
  geom_violin() + 
  labs(x = "Distance to apex (binned)", y = "Expression", title = "All")
grid.arrange(p1, p2)


# qd %>% filter(plant != 18) %>% ggplot(aes(x=d2t, y=expr)) +  stat_binhex(bins = 40)
# qd %>% ggplot(aes(x=d2t, y=expr)) +  stat_binhex() + facet_wrap(~layer)

########################################
########################################
########################################
########################################
par(mfcol = c(2, 2))
for(ii in 1:4) {
  plant.2 = qd %>% filter(layer == "L1") %>% group_by(plant, t) %>%
    summarise(nCells = n()) %>% group_by(plant, t) %>%
    spread(t, nCells) %>% .[,-1] %>% .[ii,] %>% unlist() 
  # plant.2 = plant.2[!is.na(plant.2)]
  times = as.integer(names(plant.2))
  
  # trend = glm(log(plant.2) ~ times)
  # trend = loess(plant.2 ~ times, degree = 1)
  trend = loess(plant.2 ~ times, degree = 2)
  # detrended.trajectory = (trend$residuals / lag(trend$residuals, 1) - 1)[-1]
  # detrended.trajectory = (trend$residuals) %>% sign()
  detrended.trajectory = (trend$residuals)
  # detrended.trajectory = plant.2
  
  # plot(plant.2, type = "b")
  # lines(predict(trend), col = "red")
  plot(detrended.trajectory, type = "b", xaxt="n", yaxt="n")
  # lines(mean(plant.2) + 10*sin(2*pi/3*1:22), col = "red")
  f.data <- GeneCycle::periodogram(detrended.trajectory)
  f.data = data.frame(do.call(cbind, f.data)); colnames(f.data) = c("spec", "freq")
  
  # harmonics <- 1:8
  p = ggplot(f.data, aes(x = freq*length(plant.2) * 4,
                         y = spec / sum(spec), fill = factor(freq * length(plant.2) * 4))) + 
    geom_bar(stat = "identity") +
    # geom_point() +
    labs(x = "Harmonic [h]", y = "Amplitude density")  +
    scale_fill_manual(breaks = c(".2", ".5"), values=c(rep(2, 3), 1, rep(2, length(plant.2)-4))) + theme_bw()
  assign(paste0("p", ii), p, inherits = FALSE)
  print(randtests::bartels.rank.test(plant.2, alternative = "r")$p.value)
}
grid.draw(cbind(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"), rbind(ggplotGrob(p3), ggplotGrob(p4), size = "last"),
                size="last")  )
# grid.arrange(p1, p3, p2, p4, nrow = 2)

#################
### FINAL VIOLIN
#################
sd %>% filter(!is.na(expr), layer == "L1") %>% ggplot(aes(x=factor(circ), y = expr)) + geom_jitter(alpha=.1, colour="red") + geom_violin() + facet_wrap(~plant) + labs(title = "L1", x = "Distance in cells from top-expressing cell", y = "Expression")

p1 = sd %>% filter(plant != 1, plant != 18, !is.na(circ), !is.na(m.vol), layer == "L1") %>% 
  ggplot(aes(x = factor(circ), y = m.vol)) + 
  geom_jitter(alpha = .1, colour = "red") + 
  geom_violin() + 
  facet_wrap( ~ plant)
p2 = sd %>% filter(plant != 1, plant != 18, !is.na(circ), !is.na(m.vol), layer == "L1") %>% 
  ggplot(aes(x = factor(circ), y = m.vol)) + 
  geom_jitter(alpha = .1,  colour = "red") + 
  geom_violin() 
grid.arrange(p1, p2)

# Membrane volume violin
# sd %>% filter(layer == "L1", circ < 6, m.vol < 300, plant != 1, plant != 18, !is.na(circ), !is.na(m.vol)) %>% ggplot(aes(x=factor(circ), y = m.vol)) + geom_jitter(alpha = .1, colour="red") + geom_violin() + facet_wrap(~plant) 

p1 = sd %>% filter(layer == "L1", circ < 7, plant != 1, plant != 18, !is.na(circ), !is.na(m.vol)) %>% 
  select(t, plant, id, m.vol, circ) %>%
  group_by(circ) %>% 
  summarize(mean = mean(m.vol), sd = sd(m.vol) / sqrt(n())) %>% 
  ggplot(aes(x=circ, y = mean, ymin = mean - sd, ymax = mean + sd)) + 
  geom_errorbar(width = .2, size = 1.5) +
  geom_line(size = 1.5) +
  theme_bw() +
  labs(x = "Distance to apex [#neighbours]", 
       y = "Mean membrane volume")

# p2 =
# grid.arrange(p1, p2)


######################
# Sigmoid CLV3?
######################
qd %>% filter(layer == "L1") %>% 
  group_by(circ) %>% 
  summarize(m = mean(expr, na.rm = T), sd = sd(expr, na.rm = T)) %>%
  ungroup() %>% 
  ggplot(aes(x = circ, y = m, ymin = m - sd, ymax = m + sd)) + 
    geom_errorbar(width = .1, na.rm = T) +
    geom_point() + 
    geom_line()

qd %>%
  filter(plant == 2, layer == "L1") %>%
  group_by(circ) %>% 
  summarize(m = mean(expr, na.rm = T), sd = sd(expr, na.rm = T)) %>%
  ungroup() %>% 
  ggplot(aes(x = circ, y = m, ymin = m - sd, ymax = m + sd)) + 
    geom_errorbar(width = .1, na.rm = T) + 
    geom_point() + 
    geom_line()


##################
# Add age to sd
##################
# subs = 
# cl =  makeCluster(3, outfile = "")
# # registerDoParallel(cl)
# clusterEvalQ(cl, library(tidyverse))
# clusterExport(cl, list = ls(), envir = environment())
# ages = parSapply(cl, 1:nrow(sd), function(x) {get.age(sd$t[x], sd$id[x], results[[which(plants == sd$plant[x])]]$sublines)})

