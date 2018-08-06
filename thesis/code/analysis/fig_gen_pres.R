setwd("/home/henrik/compbio/thesis/code/newTreat/")
library(plyr)
library(tidyverse)
source("functions.R")
source("read_data.R")
source("aux.R")
source("../plot.R")
plants = c(1, 2, 4, 13, 15, 18)
p = c(2, 4, 13, 15)
theme_set(theme_bw())

############################################################
### NUMBER OF NUCLEI
############################################################
# load("data/allbut18corresp_d2t_peakexpr_4.RData")
# qd = do.call(bind_rows, lapply(results, "[[", "quant.data"))
load("data/filtered_topexpr4.RData")

# manually corrected tp: 2: 44, 52; 4: 24

nNucl = qd %>% 
  # filter(!(plant == 2 & (t == 44 | t == 52)), !(plant == 4 & t == 24)) %>%
  group_by(t, plant) %>%
  summarize(count = n(), 
            expr = mean(expr, na.rm = TRUE)) %>%
  filter(plant %in% p) %>% 
  ungroup()

ggplot(nNucl, group = plant) + 
  # geom_point(aes(x = t, y = count, colour = factor(1, levels = c(1,2))), size = 2.5, alpha = .8) + 
  # geom_line(aes(x  = t, y = count, colour = factor(1, levels = c(1,2))), size = 1,   alpha = .8) + 
  # geom_point(aes(x = t, y = expr,  colour = factor(alpha(2,1), levels = c(1,2))), size = 2.5, alpha = .8) + 
  # geom_line(aes(x  = t, y = expr,  colour = factor(alpha(2,1), levels = c(1,2))), size = 1,   alpha = .8) + 
  # geom_point(aes(x = t, y = expr,  colour = "red"), size = 2.5, alpha = .5) + 
  # geom_line(aes(x  = t, y = expr,  colour = "red"), size = 1,   alpha = .5) + 
  geom_point(aes(x = t, y = count, colour = factor(plant)), size = 2.5, alpha = .8) +
  geom_line(aes(x  = t, y = count, colour = factor(plant)), size = 1,   alpha = .8) +
  facet_wrap( ~ plant, scales = "free") +
  theme(text = element_text(size = 22),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.key.size = unit(.5, "cm")
  ) + 
  scale_color_discrete(guide = FALSE) +
  # scale_color_discrete(name="", labels=c("Number CLV3 nuclei", "Mean CLV3 expression")) +
  # scale_y_continuous(
  # sec.axis = sec_axis( ~ . * max(nNucl$expr, na.rm = TRUE) / max(nNucl$count, na.rm = TRUE), 
  # name = "Mean CLV3 expression")) +
  labs(x = "Time [hours]", y = "Number of CLV3 nuclei")

for(ii in p){
  dat= nNucl %>% filter(plant == ii) 
  print(cor.test(dat$count, dat$expr)$p.value)
}

#############
# nDivs with nNucl (layer 1)
data = 
  left_join(
    qd %>% filter(plant != 1, plant != 18) %>% group_by(t, plant) %>% 
      summarize(nNucl = n()),
    dd %>% filter(plant != 1, plant != 18) %>% group_by(t, plant) %>% 
      summarize(nDivs = n()), by = c("t", "plant")) %>% ungroup()
data = 
  left_join(data,
            sd %>% filter(plant != 1, plant != 18) %>% group_by(t, plant) %>%
              summarize(nCells = n()), by = c("t", "plant")) %>% 
  mutate(norm = nDivs / nCells) %>% mutate(normNucl = nNucl / nCells)

cor.test(data$normNucl, data$norm)

for (ii in p) {
  dat = data %>% filter(plant == ii)
  print(cor.test(dat$normNucl, dat$norm))
}

library(MASS)
data %>% ungroup() %>%
  ggplot(data = ., aes(norm, normNucl)) +
  stat_smooth(method = "rlm", se = FALSE, size = 3, colour = "darkgray") +
  geom_point(
    aes(colour = factor(plant)),
    size = 3,
    alpha = .8,
    na.rm = TRUE
  ) + theme_bw()  +
  labs(colour = "Plant", x = "Number of divisions / total cells", y = "Number of CLV3 nuclei / total cells") +
  theme(text = element_text(size = 22), legend.position = "top")


############################################################
### Longevity, age
############################################################
# load("data/d2t_comparisons.RData")
# dd = do.call(bind_rows, lapply(results[[3]], "[[", "division.data"))
load("data/filtered_topexpr4.RData")
dd %>% filter(!(plant == 18 & t == 40)) %>% 
  ggplot(aes(x = age)) +
  geom_histogram(aes(y = ..density.., colour = 5), bins = 21, na.rm = T) + 
  geom_histogram(aes(y = ..density.., fill = 1, colour = 3), bins = 21, na.rm = T) + 
  facet_grid( ~ layer) + 
  theme_bw() + 
  labs(x="Division age [hours]", y = "Density") +
  theme(legend.position="none") 

dd %>% filter(!(plant == 18 & t == 40)) %>% 
  ggplot(aes(x = age)) +
  geom_histogram(aes(y = ..density.., fill = 1), bins = 21, na.rm = T) + 
  facet_grid( plant~ layer) + 
  theme_bw() + 
  labs(x="Division age [hours]", y = "Density") +
  theme(legend.position="none") 

############################################################
### Nuclear volume. Epidermal regulation?
############################################################
# load(file = "data/allbut18corresp.RData_d2t_peakexpr_4_expr_filtering.RData")
load("data/filtered_topexpr4.RData")
p1 = qd %>% filter(plant != 18, plant == 13) %>%
  ggplot(aes(x = n.vol, y = expr, colour = factor(plant), alpha = d2t)) + 
  geom_point(size = 2) +
  facet_grid(plant ~ layer) + 
  xlim(0,1) +
  labs(x = "\nNuclear volume", y = "CLV3 expression") +
  theme(legend.position = "bottom", text = element_text(size=22)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   margin = margin(t = 0, r = 0, b = 0, l = 0))) +
  scale_color_discrete(name = "", guide=FALSE) +
  scale_alpha_continuous(name = d2ap, range=c(0.1,.5), guide = FALSE)

p2 = qd %>% filter(plant %in% p, layer == "L1") %>%
  ggplot(aes(x = d2t, y = expr, colour = factor(plant), alpha = t)) +
  geom_point( na.rm = TRUE) +
  facet_wrap(~ plant) +
  labs(x = d2ap, y = "CLV3 expression") + 
  scale_y_continuous(position = "left") +
  scale_alpha_continuous(name = "", guide = FALSE, range=c(0.1,0.5)) +
  scale_color_discrete(name = "Plant", guide = FALSE) + 
  theme(legend.position = "bottom") 

p3 = qd %>% 
  filter(plant == 13, !is.na(layer), plant %in% p, layer %in% c("L1", "L2")) %>%
  ggplot(aes(x = d2t, y = expr, colour = factor(plant), alpha = t)) +
  geom_point(na.rm = TRUE) +
  facet_grid(plant~ layer) +
  labs(x = expression(paste("Distance to apex [", mu, "m]")), y = "CLV3 expression") + 
  scale_y_continuous(position = "left") +
  scale_alpha_continuous(name = "Time [hours]",  range=c(0.1,0.5), guide = FALSE) +
  scale_color_discrete(name = "Plant", guide = FALSE) + 
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.ticks.y =element_blank(), 
        axis.text.y = element_blank(), text = element_text(size=22))
# p2
# grid.arrange(p1, p3, ncol = 2)
grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p3), size = "last"))



p1############################################################
### Quality score
############################################################
quality = lapply(plants, get.quality.score, cutoff = 0)
# load("")
quality = quality %>% 
  do.call(bind_rows, .) %>%
  left_join(sd, by = c("plant", "t", "id")) %>% 
  select(plant, t, id, score, layer) %>% 
  filter(!is.na(layer)) 
quality$score[which(quality$plant == 1)] = quality$score[which(quality$plant == 1)] / 2
quality$score[which(quality$plant == 1 & quality$t < 16)] = quality$score[which(quality$plant == 1 & quality$t < 16)] * 2

quality %>%  
  ggplot(aes(x = layer, y = score)) + 
  geom_jitter(col = 1, alpha = .1) + 
  geom_violin() + 
  facet_wrap(~plant) + 
  labs(x = "Layer", y = "Quality score") + 
  theme_bw()

quality %>%
  ggplot(aes(x = factor(t), y = score)) + 
  geom_jitter(col = 1, alpha = .1) + 
  geom_violin() + 
  facet_wrap(~plant) + 
  labs(x = "Time", y = "Quality score") + 
  theme_bw() + theme(axis.text.x=element_text(angle = 90, hjust = 1))

############################################################
### CLV3 Layer sep
############################################################
# load(file="data/allbut18corresp.RData_d2t_peakexpr_4_expr_filtering.RData")
load("data/filtered_weightexpr4.RData")

qd %>% 
  filter(!is.na(layer), layer != "L3+") %>% 
  ggplot(aes(x = d2t, y = expr, colour = factor(plant))) + 
  geom_point(alpha = .5) + 
  facet_wrap( ~ layer) + 
  scale_color_discrete(name = "Plant") + 
  theme_bw() + labs(x = "Apical distance", y = "CLV3 Expression") 


########################
# Periodograms
########################
# load(file = "data/allbut18corresp.RData_d2t_peakexpr_4_expr_filtering.RData")
load("data/filtered_topexpr4.RData")

data = qd %>% 
  filter(layer == "L1") %>% 
  group_by(plant, t) %>%
  summarise(nCells = n()) %>%
  ungroup()
trends = data %>% 
  group_by(plant) %>%
  do(residuals = loess(nCells ~ t, data = .)$residuals)

# all.trends = do.call(rbind.fill, do.call(as.data.frame, trends$residuals))
times = data$t %>% unique() %>% sort()
lens = max(sapply(trends$residuals, length))
res = trends$residuals
res = lapply(res, function(x) {
  length(x) = lens
  x
})
res = do.call(rbind, res)
colnames(res) = times
res = res %>% t() %>% cbind(times, .) %>% as.tibble()
colnames(res) = c("t", "2", "4", "13", "15")
melted = melt(res, id.vars = c("t"))
ggplot(melted, aes(x = as.integer(t), y = value, group = variable, colour = variable)) + 
  geom_line(size=1, alpha = .8) + 
  geom_point(size=3, alpha = .8) +
  facet_wrap(~variable) +
  scale_color_discrete(name="Plant", guide = FALSE) +
  labs(x = "Time [hours]", y = "Number of CLV3 nuclei\n (detrended)") +
  theme(text = element_text(size=22))


f.plants = c(2, 4, 13, 15)
lens = sapply(trends$residuals, length)
f.data = lapply(trends$residuals, function(x) data.frame(do.call(cbind,GeneCycle::periodogram(x))))
f.data = lapply(1:length(f.data), function(x)
  cbind(plant = rep(f.plants[x], nrow(f.data[[x]])), spec=f.data[[x]][, 1], freq=f.data[[x]][,2] * lens[x]*4)) %>% 
  do.call(rbind, .) %>% 
  as.tibble() %>% group_by(plant) %>%  mutate(spec = spec/sum(spec)) %>% ungroup() 
f.data = ddply(f.data, .(plant), mutate, id = seq_along(plant)) %>% as.tibble()


col1 = gg_color_hue(2)[1] #alpha("red", .7)
col2 = gg_color_hue(2)[2] #alpha("blue", 1)
ggplot(f.data, aes(
  x = f.data$freq,
  y = f.data$spec,
  fill = factor(id)
)) +
  geom_bar(stat = "identity", width = 3.5)  +
  theme_bw() +
  facet_wrap(~ plant) +
  scale_fill_manual(
    name = "",
    breaks = c(3.5, 4.5),
    values = c(col1, col1, col1, col1, col1, col1, col1, col1, col1, col1, col1)
  ) +
  labs(x = "Harmonic [hours]", y = "Amplitude density") + 
  theme(text = element_text(size=22))

########################
# CLV3 distribution
########################
# load(file = "data/allbut18corresp.RData_d2t_peakexpr_4_expr_filtering.RData")
load("data/filtered_topexpr4.RData")
data = qd %>% 
  filter(plant %in% p, layer == "L1") %>% 
  mutate(colour = factor(plant))

data %>% filter(plant == 2) %>% 
  ggplot(aes(x = d2t, y = expr, group = plant, colour = t)) +
  geom_point(na.rm = TRUE, size = 3, alpha = .7) +
  facet_wrap(~ plant) +
  theme_bw() +
  labs(x = expression(paste("Distance to apex [", mu, "m]")), y = "CLV3 expression") +
  theme(text=element_text(size=22))+
  scale_color_continuous(name = "Time")
  # scale_color_discrete(name= "Plant", guide = FALSE) +
  # scale_alpha_continuous(name = "Time [hours]")

# scale_alpha_discrete(range = c(0.1, .4))
d2ap = expression(paste("Distance to apex [", mu, "m]"))

# scale_fill_manual(labels = t, values = "red", "blue", "green", "brown")
# scale_alpha_manual(values = c(0.6, 1)) 

# scale_fill_manual(labels="plant", values=c("red","blue","green","brown")) 

# scale_color_continuous(name = "Time [hours]") + theme(legend.position = "bottom")

# scale_colour_identity() + 

########################
# Trajectories
########################
load(file = "data/filtered_topexpr4.RData")

# cld = lapply(results, "[[", "lineage.sublines.data")
cld[[2]] = cld[[2]][sapply(cld[[2]], function(lineage) 
  all(c(0, 76) %in% do.call(rbind, lineage)$t) & all(do.call(rbind, lineage)$layer == "L1"))]
cld = cld[[2]]
lineage = 1:length(cld)
cld = lapply(1:length(cld), function(x) 
  lapply(1:length(cld[[x]]), function(y) { 
    cbind(lineage = rep(lineage[x], nrow(cld[[x]][[y]])), 
          subline = rep(y, nrow(cld[[x]][[y]])), cld[[x]][[y]]) %>% as.tibble()
  }))
cld = do.call(rbind, unlist(cld, recursive = FALSE))
ggplot(cld, aes(x=t, y = expr, col = factor(subline), group = factor(subline))) + 
  geom_line(na.rm = T, alpha = 1) + 
  geom_point(na.rm = TRUE, size = 1) +
  # geom_line(aes(x=t,y=d2t, col = factor(subline)), na.rm = T, alpha = .5) + 
  # geom_point(aes(x=t,y=d2t, col = factor(subline)), na.rm = TRUE, size = .5) +
  facet_wrap(~lineage, ncol = 6) + 
  theme(legend.position="none") +
  labs(x = "Time [hours]", y = "CLV3 expression")

### behaviours
select.lines = c(5, 7)
p1 = ggplot(cld %>% filter(lineage %in% select.lines), 
            aes(x = t, y = expr, col = factor(subline), group = factor(subline))) + 
  geom_line(na.rm = T, alpha = 1, size = 1) + 
  geom_point(na.rm = TRUE, size = 3) +
  facet_wrap(~lineage) + 
  theme_bw() +
  theme(legend.position="none", text=element_text(size=22)) + 
  labs(x = "", y = "CLV3 expression") +   theme(axis.title.x=element_blank(),
                                                axis.text.x=element_blank(),
                                                axis.ticks.x=element_blank()) +
  scale_color_discrete(guide=FALSE)
p2 = ggplot(cld %>% filter(lineage %in% select.lines), 
            aes(x = t, y = d2t, col = factor(subline), group = factor(subline))) + 
  geom_line(na.rm = T, alpha = 1, size = 1) + 
  geom_point(na.rm = TRUE, size = 3, alpha = 1) +
  facet_wrap(~lineage) + 
  theme(legend.position="none", strip.text.x = element_blank()) +
  labs(x = "Time [hours]", y = "Dist. to apex") + theme_bw() +
  scale_color_discrete(guide=FALSE) + 
  theme(legend.position="none", text=element_text(size=22)) + 
  theme(strip.background = element_blank(), strip.text.x = element_blank())

grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last")) 
# grid.arrange(p1,p2, ncol = 1)

########################
# Decay rate
########################
# load(file = "data/allbut18corresp.RData_d2t_peakexpr_4_expr_filtering.RData")
# load("data/filtered_topexpr4.RData")

# cld = lapply(results, "[[", "lineage.sublines.data")
pl = 2
cld[[pl]] = cld[[pl]][sapply(cld[[pl]], function(lineage) all(do.call(rbind, lineage)$layer == "L1"))]
# cld[[2]] = cld[[2]][sapply(cld[[2]], function(lineage) all(c(0, 76) %in% do.call(rbind, lineage)$t) & all(do.call(rbind, lineage)$layer == "L1"))]
cld = cld[[pl]]
lineage = 1:length(cld)
cld = lapply(1:length(cld), function(x) 
  lapply(1:length(cld[[x]]), function(y) { 
    cbind(lineage = rep(lineage[x], nrow(cld[[x]][[y]])), 
          subline = rep(y, nrow(cld[[x]][[y]])), cld[[x]][[y]]) %>% as.tibble()
  }))
cld = do.call(rbind, unlist(cld, recursive = FALSE))

data = cld %>% 
  select(lineage, subline, plant, t, d2t, expr, circ) %>%
  filter(circ > 1) %>%
  group_by(lineage) %>%
  mutate(t = t - min(t, na.rm = TRUE), expr = expr / max(expr, na.rm = TRUE)) %>% 
  ungroup()

subs = table(interaction(data$lineage, data$subline))
longer.than = names(subs[subs > 12])
data = data %>% 
  filter(interaction(lineage, subline) %in% longer.than)  

fit = data %>%
  group_by(lineage, subline, plant) %>%  na.omit() %>%
  do(fit = lm(data=., log(expr) ~ t)) 

library(Rmisc)
est = do.call(rbind, sapply(fit$fit, "[", "coefficients"))[, 2] %>% na.omit()
est = CI(est)
fit = data %>%
  group_by(lineage, subline, plant) %>%  na.omit() %>%
  do(fit = lm(data=., expr ~ t)) 

lest = do.call(rbind, sapply(fit$fit, "[", "coefficients"))[, 2] %>% na.omit()
lest = CI(lest)

d   = data.frame(d.x = 0:76, d.y = exp(est[2] * 0:76), upper = exp(est[1] * 0:76), lower = exp(est[3] * 0:76)) %>% mutate(label = 1)
d2  = data.frame(d.x = 0:76, d.y = 1 + lest[2] * 0:76, upper = 1 + lest[1] * 0:76, lower = 1 + lest[3] * 0:76) %>% mutate(label = 2)

p1 = ggplot() +
  geom_point(data=data,aes(x = t, y = expr, colour = factor(lineage), 
                           group = interaction(lineage, subline)), na.rm = TRUE, alpha = 0.2) + 
  geom_line(data=data,aes(x = t, y = expr, colour = factor(lineage),
                          group = interaction(lineage, subline)), na.rm  = TRUE, alpha = 0.2) + 
  geom_line(data = d, aes(x = d.x, y = d.y), size = 1, colour= "blue") +
  geom_ribbon(data = d, aes(x = d.x, ymin = lower, ymax = upper), alpha = 0.3) + 
  geom_text(aes(x = 58, y = .95, label = paste("coeff =", round(est[2], 3), "\u00B1", round(est[2]-est[1], 3))), parse = FALSE, size = 4) +
  theme(legend.position = "none", strip.text.x = element_blank())  +
  labs(x = "", y = "Normalised CLV3 expression")

p2 = ggplot() +
  geom_point(data=data,aes(x = t, y = expr, colour = factor(lineage), 
                           group = interaction(lineage, subline)), na.rm = TRUE, alpha = 0.2) + 
  geom_line(data=data,aes(x = t, y = expr, colour = factor(lineage),
                          group = interaction(lineage, subline)), na.rm  = TRUE, alpha = 0.2) + 
  geom_line(data = d2, aes(x = d.x, y = d.y), size = 1, colour= "blue") +
  geom_ribbon(data = d2, aes(x = d.x, ymin = lower, ymax = upper), alpha = 0.3) + 
  geom_text(aes(x = 58, y = .95, label = paste("coeff =", round(lest[2], 3), "\u00B1", round(lest[2]-lest[1], 3))), parse = FALSE, size = 4) +
  theme(legend.position = "none", strip.text.x = element_blank(), strip.text.y = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())  +
  labs(x = "", y = "")
grid.arrange(p1, p2, ncol = 2,   bottom = textGrob("Time [hours]"))

####################################
### Topcells m.vol
####################################
load(file = "data/filtered_topexpr4.RData")
# load("data/filtered_weightexpr3.RData")
# qd = do.call(bind_rows, lapply(results, "[[", "quant.data")) 
# sd = do.call(bind_rows, lapply(results, "[[", "segm.data"))
# dd = do.call(bind_rows, lapply(results, "[[", "division.data"))
dd[which(dd$plant == 18), "d2t"] = NA
sd[which(sd$plant == 18), "d2t"] = NA
qd[which(qd$plant == 18), "d2t"] = NA

avg.pop = sd %>% 
  filter(plant != 1, plant != 18, !is.na(circ), !is.na(m.vol), layer == "L1", circ < 7) %>% 
  group_by(circ) %>% 
  summarize(mean = mean(m.vol, na.rm = TRUE), se = sd(m.vol, na.rm = TRUE)/sqrt(n()))  %>% mutate(plant="All")
avg.plant = sd %>% 
  filter(plant != 1, plant != 18, !is.na(circ), !is.na(m.vol), layer == "L1", circ < 7) %>% 
  group_by(plant, circ) %>% 
  summarize(mean = mean(m.vol, na.rm = TRUE), se = sd(m.vol, na.rm = TRUE)/sqrt(n()))   %>% arrange(plant)
violin.data = sd %>% 
  filter(plant != 1, plant != 18, !is.na(circ), !is.na(m.vol), layer == "L1", circ < 7)  %>% arrange(plant)

p1 = ggplot() + 
  geom_jitter(data = violin.data, aes(x = factor(circ), y = m.vol, colour = factor(plant)), alpha = .3) + 
  geom_violin(data = violin.data, aes(x = factor(circ), y = m.vol)) + 
  # ylim(0, 400) +
  geom_boxplot(data = violin.data, aes(x=factor(circ), y = m.vol), width = .5) + 
  facet_wrap( ~ plant) +
  ylim(0, 400) +
  labs(y = "Membrane volume", x = "") +
  scale_color_discrete(guide = FALSE)

p2 = ggplot() + 
  geom_jitter(data = transform(violin.data, plant = "All"), 
              aes(x = factor(circ), y = m.vol), alpha = .1, colour = 5) + 
  geom_violin(data = transform(violin.data, plant = "All"), 
              aes(x = factor(circ), y = m.vol)) + 
  geom_point(data = avg.pop, aes(x=factor(circ), y = mean), size = 1) +
  geom_errorbar(data = avg.pop, aes(x=factor(circ), ymin = mean-se, ymax = mean+se, width= .1)) + 
  geom_boxplot(data = transform(violin.data, plant = "All"), 
               aes(x=factor(circ), y = m.vol), width = .5) + 
  facet_wrap(~plant) +
  ylim(0, 400) +
  labs(y = "", x = "") 
grid.arrange(
  p1,
  p2,
  ncol = 2,
  bottom = textGrob("Distance to apex [#neighbours]")
)

####################################
### Topcells n.vol
####################################
load(file = "data/filtered_topexpr4.RData")
# qd = do.call(bind_rows, lapply(results, "[[", "quant.data")) 
# sd = do.call(bind_rows, lapply(results, "[[", "segm.data"))
# dd = do.call(bind_rows, lapply(results, "[[", "division.data"))
dd[which(dd$plant == 18), "d2t"] = NA
sd[which(sd$plant == 18), "d2t"] = NA
qd[which(qd$plant == 18), "d2t"] = NA

avg.pop = sd %>% 
  filter(plant != 1, plant != 18, !is.na(circ), !is.na(n.vol), layer == "L1", circ < 7) %>% 
  group_by(circ) %>% 
  summarize(mean = mean(n.vol, na.rm = TRUE), se = sd(n.vol, na.rm = TRUE)/sqrt(n()))  %>% mutate(plant="All")
avg.plant = sd %>% 
  filter(plant != 1, plant != 18, !is.na(circ), !is.na(n.vol), layer == "L1", circ < 7) %>% 
  group_by(plant, circ) %>% 
  summarize(mean = mean(n.vol, na.rm = TRUE), se = sd(n.vol, na.rm = TRUE)/sqrt(n()))   %>% arrange(plant)
violin.data = sd %>% 
  filter(plant != 1, plant != 18, !is.na(circ), !is.na(n.vol), layer == "L1", circ < 7)  %>% arrange(plant)

for (ii in p) {
  # print((sd %>% filter(plant == ii, layer == "L1", circ == 0) %>% .$m.vol %>% mean(na.rm=TRUE) -
  # sd %>% filter(plant == ii, layer == "L1", circ == 1) %>% .$m.vol %>% mean(na.rm=TRUE)))
  print(t.test(sd %>% filter(plant == ii, layer == "L1", circ == 0) %>% .$n.vol, 
               sd %>% filter(plant == ii, layer == "L1", circ == 1) %>% .$n.vol)$p.value)
}

p1 = ggplot() + 
  geom_jitter(data = violin.data, aes(x = factor(circ), y = n.vol, colour = factor(plant)), alpha = .3) + 
  geom_violin(data = violin.data, aes(x = factor(circ), y = n.vol)) + 
  # ylim(0, 400) +
  geom_boxplot(data = violin.data, aes(x=factor(circ), y = n.vol), width = .3) + 
  facet_wrap( ~ plant) +
  ylim(0, 1) +
  labs(y = "Nuclear volume", x = "")  +
  scale_color_discrete(guide = FALSE)
p2 = ggplot() + 
  geom_jitter(data = transform(violin.data, plant = "All"), 
              aes(x = factor(circ), y = n.vol), alpha = .1, colour = 5) + 
  geom_violin(data = transform(violin.data, plant = "All"), 
              aes(x = factor(circ), y = n.vol)) + 
  geom_point(data = avg.pop, aes(x=factor(circ), y = mean), size = 1) +
  geom_errorbar(data = avg.pop, aes(x=factor(circ), ymin = mean-se, ymax = mean+se, width= .1)) + 
  geom_boxplot(data = transform(violin.data, plant = "All"), 
               aes(x=factor(circ), y = n.vol), width = .3) + 
  facet_wrap(~plant) +
  ylim(0, 1) +
  labs(y = "", x = "") 

grid.arrange(
  p1,
  p2,
  ncol = 2,
  bottom = textGrob("Distance to apex [#neighbours]", gp = gpar(fontsize = 11))
)

#############################################
### nDivisions at apex
#############################################
load(file = "data/allbut18corresp_d2t_peakexpr_4_expr_filtering.RData")
# load(file = "data/filtered_topexpr4.RData")
# qd = do.call(bind_rows, lapply(results, "[[", "quant.data")) 
# sd = do.call(bind_rows, lapply(results, "[[", "segm.data"))
# dd = do.call(bind_rows, lapply(results, "[[", "division.data"))
dd[which(dd$plant == 18), "d2t"] = NA
sd[which(sd$plant == 18), "d2t"] = NA
qd[which(qd$plant == 18), "d2t"] = NA

nCirc = sd %>% 
  filter(circ < 7, layer == "L1", plant %in% p) %>%
  group_by(plant, circ, t) %>% summarize(nCirc = n())
nDivs = dd %>% 
  filter(circ < 7, layer == "L1", plant %in% p) %>%
  group_by(plant, circ, t) %>% summarize(nDivs = n())
data = left_join(nCirc, nDivs, by = c("plant", "circ", "t"), all = TRUE) %>% ungroup()
data$nDivs[which(is.na(data$nDivs))] = 0
data$nCirc[which(is.na(data$nCirc))] = 0
p1 = data %>% 
  mutate(norm.nDivs = nDivs / nCirc) %>% 
  group_by(plant, circ) %>%
  summarize(mean = mean(norm.nDivs, na.rm = TRUE), se = sd(norm.nDivs, na.rm = TRUE) / sqrt(n())) %>%
  ggplot() +
  geom_point(aes(x = factor(circ), y = mean, colour = factor(plant)), size = 3) +
  geom_errorbar(aes(x = factor(circ), ymin = mean-se, ymax = mean + se,colour = factor(plant)), width = .2, size = 1) +
  facet_wrap(~plant) + 
  labs(y = "Normalised number of divisions", x = "") +
  theme(legend.position = "none")
p2 = data %>% 
  mutate(norm.nDivs = nDivs / nCirc) %>% 
  group_by(circ) %>%
  summarize(mean = mean(norm.nDivs, na.rm = TRUE), se = sd(norm.nDivs, na.rm = TRUE) / sqrt(n())) %>% 
  transform(plant="All") %>%
  ggplot() +
  geom_point(aes(x = factor(circ), y = mean), size = 3, alpha = .7) +
  geom_errorbar(aes(x = factor(circ), ymin = mean-se, ymax = mean + se), width = .2, size = 1, alpha = .7) +
  facet_wrap( ~ plant) + 
  labs(y = "", x = "")
grid.arrange(p1, p2, ncol = 2,bottom = textGrob("Distance to apex [#neighbours]", gp = gpar(fontsize = 11)))

# dd %>% filter(circ < 7, layer == "L1", plant %in% p, !is.na(age)) %>% 
#   group_by(plant, circ, t) %>% 
#   mutate(nInCirc = n()) %>%
#   select(plant, t, age, circ, nInCirc) %>% 
#   group_by(circ) %>%

# ggplot(dd, )
###########################
### Age plots
###########################
data = read.table("data/all_data_with_age_tojose.dat", header = TRUE)
load("data/filtered_topexpr4.RData")

data %>% filter(circ < 7) %>% 
  filter(layer != "L1", t > 48) %>% ggplot(aes(
    x = m.x,
    y = m.y,
    colour = age
  )) + geom_point(na.rm = T, size = 3, alpha = .3) +  geom_point(
    data = dd %>% filter(layer != "L1", t > 48, age > 48),
    aes(x = m.x, y = m.y),
    na.rm = T, colour = 2,
    size  = 3, alpha = .7,
    shape = 25,
    fill = "red"
  ) + facet_grid(plant ~ t) + labs(x = "x", y = "y", shape = "Division event") +   scale_color_continuous(name = "Age") +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "bottom"
  ) 
