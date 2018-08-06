setwd("/home/henrik/compbio/thesis/code/newTreat/")
library(plyr)
library(tidyverse)
source("functions.R")
source("read_data.R")
source("aux.R")
source("../plot.R")
plants = c(1, 2, 4, 13, 15, 18)

files  = list.files("~/Desktop/results/", pattern = "75", full.names = TRUE)
# files  = files[!grepl(files, pattern="top25|top50|top150|top125")]
times  = extract.numbers(files)[1, ]
plants = extract.numbers(files)[2, ]
slices = extract.numbers(files)[3, ]

files = files[order(plants, times, slices)]

data = lapply(files, read.csv)
data = lapply(1:length(data), function(x) {
  cbind(plant      = rep(plants[x], nrow(data[[x]])),
        t          = rep(times[x],  nrow(data[[x]])),
        slices     = rep(slices[x],  nrow(data[[x]])),
        data[[x]])
})
data = do.call(rbind, data) %>% as.tibble() %>% filter(X==0) %>% select(-X)

##############

params = data %>% select(plant, t, slices, para_p1, para_p2)
params = abs(params) #%>%
  # filter(para_p1 < 0.0045 & para_p1 > 0.0015,
         # para_p2 < 0.0045 & para_p2 > 0.0015)
params = params %>% 
  group_by(plant, t, slices) %>% 
  mutate(a = max(para_p1, para_p2, na.rm = T),
         b = min(para_p1, para_p2, na.rm = T)) %>%
  # mutate(a = para_p1, b = para_p2) %>%
  ungroup() %>%
  group_by(plant, t) %>%
  summarize(a = mean(a, na.rm=T), b = mean(b, na.rm = T))

p = ggplot() + 
  geom_point(data = params, aes(x = t, y = a, colour = "blue")) +
  geom_point(data = params, aes(x = t, y = b, colour = "red")) +
  theme_bw() +
  facet_wrap(~plant) + 
  scale_color_discrete(name = "Parameter", labels = c("a", "b")) +
  labs(x="Time [hours]", y = "a, b")
print(p)


combined = left_join(params, nNucl, by = c("t", "plant"))
combined = combined %>% filter(plant == 2)
cat(cor.test(combined$a, combined$count)$p.value, " ",sep = "\t")
cat(cor.test(combined$b, combined$count)$p.value, " ",sep = "\t")
cat(cor.test(combined$a, combined$expr)$p.value, " ",sep = "\t")
cat(cor.test(combined$b, combined$expr)$p.value, "\n",sep = "")

combined = left_join(params, nNucl, by = c("t", "plant"))
combined = combined %>% filter(plant == 4)
cat(cor.test(combined$a, combined$count)$p.value, " ", sep="\t")
cat(cor.test(combined$b, combined$count)$p.value, " ",sep="\t")
cat(cor.test(combined$a, combined$expr)$p.value, " ",sep="\t")
cat(cor.test(combined$b, combined$expr)$p.value, "\n",sep = "")

combined = left_join(params, nNucl, by = c("t", "plant"))
combined = combined %>% filter(plant == 13)
cat(cor.test(combined$a, combined$count)$p.value, " ",sep="\t")
cat(cor.test(combined$b, combined$count)$p.value, " ",sep="\t")
cat(cor.test(combined$a, combined$expr)$p.value, " ",sep="\t")
cat(cor.test(combined$b, combined$expr)$p.value, "\n",sep = "")

combined = left_join(params, nNucl, by = c("t", "plant"))
combined = combined %>% filter(plant == 15)
cat(cor.test(combined$a, combined$count)$p.value, " ",sep="\t")
cat(cor.test(combined$b, combined$count)$p.value, " ",sep="\t")
cat(cor.test(combined$a, combined$expr)$p.value, " ",sep="\t")
cat(cor.test(combined$b, combined$expr)$p.value, "\n",sep = "")

params15 = params %>% filter(plant == 15)
new.data = sapply(1:nrow(params15), function(x) {
  return(unlist(params15[x, "b"]) * (-10:10) ** 2)
}) %>% t() %>% cbind(t=params15$t, .)
colnames(new.data)[-1] = -10:10
melted = new.data %>% as.tibble() %>%  melt(id.vars = c("t"))

ggplot(melted, aes(x=as.integer(variable), y = -value, group = t, colour =t)) + 
  geom_line(size = 1)  +
  labs(x = "", y  = "") +
  scale_color_continuous(name = "Time") +
  theme(axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()
        )
  

ggplot(params %>% filter(plant == 15), 
       aes(x = -10:10, y = a * (-10:10) ** 2)) + 
  geom_line()
