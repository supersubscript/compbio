data = cld[[1]]
# data = do.call(c, cld[[1]])

new.data = lapply(data, function(x)
  lapply(x,
         function(y)
           y %>% filter(circ > 1)))

par(mfrow = c(1, 1), mar =c(4,4,4,4))
plot(NA, ylim = c(0, 160), xlim = c(0, 84), xlab = "t", ylab = "E")
for(ii in 1:length(new.data))
  for (jj in 1:length(new.data[[ii]])) {
    if (nrow(new.data[[ii]][[jj]]) > 8)
      lines(new.data[[ii]][[jj]]$t - min(new.data[[ii]][[jj]]$t, na.rm=T), new.data[[ii]][[jj]]$expr, col = sample(colors(10)))
  }
lines(all.data$t, predict(fit), col = "red", lwd = 10)

all.lines = unlist(new.data, recursive = F)
all.lines = unique(all.lines[sapply(all.lines, nrow) > 8])
all.lines = lapply(all.lines, function(x) {
  x$t = x$t - min(x$t, na.rm=T)
  x
})

all.data = do.call(rbind, all.lines)
plot(all.data$t, all.data$expr)

library(MASS)
all.data = all.data %>% arrange(t)
all.data = all.data %>% filter(!is.na(t) & !is.na(expr))
fit = rlm(all.data$expr~all.data$t)
lines(all.data$t, predict(fit), col = "red", lwd = 10)

