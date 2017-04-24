setwd("~/compbio/text/scientific_programming/assignment_2/")
.libPaths(c(.libPaths(), '/local/data/public/hpa22/R/x86_64-pc-linux-gnu-library/3.1/'))

### Read in data, avoiding unnecessary columns
files = list.files(path = "/local/data/public/sje30/weather/2016/daily-text", full.name = TRUE)
files = files[which(file.info(files)$size >= 1000)]
weather.data        = lapply(files, function(x) read.table(x, header = FALSE, col.names = c("Time", "Temp", "Humid", "DewPt", "Press", "WindSp", "WindDr", "Sun", "Rain", "Start", "MxWSpd"), fill = TRUE)) 
names(weather.data) = unname(gsub("_", "-", sapply(files, function(x) tail(strsplit(x, "/")[[1]], n = 1))))
weather.data        = lapply(weather.data, function(x) x[, c("Time", "Temp", "Sun", "Rain")])
filtered.out        = c()
par(lwd = 3, cex = 1.5)

### Filter away days with too few data points
round.times = function(x){
  times = strptime(x[, "Time"], format = "%H:%M")
  mins  = round(times$min / 30) * 30
  mins[which(mins == 60)] = 59
  times$min  = mins
  x[,"Time"] = paste(times$hour, ifelse(times$min == 0, "00", times$min), sep = ":")
  return(x)
}
weather.data = lapply(weather.data, round.times) 

### Remove row and timestamp duplicates. Only keep last (most recent) timestamp.
weather.data = lapply(weather.data, function(x) x[order(as.Date(x[,1], format = "%H:%M")), ]) 
weather.data = lapply(weather.data, function(x) x[!duplicated(x), ])
weather.data = lapply(weather.data, function(x) x[!rev(duplicated(rev(x[,1]))), ]) 

### Remove files with too little data
large.enough = sapply(weather.data, function(x) nrow(x) == 49)
filtered.out = append(filtered.out, weather.data[-large.enough])
weather.data = weather.data[large.enough]

### Filter out days which have been corrected 
corrected    = which(sapply(names(weather.data), function(x) paste0(x,"~") %in% names(weather.data)))
filtered.out = append(filtered.out, corrected)
weather.data = weather.data[-corrected]

###### 2.
times = strptime(weather.data[["2012-12-25"]][, "Time"], format = "%H:%M")
plot(times, weather.data[["2012-12-25"]][, "Temp"], bty='n', ylab = "Temperature", xlab = "Time", type= 'o', pch = 16)

###### 3.
### Remove things that contain non-numerics
remove.faulty = function(days, column)
{
  a = days[names(which(sapply(days, function(x) !is.numeric(x[, column]))))]
  a = append(a, days[names(which(sapply(weather.data, function(x) anyNA(as.numeric(x[, column])))))])
  days[which(names(days) %in% names(a))] = NULL
  return(days)  
}
temperature.days = weather.data
temperature.days = remove.faulty(temperature.days, "Temp")

### Remove weird temperatures
temperature.days = temperature.days[-which(sapply(temperature.days, function(x) max(x[,"Temp"])) > 35)]
temperature.days = temperature.days[-which(sapply(temperature.days, function(x) min(x[,"Temp"])) < -25)]

### Find mean temperature over year
mean.temperatures = sapply(temperature.days, function(x) mean(as.numeric(x[, "Temp"])))
max.name = names(which(mean.temperatures == max(mean.temperatures)))
min.name = names(which(mean.temperatures == min(mean.temperatures)))
max.mean = max(mean.temperatures)
min.mean = min(mean.temperatures)

plot(as.Date(names(mean.temperatures), "%Y-%m-%d"), mean.temperatures, cex = .2, xlab = "Year", ylab = expression(paste("Temperature (",degree,"C)")), ylim=c(-10,30), xlim=c(as.numeric(as.Date("1995-01-01")), as.numeric(as.Date("2016-12-31"))))
points(as.Date(min.name), min.mean, cex = 1.5, col = 'slateblue3', pch = 16)
points(as.Date(max.name), max.mean, cex = 1.5, col = 'violetred3', pch = 16)

#### 4. Fit mean temperature over whole data set
mean.mean.temps  = rep(0, 366)
mean.mean.counts = rep(0, 366)
names(mean.mean.temps)  = seq.Date(as.Date("2016-01-01"), as.Date("2016-12-31"), "days")
names(mean.mean.counts) = seq.Date(as.Date("2016-01-01"), as.Date("2016-12-31"), "days")
### Calculate the means
for (ii in names(temperature.days))
{
  date = as.character(as.Date(format(as.Date(ii, format = "%Y-%m-%d"), format = "%m-%d"), format = "%m-%d"))
  mean.mean.temps[date]  = mean.mean.temps[date]  + mean(temperature.days[as.character(ii)][[1]][,"Temp"])
  mean.mean.counts[date] = mean.mean.counts[date] + 1
}  
mean.mean.temps = mean.mean.temps / mean.mean.counts

### Per year
year.averages = function(x, FUN)
{
  mean.temps  = rep(NA,366)
  names(mean.temps)  = seq.Date(as.Date("2016-01-01"), as.Date("2016-12-31"), "days")
    
  ### Loop through days in year
  for(ii in names(x))
  {
    date = as.character(as.Date(format(as.Date(ii, format = "%Y-%m-%d"), format = "%m-%d"), format = "%m-%d"))
    mean.temps[date] = FUN(as.numeric(temperature.days[as.character(ii)][[1]][, "Temp"]))
  }
  return (mean.temps[which(!is.na(mean.temps))])
}


### Find best year
min.val = 99999999
max.val  = 0
min.year = ""
max.year = ""

for(ii in seq(1996, 2016))
{
  days          = temperature.days[grep(ii, names(temperature.days))]
  mean.temps    = year.averages(days, mean)
  data.for      = which(!is.na(mean.temps))
  sq.difference = sum((mean.temps[data.for] - mean.mean.temps[data.for])^2) #sum((mean.temps[data.for] - mean.mean.temps[data.for]))^2
  
  if(sq.difference > max.val)
  {
    max.val  = sq.difference
    max.year = ii
  }
  if(sq.difference < min.val)
  {
    min.val  = sq.difference
    min.year = ii
  }
}
# Most similar
days                   = temperature.days[grep(min.year, names(temperature.days))]
most.alike.mean.temps  = year.averages(days, mean)
most.alike.mean.temps  = most.alike.mean.temps[which(!is.na(most.alike.mean.temps))]

# Least similar
days                   = temperature.days[grep(max.year, names(temperature.days))]
least.alike.mean.temps = year.averages(days, mean)
least.alike.mean.temps = least.alike.mean.temps[which(!is.na(least.alike.mean.temps))]

plot(as.Date(names(mean.mean.temps), format = "%Y-%m-%d"), mean.mean.temps, type = 'l', ylab = expression(paste("Mean mean temperature (", degree, "C)")), xlab = "Month", ylim = c(0,25), cex = 3, lwd = 3, lty = 2)
lines(as.Date(names(least.alike.mean.temps), format = "%Y-%m-%d"), least.alike.mean.temps,col= 'slateblue3', lwd = 2)
lines(as.Date(names(most.alike.mean.temps),  format = "%Y-%m-%d"), most.alike.mean.temps, col= 'violetred3', lwd = 2)
legend("topleft", c("Monthly average", paste0("Least alike: ", max.year), paste0("Most alike:  ", min.year)), lty = c(2,1,1), col = c(1, 'slateblue3', 'violetred3'), lwd=c(3,2,2), cex = .8)

### Take out dates of choice
removal.dates = function(x, from, to)
{
  remove = c()
  for (ii in 1:length(to))
  {
    remove = append(remove, seq.Date(as.Date(from[ii], "%Y-%m-%d"), as.Date(to[ii], "%Y-%m-%d"), "days"))
  }
  return (remove)
}

rain.from = c(
  "2001-07-01", "2004-07-09", "2005-01-21", "2005-06-28", "2005-10-20", "2006-07-27", 
  "2006-10-07", "2007-10-09", "2007-11-22", "2008-07-08", "2008-08-12", "2008-08-19",
  "2008-08-03", "2008-08-04", "2008-10-25", "2009-04-10", "2010-02-06", "2010-08-23", 
  "2011-02-26", "2011-08-19", "2011-09-19", "2011-02-26", "2011-09-06"
)
rain.to   = c(
  "2001-08-31", "2004-07-18", "2005-03-02", "2005-07-01", "2005-11-07", "2006-07-31",
  "2006-10-09", "2007-10-16", "2007-11-22", "2008-07-08", "2008-08-12", "2008-08-27",
  "2008-08-03", "2008-08-04", "2008-11-04", "2009-04-13", "2010-02-28", "2010-08-23",
  "2011-02-26", "2011-08-23", "2020-12-31", "2011-02-26", "2011-09-06"
)
rain.malfunction.dates = removal.dates(weather.data, rain.from, rain.to)

### Remove unwanted rain days
rainy.days = weather.data
rainy.days = remove.faulty(rainy.days, "Rain")
rainy.days[which(names(rainy.days) %in% as.character(rain.malfunction.dates))] = NULL
rain.and.temperature.data  = temperature.days
rain.and.temperature.data  = rain.and.temperature.data[which(names(rain.and.temperature.data) %in% names(rainy.days))]

temperatures = sapply(sapply(rain.and.temperature.data, "[", "Temp"), mean)
rain         = sapply(sapply(rain.and.temperature.data, "[", "Rain"), max)
plot(temperatures, rain, bty ='n', type = 'p', xlab = expression(paste("Temperature (", degree, "C)")), ylab = "Rain (mm)", lwd = .2, pch = 16, cex=.3)
cat("Temperature-Rain correlation: ", cor(temperatures, rain), "\n")

######## 6.
### Remove faulty and unwanted sun days
sunny.days = weather.data
sunny.days = remove.faulty(sunny.days, "Sun")
sun.from   = c("2009-04-03", "2007-06-19", "2007-06-04", "2007-08-04", "2007-12-31")
sun.to     = c("2009-04-03", "2007-06-19", "2007-06-04", "2007-08-06", "2008-01-02")
sunshine.malfunction.dates = removal.dates(weather.data, rain.from, rain.to)

### > 16 h sun? Remove.
sunny.days = weather.data[lapply(sunny.days, function(x) max(x[, "Sun"])) <= 16]

### Wanted statistics
sunlight.and.temperature.data = sunny.days[which(names(sunny.days) %in% names(temperature.days))]
sun          = sapply(sapply(sunlight.and.temperature.data, "[", "Sun"),  max)
temperatures = sapply(sapply(sunlight.and.temperature.data, "[", "Temp"), mean)
plot(temperatures, sun, bty = 'n', type ='p', xlab = expression(paste("Temperature (",degree,"C)")), ylab = "Sunlight (hours)", cex = .3, pch = 16, ylim = c(0,16))
cat("Temperature-Sunlight correlation: ", cor(temperatures, sun), "\n")

### 7.

mins   =  vector(mode = "double", length = 12)
maxes  =  vector(mode = "double", length = 12)
means  =  vector(mode = "double", length = 12)
counts =  vector(mode = "double", length = 12)

### Mind mean, min, max over months
for(ii in 1995:2016)
{
  year = temperature.days[grep(ii, names(temperature.days))]
  for(jj in 1:12)
  {
    month = year[as.numeric(lapply(names(year), function(x) strsplit(x, "-")[[1]][2])) == jj]
    if(length(month) > 0)
    {
      mins[jj]   = mins[jj]   + mean(sapply(sapply(month, "[", "Temp"), min))
      maxes[jj]  = maxes[jj]  + mean(sapply(sapply(month, "[", "Temp"), max))
      means[jj]  = means[jj]  + mean(sapply(sapply(month, "[", "Temp"), mean))
      counts[jj] = counts[jj] + 1
    }
  }
}
mins  = mins  / counts
maxes = maxes / counts
means = means / counts

plot(seq.Date(as.Date("2016-01-01"), as.Date("2016-12-31"), "months"),  maxes, type='o', col = 'violetred3', ylim = c(-10,30), xlab = "Month", ylab = expression(paste("Temperature (", degree, "C)")), pch = 16)
grid(nx = 12, ny = 4, col = "gray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
lines(seq.Date(as.Date("2016-01-01"), as.Date("2016-12-31"), "months"), means, type='o', pch = 16)
lines(seq.Date(as.Date("2016-01-01"), as.Date("2016-12-31"), "months"), mins,  type='o', col = 'slateblue3', pch = 16)

### 8.

### Find longest consecutive interval with rain > wet.day.threshold
wet.day.threshold = .2
from = ""
to   = ""
consequtive = 0
max.consequtive = 0
best.from = ""

### First day rainy? 
if(max(as.numeric(rainy.days[[1]][, "Rain"])) > wet.day.threshold)
{
  consequtive = 1
}
for(ii in 2:length(rainy.days))
{
  # are the dates consequtive?
  # is this a rainy day?
  this.rainy = max(as.numeric(rainy.days[[ii]][, "Rain"])) > wet.day.threshold 
  if((as.Date(names(rainy.days[ii])) !=  (as.Date(names(rainy.days[ii - 1])) + 1)) | !this.rainy)
  {
    consequtive = 0
  }
  
  # is this one rainy?
  if(this.rainy)
  {
    if(consequtive == 0)
    {
      from = names(rainy.days[ii])
    }
    consequtive = consequtive + 1
    if(consequtive > max.consequtive)
    {
      max.consequtive = consequtive
      to = names(rainy.days[ii])
      best.from = from
    }
  }
}
### Which dates had longest streak? 
streak.winner = seq.Date(as.Date(best.from), as.Date(to), "days")
print(streak.winner)
print(max.consequtive)

### What have we taken out? 
filtered.rain = weather.data[which(!names(weather.data) %in% names(rainy.days))]
filtered.sun  = weather.data[which(!names(weather.data) %in% names(sunny.days))]
filtered.temp = weather.data[which(!names(weather.data) %in% names(temperature.days))]
filtered.all  = weather.data[which(!names(weather.data) %in% c(names(rainy.days), names(sunny.days), names(temperature.days)))]

complete.sets = weather.data[which(!names(weather.data) %in% c(names(filtered.rain), names(filtered.sun), names(filtered.temp))) ]
print(length(complete.sets))
print(length(complete.sets) / length(seq.Date(as.Date("1995-07-01"), as.Date("2016-10-25"), "days")))

### Hi, Stephen!