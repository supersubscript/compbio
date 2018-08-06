###############################################################################
### GENERAL IMPORT AND PLOT SETTINGS
###############################################################################
setwd("/home/henrik/compbio/thesis/code/")

# Load packages
packages = c("reshape", "scales", "stringr", "plyr", "parallel", "limma", 
             "affy", "RColorBrewer", "ggplot2", "grid", "gridExtra")
try(lapply(packages, library, character.only = TRUE), silent = TRUE)

# Set general colour palette
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 2
par(lwd = lw.s, cex = 1.5, ps = 15)

###############################################################################
### AUXILLIARY FUNCTIONS
###############################################################################
order.na = function(frame.list) {
  order(unlist(lapply(frame.list, function(x)
    sum(is.na(
      unlist(x)
    )))))
}

distance = function(v1, v2) {
  sqrt(sum((as.vector(v1) - as.vector(v2)) ** 2))
}

order.by.columns = function(data) data[do.call(order, lapply(1:NCOL(data), function(ii) data[, ii])), ]
order.by.rows    = function(data) data[do.call(order, lapply(1:NROW(data), function(ii) data[ii, ])), ]

init.counter = function() {
  x = 0
  function() {
    x <<- x + 1
    x
  }
}

# Is not NA
nna = function(data){!is.na(unlist(data))}

### Extract numbers from string
extract.numbers = function(strs) {
  return(unname(sapply(strs, function(str) {
    matches = gregexpr('[0-9]+', str)
    sapply(regmatches(str, matches), as.numeric)
  })))
}

