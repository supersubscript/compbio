###############################################################################
### GENERAL IMPORT AND PLOT SETTINGS
# Load packages
packages = c(
  "plyr",
  "reshape2",
  "scales",
  "stringr",
  "dplyr",
  "parallel",
  "limma",
  "affy",
  "RColorBrewer",
  "ggplot2",
  "grid",
  "gridExtra",
  "tidyverse",
  "snow",
  "doParallel"
)
try(lapply(packages, library, character.only = TRUE), silent = TRUE)

# Set general colour palette
palette(brewer.pal(n = 8, name = "Set1"))
lw.s = 2
par(lwd = lw.s, cex = 1.5, ps = 15)

#' gg_color_hue
#'
#' @param n Numbers of colours
#'
#' @return Array with rgb colours corresponding to the base ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

###############################################################################
### AUXILLIARY FUNCTIONS
###############################################################################

#' Order NAa
#'
#' @param frame.list 
#'
#' @return
#' @export
#'
#' @examples
order.na = function(frame.list) {
  order(unlist(lapply(frame.list, function(x)
    sum(is.na(
      unlist(x)
    )))))
}

distance = function(v1, v2) {
  sqrt(sum((as.vector(v1) - as.vector(v2)) ** 2))
}

order.by.columns = function(data)
  data[do.call(order, lapply(1:NCOL(data), function(ii)
    data[, ii])),]
order.by.rows    = function(data)
  data[do.call(order, lapply(1:NROW(data), function(ii)
    data[ii,])),]

init.counter = function() {
  x = 0
  function() {
    x <<- x + 1
    x
  }
}

# Is not NA
nna = function(data) {
  !is.na(unlist(data))
}

### Extract numbers from string
extract.numbers = function(strs) {
  return(unname(sapply(strs, function(str) {
    matches = gregexpr('[0-9]+', str)
    sapply(regmatches(str, matches), as.numeric)
  })))
}

detach.all.packages <- function() {
  basic.packages =
    c(
      "package:stats",
      "package:graphics",
      "package:grDevices",
      "package:utils",
      "package:datasets",
      "package:methods",
      "package:base"
    )
  package.list <-
    search()[ifelse(unlist(gregexpr("package:", search())) == 1, TRUE, FALSE)]
  package.list <- setdiff(package.list, basic.packages)
  if (length(package.list) > 0)
    for (package in package.list)
      detach(package, character.only = TRUE)
}