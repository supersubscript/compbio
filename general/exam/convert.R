#!/usr/bin/env Rscript
setwd("~/questions/")

args = commandArgs(trailingOnly=TRUE)
data = readLines(args[1])

empty.lines = grepl('^-|^Lecture', data)
data = data[!empty.lines]
empty.lines = grepl('^\\s*$', data)
data[empty.lines] = "WHATWHAT"
#data = data[!empty.lines]
questions = !grepl('^ |WHATWHAT', data) # get question lines

data[questions] = sapply(data[questions], function(x)  paste0(x, "//")) # Append //
questions = !grepl('^ ', data) # get question lines
data = trimws(data, which = c("both", "left", "right")) # trim whitespace
data[!questions] = unlist(sapply(data[!questions], function(x) paste0(x, "__"))) # sub-lines by __
data = gsub("\n","", data)

out.datav = paste(data, collapse = " ") # collapse
out.datav = strsplit(out.datav, "WHATWHAT")[[1]]
out.datav = trimws(out.datav, which = c("both", "left", "right")) # trim whitespace
out.datav = gsub("__$", "", out.datav)
out.datav = unname(sapply(out.datav, function(x) gsub("__ ", "__", x)))

outfile = gsub("\\.txt", "_flash\\.txt", args[1])
write(out.datav, file = outfile, sep ="\n")
