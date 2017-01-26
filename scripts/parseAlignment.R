#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

inFile = readLines(args[1])
trim.leading = function (x)  gsub("^\\s+", "", x)

### Do the deed.
queries = grep("Query: ", inFile)
pairs = paste(trim.leading(inFile[queries]), "\n", trim.leading(inFile[queries + 1]), "\n\n", sep = "")

if(length(args) < 2)
{
   write(pairs, "parsing_out.txt")  
}else
{
   write(pairs, args[2])  
}


