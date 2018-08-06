
get.plant.data = function(plant.no, map.quant = FALSE, only.l1 = FALSE, only.l2 = FALSE, quant.missing = -1, mapping.missing = -1) {
  # Load the necessary packages
  source("aux.R",          local = TRUE)
  source("read_data.R",    local = TRUE)
  source("process_data.R", local = TRUE)
  packages = c("reshape", "scales", "stringr", "plyr", "parallel", "limma", "affy")
  try(lapply(packages, library, character.only = TRUE), silent = TRUE)
  
  # Read in data
  all.data = read.data(
    plant.no,
    map.quant       = map.quant,
    only.l1         = only.l1,
    only.l2         = only.l2,
    quant.missing   = quant.missing,
    mapping.missing = mapping.missing
  )
  
  mapping.data = all.data$mapping.data
  quant.data   = all.data$quant.data
  timepoints   = all.data$timepoints
  
  # Extract lineages
  sublines = get.sublines(mapping.data, timepoints)
  lineages = get.lineages(sublines)
  
  # Get the data for all the sublines 
  lineage.sublines.data = list()
  if (plant.no != 1)
    lineage.sublines.data = lapply(lineages,
                                   get.lineage.sublines.data,
                                   quant.data = quant.data,
                                   timepoints = timepoints)
  
  # ... and for the whole lineages 
  # cell.lines.data = lapply(lineage.sublines.data, collapse.lineage)
  
  # Retrieve table with all division events and their corresponding information
  all.data = get.division.events(
    mapping.data,
    quant.data,
    sublines,
    plant.no)
  return(list(all.data = all.data, lineage.sublines.data = lineage.sublines.data))
}