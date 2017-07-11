setwd("/home/henrik/compbio/thesis/code/newTreat/")
source("functions.R")
plants = c(1, 2, 4, 13, 15, 18)

###############################################################################
get.plant.data = function(plant) {
  setwd("/home/henrik/compbio/thesis/code/newTreat/")
  source("functions.R", local = TRUE)  
  source("../aux.R",    local = TRUE)
  library(tidyverse)
  
  q.missing = list(seq(0, 76, 4), c(), 24, c(), c(12), c()); names(q.missing) = plants
  m.missing = list(c(), c(), c(), c(), c(), c(40)); names(m.missing) = plants
  base.path  = "/home/henrik/compbio/thesis"
  corr.path  = "/data/correspondences/plant"
  data.path  = "/data/newData/plant"
  code.path  = "/code"
  quant.path = "/data/clv3_complete/plant"
  
  # List files
  quant.files   = list.files(paste0(base.path, quant.path, plant, "/Results"),           pattern = ".txt", full.names = TRUE)
  mapping.files = list.files(paste0(base.path, data.path,  plant, "/tracking_data"),     pattern = ".pkl", full.names = TRUE)
  segm.files    = list.files(paste0(base.path, data.path,  plant, "/segmentation_data"), pattern = ".pkl", full.names = TRUE)
  corr.files    = list.files(paste0(base.path, corr.path,  plant, ""),                   pattern = "corr", full.names = TRUE)
  
  # Put in order
  quant.files   = if(length(quant.files) > 0) quant.files[order(extract.numbers(quant.files)[3, ])]
  mapping.files = mapping.files[order(extract.numbers(mapping.files)[2, ])]
  segm.files    = segm.files[order(extract.numbers(segm.files)[2, ])]
  corr.files    = if(length(corr.files) > 0) corr.files[order(extract.numbers(corr.files)[3, ])]
  timepoints    = seq(0, max(if(length(quant.files) > 0) extract.numbers(quant.files)[3, ] else extract.numbers(mapping.files)[3, ]), by = 4)
  
  ###############################################################################
  ###############################################################################
  ###############################################################################
  
  library(snow)
  library(doParallel)
  cl = makeCluster(3)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(tidyverse))
  # clusterEvalQ(cl, library(tidyverse))
  clusterExport(cl, list = ls(), envir = environment())
  mapping.data = parLapply(cl, mapping.files, get.mapping.data)
  quant.data   = parLapply(cl, quant.files,   get.quant.data)
  corr.data    = parLapply(cl, corr.files,    get.corresp)
  volume.data  = parSapply(cl, segm.files,    get.volume.data)
  centers.data = parSapply(cl, segm.files,    get.centers.data)
  neigh.data   = parLapply(cl, segm.files,    get.neigh.data)
  segm.data    = rbind(volume.data, centers.data, neigh.data)
  
  mapping.data = fix.mapping(mapping.data, m.missing[[as.character(plant)]], timepoints = timepoints)
  segm.data    = fix.segm(segm.data,       m.missing[[as.character(plant)]], timepoints = timepoints)
  quant.data   = fix.quant(quant.data,     q.missing[[as.character(plant)]], timepoints = timepoints)
  corr.data    = fix.corresp(corr.data,    q.missing[[as.character(plant)]], timepoints = timepoints)
  
  sublines     = get.sublines(mapping.data, timepoints)
  lineages     = get.lineages(sublines)
  
  quant.data = parLapply(cl, 1:length(quant.data), quant2segm.mapping, quant.data = quant.data, corr.data = corr.data)
  names(quant.data) = timepoints
  
  # TODO: get.l1 and get.l2 should be aware of missing data (in mapping)
  # TODO: remove duplicates in quant.data? Prolly not.
  # quant.data[[4]] in plant 4 looks super weird. There are like 5 nuclei mapping to membrane no. 92
  segm.data = add.layers(segm.data, segm.files, neigh.data, cl = cl)
  
  # Get division events
  no.division.events = sapply(mapping.data, function(mapping.event)
    mapping.event$parent %>%
      duplicated         %>%
      sum)               %>%
    sum
  
  divisions = parLapply(cl, mapping.data, function(mapping.event) {
    mapping.event$parent %>% duplicated %>% which
  })
  
  division.data = foreach(ii = 1:(length(timepoints) - 1), .combine = rbind) %:%
    foreach(jj = divisions[[ii]], .combine = rbind) %dopar% 
    get.division.event(ii, jj, sublines, mapping.data, segm.data, quant.data, 
                       timepoints, lineages, cl)
  
  division.data[division.data$m.vol == 0, "m.vol"] = NA
  division.data[division.data$n.vol == 0, "n.vol"] = NA
  
  # Get cell line data
  cell.lines = get.cell.lineages(sublines, lineages)
  
  lineage.sublines.data = list()
  lineage.sublines.data = parLapply(cl, 
                                    cell.lines,
                                    get.lineage.sublines.data,
                                    quant.data = quant.data,
                                    timepoints = timepoints
  )
  cell.lines.data = parLapply(cl, lineage.sublines.data, collapse.lineage) 
  
  # Sort by number of members?
  lineage.sublines.data = lineage.sublines.data[parSapply(cl, cell.lines.data, function (x)
    x[, "d2t"] %>% nna %>% sum) %>% unlist %>% order(decreasing = TRUE)]
  
  list(
    mapping.data          = mapping.data,
    quant.data            = quant.data,
    segm.data             = segm.data,
    division.data         = division.data,
    cell.lines            = cell.lines,
    lineage.sublines.data = lineage.sublines.data,
    cell.lines.data       = cell.lines.data,
    timepoints            = timepoints
  ) %>% return
}

