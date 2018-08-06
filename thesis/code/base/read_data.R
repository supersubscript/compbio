
get.plant.data = function(plant, remove.quant.dups = TRUE, map.quant = TRUE, d2t.method = "peak.expr", d2t.n = 4) {
  # remove.quant.dups = TRUE; map.quant = TRUE; d2t.method = "w.expr"; d2t.n = 4
  source("functions.R", local = TRUE)  
  source("aux.R",       local = TRUE)
  plants = c(1, 2, 4, 13, 15, 18)
  cat("Running plant ", plant, "...\n") 
  cat("\tSetting up files\n")
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
  
  cl = makeCluster(3)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(tidyverse))
  clusterExport(cl, list = ls(), envir = environment())
  cat("\tRetrieving bulk data\n") 
  mapping.data = parLapply(cl, mapping.files, get.mapping.data, remove.score = TRUE)
  quant.data   = parLapply(cl, quant.files,   get.quant.data)
  corr.data    = parLapply(cl, corr.files,    get.corresp)
  volume.data  = parSapply(cl, segm.files,    get.volume.data)
  centers.data = parSapply(cl, segm.files,    get.centers.data)
  neigh.data   = parLapply(cl, segm.files,    get.neigh.data)
  segm.data    = rbind(volume.data, centers.data, neigh.data)
  
  cat("\tFixing bulk data\n")
  mapping.data = fix.mapping(mapping.data,  m.missing[[as.character(plant)]], timepoints = timepoints)
  segm.data    = fix.segm(segm.data,        m.missing[[as.character(plant)]], timepoints = timepoints)
  neigh.data   = fix.neigh.data(neigh.data, m.missing[[as.character(plant)]], timepoints = timepoints)
  quant.data   = fix.quant(quant.data,      q.missing[[as.character(plant)]], timepoints = timepoints, plant = plant)
  corr.data    = fix.corresp(corr.data,     q.missing[[as.character(plant)]], timepoints = timepoints)
  
  cat("\tLoading sublines / lineages\n")
  sublines     = get.sublines(mapping.data, timepoints)
  lineages     = get.lineages(sublines)
  
  cat("\tMapping nuclear ID's\n")
  # TODO: Make this neat
  if ((length(corr.files) == 0 | !map.quant) & length(quant.files) > 0)  {
    quant.data = parLapply(cl, quant.data, function(qd) { qd$id = qd$id + 10000; qd })
  } else if (map.quant) {
    quant.data = parLapply(cl, 1:length(quant.data), quant2segm.mapping, quant.data = quant.data, corr.data = corr.data)
  } else {
    cat("\t\tNo nuclear quantification files found. Mapping unsuccessful.")
  }
  
  ### Compress multiple nuclei mapping to the same cell membrane
  if (remove.quant.dups & length(quant.files) != 0) {
    cat("\tRemoving duplicates in nuclear data\n")
    quant.data = parLapply(cl, quant.data, consolidate.duplicates.quant)
  }
  
  # quant.data[[4]] in plant 4 looks super weird. There are like 5 nuclei mapping to.dat membrane no. 92
  cat("\tAdding periclinal layers\n")
  segm.data  = add.layers.segm(segm.data, segm.files, neigh.data, cl = cl, missing = m.missing[[as.character(plant)]])
  quant.data = add.layers.quant(quant.data, segm.data, cl)
  
  ### Add distance to the top cell
  cat("\tAdding d2t measure\n")
  quant.data = parLapply(cl, quant.data, add.d2t, method = d2t.method, n = d2t.n)

  
  ########################
  ### Get division events
  ########################
  cat("\tCalculating division events\n")
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
                       timepoints, lineages, cl, plant)
  
  # TODO: Make this nice
  segm.data = rbind(plant = sapply(1:ncol(segm.data), function(x)
    rep(plant, length(segm.data[[1, x]]))), segm.data)
  
  #############################
  ### Make data proper tibbles
  #############################
  # segm.data
  segm.data.n = colnames(segm.data)
  segm.data = lapply(1:ncol(segm.data), function(x)
    as.tibble(segm.data[, x]))
  segm.data = lapply(1:length(segm.data), function(x) {
    add_column(segm.data[[x]],
               t = rep(as.integer(segm.data.n[x]), nrow(segm.data[[x]])),
               .before = 1)
  }) %>% do.call(rbind, .)
  # quant.data
  quant.data = do.call(bind_rows, quant.data)
  if(tolower(d2t.n) == "all")
    d2t.n = nrow(quant.data)
  
  ### Filter a little
  cat("\tFiltering data\n")
  MVOL_MIN      = 0
  NVOL_MIN      = 0
  
  MVOL_MAX = mean(segm.data$m.vol, na.rm = TRUE) + sd(segm.data$m.vol, na.rm = T) * 3
  NVOL_MAX = mean(quant.data %>% filter(layer == "L1") %>% .$n.vol, na.rm = TRUE) + 
               sd(quant.data %>% filter(layer == "L1") %>% .$n.vol, na.rm = TRUE) * 5 
  D2T_MAX = 
    mean(quant.data %>%
           filter(layer == "L1") %>%
           .$d2t , na.rm = TRUE) +
    sd(quant.data %>% 
         filter(layer == "L1") %>% 
         .$d2t, na.rm = T) * 3
  
  # Divisions
  if(length(mapping.files) != 0) {
    division.data[which(division.data$m.vol        < MVOL_MIN), "m.vol"] = NA
    division.data[which(division.data$m.vol        > MVOL_MAX), "m.vol"] = NA
    division.data[which(division.data$n.vol        > NVOL_MAX), "n.vol"] = NA
    division.data[which(division.data$d2t          > D2T_MAX),  "d2t"]   = NA
    division.data[which(division.data$d2t          > D2T_MAX),  "expr"]  = NA
    division.data[which(division.data$d2t          > D2T_MAX),  "age"]   = NA
    division.data[which(division.data$d2t          > D2T_MAX),  "n.vol"] = NA
    division.data[which(division.data$d2t          > D2T_MAX),  "m.vol"] = NA
    division.data[which(division.data$child1.m.vol < MVOL_MIN), "m.vol"] = NA
    division.data[which(division.data$child1.m.vol > MVOL_MAX), "m.vol"] = NA
    division.data[which(division.data$child1.n.vol > NVOL_MAX), "n.vol"] = NA
    division.data[which(division.data$child1.d2t   > D2T_MAX),  "d2t"]   = NA
    division.data[which(division.data$child2.m.vol < MVOL_MIN), "m.vol"] = NA
    division.data[which(division.data$child2.m.vol > MVOL_MAX), "m.vol"] = NA
    division.data[which(division.data$child2.n.vol > NVOL_MAX), "n.vol"] = NA
    division.data[which(division.data$child2.d2t   > D2T_MAX),  "d2t"]   = NA
  }
  
  ### Append the segmented info we're lacking to the respective tables
  cat("\tMerging membrane and nuclear data\n")
  quant.data = quant.data %>% left_join(segm.data,
                                        by = intersect(colnames(quant.data), colnames(segm.data)),
                                        all.x = TRUE)
  
  segm.data  = segm.data %>% left_join(quant.data,
                                       by = intersect(colnames(quant.data), colnames(segm.data))[-9],
                                       all.x = TRUE)
  segm.data = segm.data %>% select(-one_of("neighs.y"))
  segm.data = dplyr::rename(segm.data, neighs = neighs.x)
  
  # Membranes
  if(length(segm.files) != 0) {
    segm.data$m.vol[which(segm.data$m.vol > MVOL_MAX | segm.data$m.vol < MVOL_MIN)] = NA
    segm.data$n.vol[which(segm.data$n.vol > NVOL_MAX | segm.data$n.vol < NVOL_MIN)] = NA
    segm.data$m.vol[which(segm.data$n.vol > MVOL_MAX | segm.data$m.vol < MVOL_MIN)] = NA
    segm.data$m.vol[which(segm.data$d2t   > D2T_MAX)] = NA
    segm.data$d2t[which(segm.data$d2t   > D2T_MAX)] = NA
    segm.data$n.vol[which(segm.data$d2t > D2T_MAX)] = NA
    segm.data$expr[which(segm.data$d2t  > D2T_MAX)] = NA
  }
  
  # tmp = quant.data
  # Nuclei
  if(length(quant.files) != 0) {
    quant.data$n.vol[which(quant.data$n.vol > NVOL_MAX | quant.data$n.vol < NVOL_MIN)] = NA
    quant.data$m.vol[which(quant.data$n.vol > MVOL_MAX | quant.data$m.vol < MVOL_MIN)] = NA
    quant.data$m.vol[which(quant.data$d2t   > D2T_MAX)] = NA
    quant.data$d2t[which(quant.data$d2t     > D2T_MAX)] = NA
    quant.data$n.vol[which(quant.data$d2t   > D2T_MAX)] = NA
    quant.data$expr[which(quant.data$d2t    > D2T_MAX)] = NA
  }
  
  cat("\tAdding circ information\n")
  segm.data = segm.data %>% 
    group_by(plant, t) %>% 
    mutate(circ = ifelse(d2t %in% sort(d2t)[1:d2t.n] & layer == "L1", integer(1), NA_integer_))
  
  prev.unmapped = vector(mode="numeric", length = 9000000)
  add.circles = function(data, prev.unmapped) {
    unmapped = which(is.na(data$circ))
    if (length(unmapped) == length(prev.unmapped))
      return(data)
    cat("\t\t", length(unmapped), " still unmapped\n")
    max.circ = max(data$circ, na.rm = TRUE)
    
    tplant = data %>%
      filter(circ == max.circ) %>%
      group_by(t, plant) %>%
      summarize(n = list(neighs))
    
    new = right_join(tplant, data, by = c("t", "plant"))
    new[which(apply(new, 1, function(x)
      x$id %in% unlist(x$n) & is.na(x$circ))), "circ"] = max.circ + 1L
    data = new %>% select(-n)
    
    if (length(which(is.na(data$circ))) == 0)
      return(data)
    else
      return(add.circles(data, unmapped))
  }
  if (sum(!is.na(segm.data$circ)) > 1)
    segm.data = add.circles(segm.data, prev.unmapped)
  
  quant.data = quant.data %>% left_join(segm.data %>% select(-one_of("neighs")),
                                        by = intersect(colnames(quant.data), colnames(segm.data %>% select(-one_of("neighs")))),
                                        all.x = TRUE)
  division.data = division.data %>% left_join(segm.data %>% select(-one_of("neighs")),
                                              by = intersect(colnames(division.data), 
                                                             colnames(segm.data %>% select(-one_of("neighs")))),
                                              all.x = TRUE)
  
  quant.data = quant.data %>% 
    select(plant, t, id, m.x, m.y, m.z, n.x, n.y, n.z, m.vol, n.vol, expr, d2t, circ, layer, neighs) %>% 
    arrange(plant,t,id)
  segm.data = segm.data %>% 
    select(plant, t, id, m.x, m.y, m.z, n.x, n.y, n.z, m.vol, n.vol, expr, d2t, circ, layer, neighs) %>% 
    arrange(plant,t,id)
  division.data = division.data %>% 
    select(plant, t, id, lineage, age, m.x, m.y, m.z, n.x, n.y, n.z, m.vol, n.vol, expr, d2t, circ, layer, neighs, everything()) %>% 
    arrange(plant, t, id, lineage, age)
  
  EXPR.MIN.DIV.FRACT = 0.7
  if(length(quant.files) > 0){  
    to.filter = c() # should be few
    # Backwards
    for(ii in 1:nrow(quant.data)) {
      index = which(timepoints == unlist(quant.data[ii, "t"]))
      if (index == 1 )
        next
      
      subline  = unlist(sublines[match(quant.data[ii, "id"], unlist(sublines[, index])), ])
      parent   = unlist(subline[index - 1])
      parent.t = timepoints[index - 1]
      
      parent.expr = unlist(quant.data[which(quant.data$id == parent & quant.data$t == parent.t), "expr"])
      
      if(length(parent.expr) > 0 && unlist(quant.data[ii, "expr"]) / parent.expr < EXPR.MIN.DIV.FRACT)
        to.filter = rbind(to.filter, quant.data[ii, ])
    }
    to.filter = interaction(to.filter$t, to.filter$id)
    quant.data[which(interaction(quant.data$t, quant.data$id) %in% to.filter),          "expr"] = NA
    segm.data[which(interaction(segm.data$t, segm.data$id) %in% to.filter),             "expr"] = NA
    division.data[which(interaction(division.data$t, division.data$id) %in% to.filter), "expr"] = NA
    
    # Forward
    to.filter = division.data[which(
      division.data$expr /
        mean(c(division.data$child1.expr, division.data$child2.expr), na.rm = TRUE) < EXPR.MIN.DIV.FRACT
    ), ]
    to.filter = interaction(to.filter$t, to.filter$id)
    quant.data[which(interaction(quant.data$t, quant.data$id) %in% to.filter),          "expr"] = NA
    segm.data[which(interaction(segm.data$t, segm.data$id) %in% to.filter),             "expr"] = NA
    division.data[which(interaction(division.data$t, division.data$id) %in% to.filter), "expr"] = NA
  }
  
  cat("\tRetrieving cell line / lineage data\n")
  # Get cell line data
  cell.lines = get.cell.lineages(sublines, lineages)
  
  # lineage.sublines.data = list()
  lineage.sublines.data = parLapply(cl, 
                                    cell.lines,
                                    get.lineage.sublines.data,
                                    quant.data = quant.data,
                                    timepoints = timepoints
  )
  
  # Take away empty fellas 
  lineage.sublines.data = 
    lineage.sublines.data[parLapply(cl, lineage.sublines.data, function(x) 
      length(unlist(x))) > 0]
  cell.lines.data = parLapply(cl, lineage.sublines.data, collapse.lineage) 
  
  # Sort by number of members?
  if(length(cell.lines.data) > 0){
    lineage.sublines.data = lineage.sublines.data[parSapply(cl, cell.lines.data, function (x)
      x[, "d2t"] %>% nna %>% sum) %>% unlist %>% order(decreasing = TRUE)]
  }
  stopCluster(cl)
  
  list(
    quant.data            = quant.data,
    segm.data             = segm.data,
    division.data         = division.data,
    sublines              = sublines,
    cell.lines            = cell.lines,
    lineage.sublines.data = lineage.sublines.data,
    cell.lines.data       = cell.lines.data,
    timepoints            = timepoints
  ) %>%
    return()
}

# segm.data %>% group_by(plant, t) %>% 
# mutate(m.d2t = (c(m.x, m.y, m.z) - first(2:4, order_by = d2t)))

# segm.data = segm.data %>% group_by(plant, t) %>%
#   mutate(m.d2t = sqrt((m.x - m.x[d2t == min(d2t, na.rm = TRUE)]) ** 2 + 
#                       (m.y - m.y[d2t == min(d2t, na.rm = TRUE)]) ** 2 + 
#                       (m.z - m.z[d2t == min(d2t, na.rm = TRUE)]) ** 2))
# quant.data = quant.data %>% group_by(plant, t) %>%
#   mutate(m.d2t = sqrt((m.x - m.x[d2t == min(d2t, na.rm = TRUE)]) ** 2 + 
#                       (m.y - m.y[d2t == min(d2t, na.rm = TRUE)]) ** 2 + 
#                       (m.z - m.z[d2t == min(d2t, na.rm = TRUE)]) ** 2))
# hist(quant.data$m.d2t, breaks = 1000)
# hist(segm.data$m.d2t, breaks = 1000)


# segm.data %>% group_by(plant, t) %>% arrange(d2t) %>% mutate(
# m.d2t = sqrt(sum((c(m.x, m.y, m.z) - c(m.x[1], m.y[1], m.z[1]))**2))
# )


# if(length(to.filter) > 0) {
#   div.filter.stats = division.data[to.filter,]  %>% select(plant, t, id)
#   division.data[to.filter, ]$expr = NA
#   quant.data[match(
#     interaction(div.filter.stats$plant, div.filter.stats$t, div.filter.stats$id), 
#     interaction(quant.data$plant, quant.data$t, quant.data$id), nomatch = integer()), ]$expr = NA
#   segm.data[match(interaction(div.filter.stats$plant, div.filter.stats$t, div.filter.stats$id),
#                   interaction(segm.data$plant, segm.data$t, segm.data$id), nomatch = integer()), ]$expr = NA
# }

