map.parent = function(mapping.data, id.maps, m.time, i.time, mapping, shift = 10000){
  parent.id = mapping.data[[m.time]][[mapping]][1] 
  row.index = which(id.maps[[i.time]][, 1] == parent.id)[1]
  not.found = is.na(row.index)
  
  if (not.found)
    mapping.data[[m.time]][[mapping]][1] = parent.id + shift
  else
    mapping.data[[m.time]][[mapping]][1] = id.maps[[i.time]][row.index, 2]
  
  return(mapping.data[[m.time]][[mapping]][1])
}

map.children = function(mapping.data, id.maps, m.time, i.time, mapping, shift = 10000){
  result = sapply(mapping.data[[m.time]][[mapping]][-1], function(child.id) {
    row.index = which(id.maps[[i.time]][, 1] == child.id)[1]
    not.found = is.na(row.index)
    
    if (not.found) id = child.id + shift
    else           id = id.maps[[i.time]][row.index, 2]
    id
  })
  return(result)
}

###############################################################################
# READ DATA
###############################################################################
read.data = function(plant.no,
                     map.quant       = TRUE,
                     only.l1         = FALSE,
                     only.l2         = FALSE,
                     quant.missing   = -1,
                     mapping.missing = -1,
                     interval        = 4,
                     MQ_CUT          = 0.3) {
  # if(plant.no == 4) quant.missing = c(24)
  # if(plant.no == 15) quant.missing = c(12)
  # 
  print(quant.missing)
  if (quant.missing[1] < 0 | is.na(quant.missing[1])){
    quant.missing   = integer()
  } else {
    quant.missing   = quant.missing
  } 
  if (mapping.missing[1] < 0 | is.na(mapping.missing[1])) {
    mapping.missing = integer()
  } else {
    mapping.missing = mapping.missing
  }
  print(quant.missing)
  
  # Read in the right files
  quant.files   = list.files(paste0("../data/clv3_complete/plant", plant.no, "/Results"), full.names = TRUE, pattern=".txt")
  mapping.files = list.files(paste0("../data/newData/plant", plant.no, "/tracking_data"), full.names = TRUE, pattern=".dat")
  l1.files      = list.files(paste0("../data/newData/plant", plant.no, "/segmentation_data"), full.names = TRUE, pattern="L1.dat")
  neigh.files   = list.files(paste0("../data/newData/plant", plant.no, "/segmentation_data"), full.names = TRUE, pattern="neighbourhood.dat")
  vol.files     = list.files(paste0("../data/newData/plant", plant.no, "/segmentation_data"), full.names = TRUE, pattern="volumes.dat")
  
  timepoints    = seq(0, max(max(extract.numbers(mapping.files)[3,]), ifelse(plant.no == 1, -1,
                             extract.numbers(quant.files)[3,])), by = interval)

  no.track.files = length(mapping.files)
  no.quant.files = length(quant.files)
  
  # Put in right order
  if(no.quant.files > 0) {
    order.quant = extract.numbers(quant.files)[3, ]
    quant.files = quant.files[order(order.quant)]
  }
  
  order.plant = extract.numbers(mapping.files)[3, ]
  order.l1    = extract.numbers(l1.files)[2, ]
  order.neigh = extract.numbers(neigh.files)[2, ]
  
  mapping.files = mapping.files[order(order.plant)]
  l1.files      = l1.files[order(order.l1)]
  neigh.files   = neigh.files[order(order.neigh)]
  
  # TODO: this is the ugliest piece of code I have ever seen
  mapping.data = lapply(1:length(mapping.files), function(x) {
    data        = read.table(mapping.files[x], header = FALSE)
    data[, 1:2] = data[, 2:1] # idk why this is
    data        = data[data[, 3] > MQ_CUT, -3]
    
    data = split(data, seq(nrow(data)))
    data = lapply(data, function(x) unname(unlist(as.vector(x), recursive = F)))
    data = unname(data)

    parents = sapply(data, "[[", 1)    
    dups    = which(duplicated(parents))
    for(ii in data[dups]) {
      # print(ii)
      id = match(ii[1], parents)
      data[[id]] = c(data[[id]], ii[2])
    }
    data = data[-dups]
    data
  })
  
  # Read in quantified data, add the missing ones
  quant.data = lapply(quant.files, function(x) read.table(x, header = T, sep = "\t"))
  for(ii in (quant.missing/4)) 
    quant.data = c(quant.data[1:ii], list(
      matrix(NA, nrow=0, ncol = ncol(quant.data[[1]]), 
      dimnames = list(NULL, colnames(quant.data[[1]])))), 
      quant.data[(ii + 1):length(quant.data)])
  
  # Read in L1 cells
  l1.data = lapply(l1.files, scan)
  
  # Read in neighbourhood data
  neigh.data = lapply(neigh.files, function(x) {
    neighs = scan(x, what = "", sep = "\n")
    neighs = sapply(strsplit(neighs, "\\s+"), as.integer)
    neighs = neighs[sapply(neighs, function(x) x[1]) > 1] # remove background neighborhood
    # neighs = lapply(neighs, function(x) x[x > 1] - 2) # remove background, shift to get right correspondence
    neighs
  })
  
  # Read in L2 cells (these should be the neighbours to L1)
  l2.data = lapply(1:length(neigh.data), function(x) {
    l1.neighs = neigh.data[[x]][lapply(neigh.data[[x]], "[[", 1) %in% l1.data[[x]]]
    l1l2 = unique(unlist(l1.neighs))
    l2 = l1l2[!l1l2 %in% l1.data[[x]]]
    l2  
  })
    
  if (only.l1) {
    if (no.quant.files > 0) {
      quant.data = lapply(1:length(quant.data), function(x) 
        quant.data[[x]][quant.data[[x]][, "Cell.id"] %in% l1.data[[x]], ])
    }
    # TODO: get mapping data only for l1
    # if(no.mapping.files > 0) ...
  } else if (only.l2) {
    if (no.quant.files > 0) {
      quant.data = lapply(1:length(quant.data), function(x) 
        quant.data[[x]][quant.data[[x]][, "Cell.id"] %in% l2.data[[x]], ])
    } 
    # TODO: get mapping data only for l2
    # if(no.mapping.files > 0) ...
  }

  if (map.quant && 
      plant.no != 18 && 
      plant.no != 1) {
    id.files       = list.files(paste0("../data/plant", plant.no, "_correspondence"), full.names = TRUE, pattern = "correspondence")
    id.files.order = extract.numbers(id.files)[3,]
    id.files       = id.files[order(id.files.order)]
    id.maps        = lapply(id.files, function(x) read.table(x, header = FALSE, sep = ","))
    for(ii in (quant.missing/4)) {
      id.maps = c(id.maps[1:ii], list(
        matrix(NA, nrow=0, ncol = ncol(id.maps[[1]]), 
               dimnames = list(NULL, colnames(id.maps[[1]])))), 
        id.maps[(ii + 1):length(id.maps)])
    }
    
    # TODO: Some cells have multiple nuclei. This is due to errors in the membrane segmentation.
    ### Adjust cell id according to whether they are in quantified data or not.
    no.mapping.data = length(mapping.data)
    for (time in 1:(no.mapping.data)) {
      no.mapping.events = length(mapping.data[[time]])
      for (mapping in 1:no.mapping.events) {
        mapping.data[[time]][[mapping]][1]  = map.parent(mapping.data,   id.maps, time, time,     mapping)
        mapping.data[[time]][[mapping]][-1] = map.children(mapping.data, id.maps, time, time + 1, mapping)
      }
    }
  }
  
  # Add distance to the top
  if(plant.no != 1){
    top.coords = lapply(1:length(quant.data), function(time) get.top.coordinates(quant.data, time))
    dist2tops  = lapply(1:length(quant.data), function(time)
      apply(quant.data[[time]], 1, function(x)
        distance(x[2:4], top.coords[[time]])))
    quant.data = lapply(1:length(quant.data), function(x)
      cbind(quant.data[[x]], dist2top = dist2tops[[x]]))
  }
  
  for(ii in quant.data) colnames(ii) = c("Cell.id", "x", "y", "z", "Boa.volume", "Mean.cell.intensity", "dist2top")
  
  return(list(
    mapping.data = mapping.data,
    quant.data   = quant.data,
    timepoints   = timepoints
  ))
}

get.q.timepoints = function(plant.no) {
  quant.files = list.files(paste0("../data/clv3_complete/plant", plant.no, "/Results"), full.names = TRUE, pattern=".txt")
  times = extract.numbers(quant.files)[3,]
  sort(times)  
}

