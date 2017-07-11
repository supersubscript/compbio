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
read.data = function(plant.no, map.quant = TRUE, only.l1 = FALSE, only.l2 = FALSE) {
  
  # Read in the right files
  mapping.files = list.files(paste0("../data/PNAS/plant", plant.no, "/tracking_data"), full.names = TRUE)
  quant.files   = list.files(paste0("../data/clv3_complete/plant", plant.no, "/Results"), full.names = TRUE, pattern=".txt")
  
  # TODO: Should we couple this with the quantified files, as some lack points in tracking, some in quantification?
  timepoints    = extract.numbers(mapping.files)[3,]
  
  no.track.files = length(mapping.files)
  no.quant.files = length(quant.files)
  
  # Put in right order
  if(plant.no != 1) {
    order.quant = extract.numbers(quant.files)[3,]
    quant.files = quant.files[order(order.quant)]
  }
  
  order.plant = extract.numbers(mapping.files)[3,]
  mapping.files = mapping.files[order(order.plant)]
  timepoints    = c(0, timepoints[order(order.plant)])
  
  # Read in the mapping data between timepoints
  mapping.data = lapply(1:length(mapping.files), function(x) {
    data = readLines(textConnection(gsub(":", ",", readLines(mapping.files[x]))))
    data = lapply(data, function(y) gsub(",", "", y))
    return(lapply(data, function(y) as.integer(str_split(y," ")[[1]])))
  })
  
  # Read in quantified data
  quant.data = lapply(quant.files, function(x) read.table(x, header = T, sep = "\t"))
  

  
  if (map.quant && 
      plant.no != 18 && 
      plant.no != 1) {
    id.files       = list.files(paste0("../data/plant", plant.no, "_correspondence"), full.names = T, pattern="correspondence")
    id.files.order = extract.numbers(id.files)[3,]
    id.files       = id.files[order(id.files.order)]
    id.maps        = lapply(id.files, function(x) read.table(x, header = F, sep = ","))
    
    # TODO: Some cells have multiple nuclei. This is due to errors in the membrane segmentation.
    ### Adjust cell id according to whether they are in quantified data or not.
    no.mapping.data = length(mapping.data)
    for (time in 1:(no.mapping.data)) {
      no.mapping.events = length(mapping.data[[time]])
      
      ### Account for the missing data points in quantification for plant 4 and 15
      if ((plant.no == 4  && timepoints[time] == 20) ||
          (plant.no == 15 && timepoints[time] == 8)) {
        
        # Add only parents
        for (mapping in 1:length(mapping.data[[time]]))
          mapping.data[[time]][[mapping]][1] = map.parent(mapping.data, id.maps, time, time, mapping)
        
        # Add only children (in next)
        for (mapping in 1:length(mapping.data[[time + 1]]))
          mapping.data[[time + 1]][[mapping]][-1] = map.children(mapping.data, id.maps, time + 1, time + 1, mapping)
        
        # Add the rest, but nicely
        for (new.time in (time + 2):no.mapping.data) {
          for (mapping in 1:length(mapping.data[[new.time]])) {
            mapping.data[[new.time]][[mapping]][1]  = map.parent(mapping.data,   id.maps, new.time, new.time - 1, mapping)
            mapping.data[[new.time]][[mapping]][-1] = map.children(mapping.data, id.maps, new.time, new.time,     mapping)
          }
        }
        break # End this mess
      } else {
        ### Normal mapping, when data points not missing
        for (mapping in 1:no.mapping.events) {
          mapping.data[[time]][[mapping]][1]  = map.parent(mapping.data,   id.maps, time, time,     mapping)
          mapping.data[[time]][[mapping]][-1] = map.children(mapping.data, id.maps, time, time + 1, mapping)
        }
      }
    }
  }
  
  if(plant.no != 1 && only.l1){
    quant.data[[1]] = quant.data[[1]][quant.data[[1]][, 1] %in% get.mom(mapping.data[[1]]), ]
    for(ii in 2:length(mapping.data)) {
      
      # TODO: Make this look nice
      if (plant.no == 18 & timepoints[ii] == 44) {
        quant.data[[ii]] = matrix(NA, nrow = 0, ncol = ncol(quant.data[[ii]]))
        colnames(quant.data[[ii]]) = colnames(quant.data[[ii]] - 1)
      } else if (plant.no == 18 & timepoints[ii] == 48) {
        quant.data[[ii]] = quant.data[[ii]][quant.data[[ii]][,1] %in% 
                                              get.mom(mapping.data[[ii - 1]]), ]
      } else if (plant.no == 18 & timepoints[ii] >= 52) {
        quant.data[[ii]] = quant.data[[ii]][quant.data[[ii]][,1] %in% 
                                              union(get.mom(mapping.data[[ii - 1]]), 
                                                    get.dau(mapping.data[[ii - 2]])), ]
      } else if ((plant.no == 4  & timepoints[[ii]] == 24) ||
                 (plant.no == 15 & timepoints[[ii]]  == 12)) {
        quant.data[[ii]] = quant.data[[ii]][quant.data[[ii]][, 1] %in% 
                                              get.dau(mapping.data[[ii]]), ]
      } else if ((plant.no == 4  & timepoints[[ii]] == timepoints[length(timepoints) - 1]) ||
                 (plant.no == 15 & timepoints[[ii]] == timepoints[length(timepoints) - 1])) {
        quant.data[[ii]] = quant.data[[ii]][quant.data[[ii]][, 1] %in% 
                                              get.dau(mapping.data[[ii]]), ]
      } else if ((plant.no == 4  & timepoints[[ii]] > 24) ||
                 (plant.no == 15 & timepoints[[ii]] > 12)) {
        quant.data[[ii]] = quant.data[[ii]][quant.data[[ii]][, 1] %in% 
                                              union(get.mom(mapping.data[[ii + 1]]), 
                                                    get.dau(mapping.data[[ii]])), ]
      } else
        quant.data[[ii]] = quant.data[[ii]][quant.data[[ii]][,1] %in% 
                                              union(get.mom(mapping.data[[ii]]), 
                                                    get.dau(mapping.data[[ii - 1]])), ]
    }
    if (plant.no == 18) {
      quant.data[[length(timepoints) + 1]] =
        quant.data[[length(timepoints) + 1]][quant.data[[length(timepoints) + 1]][, 1] %in%
                                               get.dau(mapping.data[[length(timepoints) - 1]]), ]
    }
  }  
  # Actually it's !L1
  if(plant.no != 1 && only.l2) {
    quant.data[[1]] = quant.data[[1]][!quant.data[[1]][, 1] %in% get.mom(mapping.data[[1]]), ]
    for(ii in 2:length(mapping.data)) {
      # TODO: Make this look nice
      if (plant.no == 18 & timepoints[ii] == 44) {
        quant.data[[ii]] = matrix(NA, nrow = 0, ncol = ncol(quant.data[[ii]]))
        colnames(quant.data[[ii]]) = colnames(quant.data[[ii]] - 1)
      } else if (plant.no == 18 & timepoints[ii] == 48) {
        quant.data[[ii]] = quant.data[[ii]][!quant.data[[ii]][,1] %in% 
                                              get.mom(mapping.data[[ii - 1]]), ]
      } else if (plant.no == 18 & timepoints[ii] >= 52) {
        quant.data[[ii]] = quant.data[[ii]][!quant.data[[ii]][,1] %in% 
                                              union(get.mom(mapping.data[[ii - 1]]), 
                                                    get.dau(mapping.data[[ii - 2]])), ]
      } else if ((plant.no == 4  & timepoints[[ii]] == 24) ||
                 (plant.no == 15 & timepoints[[ii]]  == 12)) {
        quant.data[[ii]] = quant.data[[ii]][!quant.data[[ii]][, 1] %in% 
                                              get.dau(mapping.data[[ii]]), ]
      } else if ((plant.no == 4  & timepoints[[ii]] == timepoints[length(timepoints) - 1]) ||
                 (plant.no == 15 & timepoints[[ii]] == timepoints[length(timepoints) - 1])) {
        quant.data[[ii]] = quant.data[[ii]][!quant.data[[ii]][, 1] %in% 
                                              get.dau(mapping.data[[ii]]), ]
      } else if ((plant.no == 4  & timepoints[[ii]] > 24) ||
                 (plant.no == 15 & timepoints[[ii]] > 12)) {
        quant.data[[ii]] = quant.data[[ii]][!quant.data[[ii]][, 1] %in% 
                                              union(get.mom(mapping.data[[ii + 1]]), 
                                                    get.dau(mapping.data[[ii]])), ]
      } else
        quant.data[[ii]] = quant.data[[ii]][!quant.data[[ii]][,1] %in% 
                                              union(get.mom(mapping.data[[ii]]), 
                                                    get.dau(mapping.data[[ii - 1]])), ]
    }
    if (plant.no == 18) {
      quant.data[[length(timepoints) + 1]] =
        quant.data[[length(timepoints) + 1]][!quant.data[[length(timepoints) + 1]][, 1] %in%
                                               get.dau(mapping.data[[length(timepoints) - 1]]), ]
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

