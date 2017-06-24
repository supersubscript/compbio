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
read.data = function(plant.no, map.quant = TRUE) {
  
  # Read in the right files
  mapping.files = list.files(paste0("../data/PNAS/plant", plant.no, "/tracking_data"), full.names = TRUE)
  quant.files   = list.files(paste0("../data/clv3_complete/plant", plant.no, "/Results"), full.names = TRUE, pattern=".txt")
  
  # TODO: Should we couple this with the quantified files, as some lack points in tracking, some in quantification?
  timepoints    = extract.numbers(mapping.files)[3,]
  
  no.track.files = length(mapping.files)
  no.quant.files = length(quant.files)
  
  # Put in right order
  order.plant = extract.numbers(mapping.files)[3,]
  order.quant = extract.numbers(quant.files)[3,]
  
  mapping.files = mapping.files[order(order.plant)]
  quant.files   = quant.files[order(order.quant)]
  timepoints    = c(0, timepoints[order(order.plant)])
  
  # Read in the mapping data between timepoints
  mapping.data = lapply(1:length(mapping.files), function(x) {
    data = readLines(textConnection(gsub(":", ",", readLines(mapping.files[x]))))
    data = lapply(data, function(y) gsub(",", "", y))
    return(lapply(data, function(y) as.integer(str_split(y," ")[[1]])))
  })
  
  # Read in quantified data
  quant.data = lapply(quant.files, function(x) 
    read.table(x, header = T, sep = "\t"))
  if(length(quant.data) != 0){
    all.quant.data = lapply(1:ncol(quant.data[[1]]), function(x) 
      t(ldply(unname(sapply(quant.data, "[", x)), rbind)))
    names(all.quant.data) = colnames(quant.data[[1]])
  }
  
  if (map.quant) {
    id.files       = list.files(paste0("../data/plant", plant.no, "_correspondence"), full.names = T, pattern="correspondence")
    id.files.order = extract.numbers(id.files)[3,]
    id.files       = id.files[order(id.files.order)]
    id.maps        = lapply(id.files, function(x) read.table(x, header = F, sep = ","))
    
    # TODO: Some cells have multiple nuclei. This is due to errors in the membrane segmentation.
    ### Adjust cell id according to whether they are in quantified data or not.
    no.mapping.data = length(mapping.data)
    for (time in 1:(no.mapping.data)) {
      no.mapping.events = length(mapping.data[[time]])
      # print(time) 
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
  return(list(
    mapping.data = mapping.data,
    quant.data   = quant.data,
    timepoints   = timepoints
  ))
}

