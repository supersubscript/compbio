###############################################################################
###############################################################################
get.mom = function(x) {unlist(lapply(x, "[[", 1))}
get.dau = function(x) {unlist(lapply(x, "[", -1))}

###############################################################################
# Read in all the cell lines that we observe throughout the data set
###############################################################################
get.cell.lines = function(all.mapping.data, remove.empty = TRUE, sort.order = "n"){
  cell.lines = cbind(cbind(t(t(cbind(lapply(all.mapping.data[[1]], "[[", 1), 
                                     lapply(all.mapping.data[[1]], "[", -1))))), 
                     matrix(list(NA), ncol = length(all.mapping.data) - 2 + 1, 
                            nrow = length(all.mapping.data[[1]])))
  for(time in 2:length(all.mapping.data)){
    this.mapping.data = all.mapping.data[[time]]
    parents           = lapply(this.mapping.data, "[[", 1)
    children          = lapply(this.mapping.data, "[", -1)
    parents.found     = which(sapply(parents, function(x) x %in% unlist(cell.lines[, time])))
    parents.not.found = !1:length(parents) %in% parents.found
    
    # Add the children of the parents' which are found
    for(ii in parents.found){
      for(line in 1:nrow(cell.lines)){
        if(parents[[ii]] %in% cell.lines[[line, time]]){
          already.there = cell.lines[[line, time + 1]]
          cell.lines[[line, time + 1]] = c(already.there[!is.na(already.there)], as.integer(children[[ii]]))
          break
        }
      }
    }
    
    ### Create new cell lines for the ones which aren't found
    # Fill out times until this one with NAs
    new.lines = cbind(cbind(matrix(
      list(NA),
      ncol = time - 1,
      nrow = length(all.mapping.data[[time]][-parents.found])
    ),
    # Glue on new line with mother and daughters in separate timepoints
    t(t(
      cbind(
        lapply(all.mapping.data[[time]][-parents.found], "[[", 1), # mother
        lapply(all.mapping.data[[time]][-parents.found], "[", -1)  # daughters
      )
    ))),
    # Fill up with NAs for the rest
    matrix(
      list(NA),
      ncol = length(all.mapping.data) - time,
      nrow = length(all.mapping.data[[time]][-parents.found])
    ))
    cell.lines = rbind(cell.lines, new.lines)
  }
  
  # Should we take away empty cell lines?
  no.members = apply(cell.lines, 1, function(x) sum(!is.na(unlist(x))))
  if (remove.empty) {
    cell.lines = cell.lines[no.members > 0, ]
  }
  # Should we sort by number of members?
  if (sort.order == "no.members") {
    cell.lines = cell.lines[order(-no.members),]
  }
  return(cell.lines)
}
###############################################################################
# Take out the cell with the highest CLV3 expression (this should be the 
# topmost one).
###############################################################################
get.top.coordinates = function(quant.data, time) {
  top.coords  = c(NA, NA, NA)
  high.to.low = order(-quant.data[[time]][, "Mean.cell.intensity"])
  top.cell    = quant.data[[time]][high.to.low,][1, ]
  top.coords  = top.cell[2:4]
  return(top.coords)
}

###############################################################################
###############################################################################
# get.cell.line.data = function(cell.lines,
#                               quant.data,
#                               timepoints,
#                               min.no.members = 0,
#                               param,
#                               no.z           = 1,
#                               no.xy          = no.z,
#                               sorting.param  = "n",
#                               sorting.fct    = mean) {
#   cell.lines.data = vector("list", nrow(cell.lines))
#   no.cell.lines   = nrow(cell.lines)
#   no.timepoints   = length(timepoints)
#   
#   ### Extract the quantified information for the cells in the cell lines.
#   for (ii in 1:no.cell.lines) {
#     lineage.data = matrix(NA, nrow = 0, ncol = 1 + ncol(quant.data[[1]])) # add time and dist2top
#     
#     # How long we will iterate for depends on whether we have missing data points or not
#     time.len = ifelse(plant.no == 4 || plant.no == 15,
#                       length(timepoints) - 1,
#                       length(timepoints))
#     
#     for (time in 1:time.len) {
#       mod.time = ifelse((plant.no == 4  && timepoints[time] >= 24) ||
#                           (plant.no == 15 && timepoints[time] >= 12), 
#                         time + 1, 
#                         time)
#       cells.in.lineage = which(quant.data[[time]]$Cell.id %in% cell.lines[[ii, mod.time]])
#       # top.coords       = get.top.coordinates(quant.data, time, param = param, no.z = no.z, no.xy = no.xy)
#       
#       # If the cell line finds any cells which have quantified information 
#       # about them, calculate distance to the center and add as a parameter.
#       no.cells.in.lineage = length(cells.in.lineage)
#       if (no.cells.in.lineage > 0) {
#         for (kk in 1:no.cells.in.lineage) {
#           # dist2tops    = distance(quant.data[[time]][cells.in.lineage[kk], 2:4], top.coords)
#           lineage.data = rbind(lineage.data,
#                                c(Time = timepoints[mod.time],
#                                  quant.data[[time]][cells.in.lineage[kk], ]))
#         }
#       } else 
#         lineage.data = rbind(lineage.data, rep(NA, ncol(lineage.data)))
#     }
#     cell.lines.data[[ii]] = lineage.data
#   }
#   
#   # Unlist
#   cell.lines.data = lapply(cell.lines.data, function(x) apply(x, 1:2, function(y) unlist(y)))
#   
#   # Remove cell lines with too few members
#   cell.lines.data = cell.lines.data[!sapply(cell.lines.data, function(x)
#     sum(!is.na(unlist(x))) < min.no.members * ncol(quant.data[[1]]))]
#   
#   # Sort the cell lines according to the set rule
#   if (sorting.param != "n") {
#     cell.lines.order = sapply(cell.lines.data, function(x) sorting.fct(x[, sorting.param], na.rm = TRUE))
#     cell.lines.data  = cell.lines.data[order(cell.lines.order)]
#   }
#   
#   return(cell.lines.data)
# }
###############################################################################
# Process which cells have which mother at what time
###############################################################################
# get.mothers = function(all.mapping.data, timepoints) {
#   mothers = lapply(1:length(all.mapping.data), function(time)
#     lapply(all.mapping.data[[time]], function(this.mapping.data) {
#       
#       mother.id = this.mapping.data[1]
#       cell.ids  = this.mapping.data[-1]
#       
#       mums      = rep(mother.id, length(cell.ids))
#       times     = rep(timepoints[time + 1], length(cell.ids))
#       
#       return(cbind(
#         t         = times,
#         mother.id = mums,
#         cell.id   = cell.ids
#       ))
#     }))
#   
#   # Bind the mothers together in a matrix
#   mothers = do.call(rbind, unlist(mothers, recursive = F))
#   
#   # Add the first generation
#   first.gen.t          = rep(0, length(all.mapping.data[[1]]))
#   first.gen.mother.ids = rep(NA, length(all.mapping.data[[1]]))
#   first.gen.cell.ids   = sapply(all.mapping.data[[1]], "[[", 1)
#   
#   mothers = rbind(cbind(t         = first.gen.t, 
#                         mother.id = first.gen.mother.ids, 
#                         cell.id   = first.gen.cell.ids), 
#                   mothers)
#   
#   return(mothers)
# }

###############################################################################
###############################################################################
get.division.events = function(mapping.data, quant.data, sublines, plant.no){
  no.quant.data.files = length(quant.data)
  
  # Preallocate space
  no.division.events = sum(unlist(lapply(mapping.data, function(x) length(x[sapply(x, length) > 2]))))
  
  output.data = matrix(NA, ncol = 17, nrow = no.division.events)
  colnames(output.data) = c("cell.id", "lineage.id", "mother.id", "daughter1.id", 
                            "daughter2.id", "age", "t","x","y","z", "vol", "expr", 
                            "dist2top", "daughter1.vol", "daughter2.vol", "daughter1.expr", "daughter2.expr")
  count = 1
  for(time in 1:length(mapping.data)){
    division.data = mapping.data[[time]]
    division.data = division.data[sapply(division.data, length) > 2] # Only keep division events
    
    if((plant.no == 4  && timepoints[time] == 20) ||
       (plant.no == 15 && timepoints[time] == 8)) {
      p.time = time
      d.time = NA
    } else if((plant.no == 4  && timepoints[time] > 20) ||
              (plant.no == 15 && timepoints[time] > 12)) {
      p.time = ifelse((plant.no == 4  && timepoints[time] == 24) ||
                      (plant.no == 15 && timepoints[time] == 12), 
                      NA, time - 1)
      d.time = time
    } else if (plant.no == 18 && timepoints[time] == 36){
      p.time = time
      d.time = time + 2
    } else if (plant.no == 18 && timepoints[time] >= 44){
      p.time = time + 1
      d.time = time + 2
    } else {
      p.time = time
      d.time = time + 1
    }
    
    # for every cell, take out the data FROM quant.data that is 
    # interesting for that and add it to a row. 
    for(ii in 1:length(division.data)){
      xyzvid = rep(NA, 6); this.cell.index = NA; daughter1.index = NA; daughter2.index = NA
      daughter1.vol = NA;  daughter2.vol   = NA; daughter1.expr  = NA; daughter2.expr  = NA
      d1.d2t = NA; d2.d2t = NA
      
      if (plant.no != 1) {
        this.cell.index = match(division.data[[ii]][1], quant.data[[p.time]][, "Cell.id"])
        daughter1.index = match(division.data[[ii]][2], quant.data[[d.time]][, "Cell.id"])
        daughter2.index = match(division.data[[ii]][3], quant.data[[d.time]][, "Cell.id"])
        
        # Does this cell exist in the quantified data?
        # print(this.cell.index)
        # print(length(this.cell.index))
        xyzvid         = c(unlist(quant.data[[p.time]][this.cell.index, -1]))
        if(is.null(xyzvid))
          xyzvid = rep(NA, 6)
        daughter1.expr = quant.data[[d.time]][daughter1.index, "Mean.cell.intensity"]
        daughter1.vol  = quant.data[[d.time]][daughter1.index, "Boa.volume"]
        daughter2.expr = quant.data[[d.time]][daughter2.index, "Mean.cell.intensity"]
        daughter2.vol  = quant.data[[d.time]][daughter2.index, "Boa.volume"]
        daughter1.vol  = ifelse(is.null(daughter1.vol),  NA, daughter1.vol)
        daughter1.expr = ifelse(is.null(daughter1.expr), NA, daughter1.expr)
        daughter2.vol  = ifelse(is.null(daughter2.vol),  NA, daughter2.vol)
        daughter2.expr = ifelse(is.null(daughter2.expr), NA, daughter2.expr)
      }
      
      ### Get the identifiers
      this.cell.id = division.data[[ii]][1]
      lineage.id   = extract.numbers(rownames(sublines)[match(this.cell.id, sublines[, time])])[1, ]
      mother.id    = sublines[match(this.cell.id, sublines[, time]), time - 1][1]
      
      ### Find how long the mother has been alive for
      age = 0
      subline = sublines[match(this.cell.id, sublines[, time]), ]
      for (kk in time:2) {
        parent      = subline[kk - 1]
        lines       = sublines[which(sublines[, kk - 1] == parent), ]
        no.siblings = length(unique(lines[nna(lines[, kk]), kk]))
        
        # print(no.siblings)
        if(no.siblings > 1) break
        if(no.siblings == 0) {age = NA; break}
        
        # If this is the only offspring
        age = age  + timepoints[kk + 1] - timepoints[kk]
      }
      
      # Print output data
      output.data[count, ] = c(
        this.cell.id   = this.cell.id,
        lineage.id     = lineage.id,
        mother.id      = mother.id,
        daughter1.id   = division.data[[ii]][2],
        daughter2.id   = division.data[[ii]][3],
        age            = age,
        t              = timepoints[time],
        xyzvid,
        daughter1.vol  = daughter1.vol,
        daughter2.vol  = daughter2.vol,
        daughter1.expr = daughter1.expr / xyzvid[5],
        daughter2.expr = daughter2.expr / xyzvid[5]
      )
      count = count + 1
    }
  }
  return(output.data)
}

###############################################################################
### Get sublines from normal lineages
###############################################################################
get.sublines = function(mapping.data, timepoints) {
  
  # Add first generation
  lines = do.call(rbind, unlist(sapply(mapping.data[[1]], function(mapping)
    lapply(mapping[-1], function(child) c(mapping[1], child))), recursive = F))
  lines = cbind(lines, matrix(NA, ncol = length(timepoints) - 2, nrow = nrow(lines)))
  
  ### Add the rest
  for (time in 2:(length(timepoints) - 1)) {
    new.lines = matrix(NA, nrow = 0, ncol = length(timepoints)) # this will hold update
    fill1     = rep(NA, time - 1)
    fill2     = rep(NA, length(timepoints) - time - 1)
    found     = c()
    
    for (mapping in 1:length(mapping.data[[time]])) {
      parent   = mapping.data[[time]][[mapping]][1] 
      children = mapping.data[[time]][[mapping]][-1]
      match.id = match(parent, lines[, time])
      
      # If the parent is found in an earlier cell line, add to this
      # otherwise, create and fill out new line
      if (!is.na(match.id)) {
        new.lines = rbind(new.lines, t(sapply(children, function(child)
          replace(lines[match.id, ], time + 1, child))))
        found = c(found, match.id) # to remove old cell line
      } else {
        new.lines = rbind(new.lines, t(sapply(children, function(child) 
          c(fill1, parent, child, fill2))))
      }
    }
    if(length(found) != 0)
      new.lines = rbind(lines[-found, ], new.lines)
    lines = new.lines
  }
  lines.sorted = lines[do.call(order, lapply(1:NCOL(lines), function(ii) lines[, ii])), ] # order by all columns
  counter1 = 0
  counter2 = 1
  prev.first = -1
  prev.time  = -1
  rownames(lines.sorted) = 1:nrow(lines.sorted)
  
  ### Set rownames according to lineage
  for (line in 1:nrow(lines.sorted)) {
    this.first = lines.sorted[line, which(!is.na(lines.sorted[line, ]))][1]
    this.time  = which(lines.sorted[line, ] == this.first)[1]
    
    # Add to lineage
    if (this.first == prev.first && this.time == prev.time) {
      counter2 = counter2 + 1
      rownames(lines.sorted)[line] = paste0(counter1, ".", counter2)
    } else {
      counter1 = counter1 + 1
      counter2 = 1
      rownames(lines.sorted)[line] = paste0(counter1, ".", counter2)
      prev.first = this.first
      prev.time  = this.time
    }
  } 
  lines.sorted
}

###############################################################################
###############################################################################
get.lineage.data = function(lineage, quant.data, timepoints) {
  # lineage = lineages[["803"]]
  no.timepoints = length(timepoints)
  
  time.len = ifelse(plant.no == 4 || plant.no == 15,
                    no.timepoints - 1,
                    no.timepoints)
  
  lineage.data = matrix(NA,
                        nrow = sum(nna(lineage) & unlist(lineage) < 10000),
                        ncol = 1 + ncol(quant.data[[1]]) + 1) # add time and dist2top
  
  counter = 1
  for (time in 1:time.len) {
    m.time = ifelse((plant.no == 4  && timepoints[time] >= 24) ||
                    (plant.no == 15 && timepoints[time] >= 12),
                    time + 1,
                    time)
    q.time = ifelse(plant.no == 18 && timepoints[time] >= 40, time + 1, time)
    
    top.coords = get.top.coordinates(quant.data, q.time, param = param, no.z = no.z, no.xy = no.xy)
    
    for (member in 1:nrow(lineage)) {
      qd.match = match(lineage[member, m.time], quant.data[[q.time]]$Cell.id)
      
      # If the member exists in the quantified data, add it
      if (nna(qd.match)) {
        # dist2tops = distance(quant.data[[time]][qd.match, 2:4], top.coords)
        lineage.data[counter,] =
          c(timepoints[m.time],
            unlist(quant.data[[q.time]][qd.match, ]))
        counter = counter + 1
      }
    }
  }
  colnames(lineage.data) = c("Time", colnames(quant.data[[1]]))
  lineage.data
}

###############################################################################
###############################################################################
get.lineage.sublines.data = function(lineage, quant.data, timepoints) {
  no.timepoints = length(timepoints)
  no.q.cols     = ncol(quant.data[[1]]) + 1 # add time and dist2top
  
  time.len = ifelse(plant.no == 4 || plant.no == 15,
                    no.timepoints - 1,
                    no.timepoints)
  
  members       = apply(lineage, 1, function(x) sum(nna(x) & unlist(x) < 10000))
  sublines.data = lapply(members, function(x)
    matrix(NA, nrow = x, ncol = no.q.cols)) # preallocate space
  
  for (subline in 1:length(members)) {
    counter = 1
    for (time in 1:time.len) {
      m.time = ifelse((plant.no == 4  && timepoints[time] >= 24) ||
                        (plant.no == 15 && timepoints[time] >= 12),
                      time + 1,
                      time)
      q.time = ifelse(plant.no == 18  && timepoints[time] >= 40,
                      time + 1,
                      time)
      
      # top.coords = get.top.coordinates(quant.data, time) 
      qd.match = match(lineage[subline, m.time], quant.data[[q.time]]$Cell.id)
      
      # If the member exists in the quantified data, add it
      if (nna(qd.match)) {
        # dist2tops = distance(quant.data[[time]][qd.match, 2:4], top.coords)
        sublines.data[[subline]][counter,] =
          c(timepoints[m.time],
            unlist(quant.data[[q.time]][qd.match, ]))
        counter = counter + 1
      }
    }
    colnames(sublines.data[[subline]]) = c("Time", colnames(quant.data[[1]]))
  }
  sublines.data
}

###############################################################################
###############################################################################
# get.sublines.data = function(sublines,
#                              mapping.data,
#                              quant.data,
#                              timepoints,
#                              min.no.members,
#                              param,
#                              no.z           = 1,
#                              no.xy          = no.z,
#                              sorting.param  = "n",
#                              sorting.fct    = mean) {
#   sublines.data = vector("list", nrow(sublines))
#   no.sublines   = nrow(sublines)
#   no.timepoints = length(timepoints)
#   
#   ### Extract the quantified information for the cells in the cell lines.
#   for (ii in 1:no.sublines) {
#     lineage.data = matrix(NA, nrow = 0, ncol = 1 + ncol(quant.data[[1]]) + 1) # add time and dist2top
#     
#     # How long we will iterate for depends on whether we have missing data points or not
#     time.len = ifelse(plant.no == 4 || plant.no == 15,
#                       length(timepoints) - 1,
#                       length(timepoints))
#     
#     for (time in 1:time.len) {
#       mod.time = ifelse((plant.no == 4  && timepoints[time] >= 24) ||
#                           (plant.no == 15 && timepoints[time] >= 12),
#                         time + 1,
#                         time)
#       
#       qd.idx     = which(quant.data[[time]]$Cell.id %in% sublines[[ii, mod.time]])[1]
#       # top.coords = get.top.coordinates(quant.data, time, param = param, no.z = no.z, no.xy = no.xy)
#       
#       if(!is.na(qd.idx)){
#         # dist2top = distance(quant.data[[time]][qd.idx, 2:4], top.coords)
#         lineage.data = rbind(lineage.data,
#                              c(Time = timepoints[mod.time],
#                                quant.data[[time]][qd.idx, ]))
#       } else
#         lineage.data = rbind(lineage.data, rep(NA, ncol(lineage.data)))
#     }
#     sublines.data[[ii]] = lineage.data
#   }
#   return(sublines.data)
# }
###############################################################################
###############################################################################
get.lineages = function(sublines){
  subline.ids = rownames(sublines)
  lineage.ids = factor(sapply(subline.ids, function(x)
    strsplit(x, "\\.")[[1]][1]))
  return(split(as.data.frame(sublines), lineage.ids))
}
###############################################################################
###############################################################################
get.subline = function(sublines, line.id, all.lineages = NULL) {
  if (!is.null(all.lineages))
    return(all.lineages[line.id])
  else
    return(get.lineages(sublines)[line.id])
}
###############################################################################
###############################################################################
collapse.lineage = function(lineage) {
  all.rows = do.call(rbind, lineage)
  return(unique(all.rows))
}

###############################################################################
###############################################################################
get.roc = function(subline){
  rocs = lapply(1:length(subline), function(x) {
    intensities = subline[[x]][, "Mean.cell.intensity"]
    shift = 4
    if(length(intensities) > shift)
      return(cbind(subline[[x]][-((length(intensities)-(shift-1)):length(intensities)), "dist2top"], diff(intensities, lag = shift) / intensities[-((length(intensities)-(shift-1)):length(intensities))]))
  })
  
  result = unique(do.call(rbind, rocs))
  result
}

