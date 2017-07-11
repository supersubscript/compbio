get.mapping.data = function(file, CUTOFF = 0.3, remove.score = TRUE, remove.bg = TRUE) {
  mapping.data = system2("python", args = paste0(base.path, code.path, "/mapping_pkl2dat.py ", file), stdout = TRUE)
  mapping.data = read.table(text = mapping.data)
  mapping.data[, 1:2] = mapping.data[, 2:1]
  colnames(mapping.data) = c("parent", "child", "score")
  mapping.data = mapping.data[mapping.data[, "score"] > CUTOFF, ]
  
  if(remove.bg) 
    mapping.data = mapping.data[which(mapping.data[, "parent"] > 1), ]
  if (remove.score)
    mapping.data = mapping.data[, -3]
  
  as_tibble(mapping.data)
}

get.volume.data = function(file) {
  volume.data = system2("python", args = paste0(base.path, code.path, "/get_volumes.py ", file, " volumes"), stdout = TRUE)
  volume.data = read.table(text = volume.data, col.names = c("id", "volume"))
  as_tibble(volume.data[-1, ])
}

get.centers.data = function(file) {
  centers.data = system2("python", args = paste0(base.path, code.path, "/get_volumes.py ", file, " barycenter"), stdout = TRUE)
  centers.data = read.table(text = centers.data, col.names = c("id", "x", "y", "z"))
  centers.data = centers.data[-1, ]
  as_tibble(centers.data[, -1])
}

get.neigh.data = function(file, remove.bg = TRUE) {
  neigh.data = system2("python", args = paste0(base.path, code.path, "/get_volumes.py ", file, " neigh"), stdout = TRUE)
  neigh.data = sapply(neigh.data, function(x) as.integer(strsplit(x, "\\s+")[[1]]))
  if (remove.bg) {
    neigh.data = sapply(neigh.data, function(x) x[x > 1])
    neigh.data = neigh.data[-1]
  }
  neigh.data %>% unname
}

get.l1.data = function(file) {
  l1.data = system2("python", args = paste0(base.path, code.path, "/get_volumes.py ", file, " L1"), stdout = TRUE)
  l1.data = read.table(text = l1.data, col.names = "L1")
  l1.data %>% unlist
}

get.l2.data = function(l1.data, neigh.data) {
  l1.data   = unlist(l1.data)
  cells     = sapply(neigh.data, "[[", 1) 
  matches   = match(l1.data, cells, nomatch = integer())
  l1.neighs = unlist(neigh.data[matches])
  l2.data   = l1.neighs[!l1.neighs %in% l1.data]
  l2.data 
}

get.quant.data = function(file, d2t = TRUE) {
  quant.data = read.table(
    file,
    col.names = c("id", "n.x", "n.y", "n.z", "n.vol", "expr"),
    skip = 1
  )
  
  # Add distance to the top cell?
  if (d2t) {
    topcell    = which.max(quant.data[, "expr"])
    d2t        = sapply(1:nrow(quant.data), function(x)
      (quant.data[x, 2:4] - quant.data[topcell, 2:4]) %>%
        .**2 %>%
        sum  %>%
        sqrt
    )
    
    # sqrt(sum(** 2)))
    quant.data = cbind(quant.data, d2t = d2t)
  }
  quant.data %>% as_tibble()
}

empty.matr = function(nrow = 0, ncol = 7, colnames = c("id", "n.x", "n.y", "n.z", "n.vol", "expr", "d2t")) {
  matr = matrix(
    data = NA,
    nrow = nrow,
    ncol = ncol)
  colnames(matr) = colnames
  matr
}

fix.quant = function(quant.data, missing = NULL, interval = 4, timepoints = NULL) {
  if (length(quant.data)  == 0)
    quant.data = replicate(length(missing), empty.matr())
  else
    for (ii in (missing / interval)) {
      quant.data = c(quant.data[1:ii], list(empty.matr()),
                     quant.data[(ii + 1):length(quant.data)])
      quant.data[[ii + 1]] %>% as_tibble()
    }
  names(quant.data) = timepoints
  quant.data
}

fix.corresp = function(corr.data, missing = NULL, interval = 4, timepoints = NULL) {
  if (length(corr.data) == 0)
    corr.data = replicate(length(timepoints), empty.matr(ncol = 2, colnames = c("segm", "quant")))
  else
    for (ii in sort(missing / interval)) {
      corr.data = c(corr.data[1:ii],
                             list(empty.matr(
                               ncol = 2, colnames = c("segm", "quant")
                             )),
                             corr.data[(ii + 1):length(corr.data)])
      corr.data[[ii + 1]] %>% as_tibble()
    }
  names(corr.data) = timepoints
  corr.data
}

get.sublines = function(mapping.data, timepoints = NULL) {
  # This is a work of art
  map.previous = function(next.time, prev.time) {
    sublines = suppressWarnings(merge(prev.time, next.time, by.x = "child", by.y = "parent", all = TRUE, incomparables = NA))
    sublines[, 1:2] = sublines[, 2:1]
    sublines = unname(sublines)
    colnames(sublines)[1] = "parent"
    sublines 
  }
  
  # Retrieve cell lines
  sublines = rev(mapping.data) %>% reduce(map.previous)
  colnames(sublines) = timepoints
  
  # Order by members
  sublines = sublines[do.call(order, as.data.frame(sublines)), ] # sort
  sublines = sublines %>% as_tibble()
  sublines
}

get.corresp = function(file) {
  read.table(file, col.names = c("segm", "quant")) %>% as_tibble() 
}

fix.mapping = function(plant.mapping.data, missing = NULL, interval = 4, timepoints = NULL) {
  if (length(mapping.data) == 0)
    mapping.data = replicate(length(missing), empty.matr(ncol = 2, colnames = c("parent", "child")))
  else
    for (ii in sort((missing / interval))) {
      mapping.data = c(mapping.data[1:ii],
                       mapping.data[ii],
                       mapping.data[(ii + 1):length(mapping.data)])
      mapping.data[[ii]][, "child"] = -1 * mapping.data[[ii]][, "parent"]
      mapping.data[[ii]] = unique(mapping.data[[ii]]) # don't ruin the cell lines
      mapping.data[[ii + 1]][, "parent"] = -1 * mapping.data[[ii + 1]][, "parent"]
    }
  names(mapping.data) = paste(timepoints[-length(timepoints)], timepoints[-1], sep = "_to_")
  mapping.data
}

get.lineages = function(sublines) {
  rnames = 1:nrow(sublines)
  prev.first = -1; prev.time = -1; counter1 = 0; counter2 = 1
  first.times = apply(sublines, 1, function(x) nna(x) %>% which %>% min)
  first.vals  = sapply(1:nrow(sublines), function(x) 
    sublines[x, first.times[x]]) %>% 
    unname %>% 
    unlist
  
  for (line in 1:nrow(sublines)) {
    if (first.vals[line] == prev.first && first.times[line] == prev.time) {
      counter2 = counter2 + 1
      rnames[line] = counter1 #paste0(counter1, ".", counter2)
    } else {
      counter1 = counter1 + 1 #; counter2 = 1
      rnames[line] = counter1 #paste0(counter1, ".", counter2)
      prev.first = first.vals[line]; prev.time  = first.times[line]
    }
  }
  rnames
}

fix.segm = function(segm.data, missing = NULL, interval = 4, timepoints = NULL) {
  if (length(mapping.data) == 0)
    segm.data = replicate(length(missing), empty.matr(ncol = 6, colnames = c("id", "m.vol", "m.x", "m.y", "m.z", "neighs")))
  else
    for(ii in sort((missing / interval))) {
      segm.data = cbind(
        segm.data[, 1:ii],
        list(integer(), numeric(), numeric(), numeric(), numeric(), list()),
        # list(NULL, NULL, NULL, NULL, NULL, NULL), 
        segm.data[, ((ii + 1):ncol(segm.data))]
      )
    }
  colnames(segm.data) = timepoints
  rownames(segm.data) = c("id", "m.vol", "m.x", "m.y", "m.z", "neighs")
  segm.data
}

quant2segm.mapping = function(quant.data, corr.data, x) {
  this.qd = quant.data[[x]]
  this.cd = corr.data[[x]]
  matches = match(this.qd[, "id"] %>% unlist(), this.cd[, "quant"] %>% unlist())
  # this.qd = this.qd[nna(matches), ]
  new.ids = this.cd[matches, "segm"]
  this.qd[,"id"] = new.ids
  this.qd %>% as_tibble()
}

add.layers = function(segm.data, segm.files, neigh.data, cl = NULL) {
  get.layer.single = function(segm.data, l1.data, l2.data, cl = NULL) {
    # TODO: Add parallelisation
    ids = segm.data$id
    sapply(ids, function(x) 
      if(x %in% l1.data) return("L1") 
      else if (x %in% l2.data) return("L2") 
      else return("L3+"))
  }
  
  if(!is.null(cl))
    l1.data = parLapply(cl, segm.files, get.l1.data) 
  else
    l1.data = lapply(segm.files, get.l1.data)
  
  # TODO: Fix l1.data
  # TODO: Fix neigh.data
  
  l2.data = lapply(1:length(l1.data), function(x) get.l2.data(l1.data[[x]], neigh.data[[x]]))
  layer   = lapply(1:ncol(segm.data), function(x) factor(get.layer.single(segm.data[, x], l1.data[[x]], l2.data[[x]]), levels=c("L1","L2","L3+")))
  rbind(segm.data, layer)
}

get.division.event = function(ii, jj, sublines, mapping.data, segm.data, quant.data, timepoints, lineages, cl = NULL) {
  divs = mapping.data[[ii]] %>% 
    filter(parent == unlist(.[jj, "parent"])) %>% 
    .[1:2, ]
  
  parent.qd = quant.data[[ii]]     %>% 
    filter(.$id == divs$parent[1]) %>% 
    summarise(id = n(), n.x = mean(n.x), n.y = mean(n.y), n.z = mean(n.z), 
              n.vol = sum(n.vol), expr = mean(expr), d2t = mean(d2t))
  child1.qd = quant.data[[ii + 1]] %>% 
    filter(.$id == divs$child[1])  %>% 
    summarise(id = n(), n.x = mean(n.x), n.y = mean(n.y), n.z = mean(n.z), 
              n.vol = sum(n.vol), expr = mean(expr), d2t = mean(d2t))
  child2.qd = quant.data[[ii + 1]] %>% 
    filter(.$id == divs$child[2])  %>% 
    summarise(id = n(), n.x = mean(n.x), n.y = mean(n.y), n.z = mean(n.z), 
              n.vol = sum(n.vol), expr = mean(expr), d2t = mean(d2t))
  
  parent.sd = segm.data[, ii]      %>% as.tibble() %>% filter(.$id == divs$parent[1])
  child1.sd = segm.data[, ii + 1]  %>% as.tibble() %>% filter(.$id == divs$child[1])
  child2.sd = segm.data[, ii + 1]  %>% as.tibble() %>% filter(.$id == divs$child[2])
  
  colnames(child1.sd) = paste0("child1.", colnames(child1.sd))
  colnames(child1.qd) = paste0("child1.", colnames(child1.qd))
  colnames(child2.sd) = paste0("child2.", colnames(child2.sd))
  colnames(child2.qd) = paste0("child2.", colnames(child2.qd))
  
  if(nrow(parent.qd) == 0) parent.qd = parent.qd %>% add_row()
  if(nrow(child1.qd) == 0) child1.qd = child1.qd %>% add_row()
  if(nrow(child2.qd) == 0) child2.qd = child2.qd %>% add_row()
  if(nrow(parent.sd) == 0) parent.sd = parent.sd %>% add_row()
  if(nrow(child1.sd) == 0) child1.sd = child1.sd %>% add_row()
  if(nrow(child2.sd) == 0) child2.sd = child2.sd %>% add_row()
  
  lineage.id = lineages[match(divs$parent[1], unlist(sublines[, ii]))]
  mother.id  = sublines[match(divs$parent[1], unlist(sublines[, ii])), ii - 1]
  
  if(length(mother.id) == 0) mother.id = NA
  
  ### Find how long the mother has been alive for
  age = 0
  subline = sublines[match(divs$parent[1], unlist(sublines[, ii])), ] %>% unlist
  for (kk in ii:1) {
    if (ii == 1) {
      age = NA; break
    }
    if (kk == 1) break
    
    parent      = subline[kk - 1]
    no.siblings = sublines[which(unlist(sublines[, kk - 1]) == parent), ] %>% nrow
    
    if(no.siblings > 1) break
    if(no.siblings == 0) {age = NA; break}
    
    # If this is the only offspring
    age = age + timepoints[kk + 1] - timepoints[kk]
  }
  
  bind_cols(parent.qd, parent.sd, child1.qd, child1.sd, child2.qd, child2.qd) %>% 
    add_column(age     = age,            .before = 1) %>%
    add_column(lineage = lineage.id,     .before = 1) %>%
    add_column(t       = timepoints[ii], .before = 1)
}

get.division.data = function(sublines, mapping.data, segm.data, quant.data, timepoints, lineages, cl = NULL) {
  no.division.events = sapply(mapping.data, function(mapping.event)
    mapping.event$parent %>%
      duplicated         %>%
      sum)               %>%
    sum
  
  divisions = parLapply(cl, mapping.data, function(mapping.event) {
    mapping.event$parent %>% duplicated %>% which
  })
  
  clusterExport(cl, list = ls(envir = environment()), envir = environment())
  print(ls(envir = environment()))
  
  division.data = foreach(ii = 1:(length(timepoints) - 1), .combine = rbind) %:%
    foreach(jj = divisions[[ii]], .combine = rbind) %dopar% 
    get.division.event(ii, jj, sublines, mapping.data, segm.data, quant.data, 
                       timepoints, lineages, cl)
  
  division.data[division.data$m.vol == 0, "m.vol"] = NA
  division.data[division.data$n.vol == 0, "n.vol"] = NA
  division.data
}

get.cell.lineages = function(sublines, lineages, sep = "\\."){a
  split(as.data.frame(sublines), as.factor(lineages))
}

get.lineage.sublines.data = function(lineage, quant.data, timepoints) {
  no.timepoints = length(timepoints)
  no.q.cols     = ncol(quant.data[[1]]) + 1 # add time as first column
  
  members       = apply(lineage, 1, function(x) nna(x) %>% sum) %>% unlist %>% unname
  sublines.data = lapply(members, function(x)
    matrix(NA, nrow = x, ncol = no.q.cols)) # preallocate space
  
  for (subline in 1:length(members)) {
    counter = 1
    for (time in 1:length(quant.data)) {
      # If the member exists in the quantified data, add it
      qd.match = match(lineage[subline, time], quant.data[[time]]$id)
      if (nna(qd.match)) {
        sublines.data[[subline]][counter, ] =
          c(timepoints[time], unlist(quant.data[[time]][qd.match, ]))
        counter = counter + 1
      }
    }
    colnames(sublines.data[[subline]]) = c("t", colnames(quant.data[[1]]))
  }
  sublines.data
}

collapse.lineage = function(lineage) {
  do.call(rbind, lineage) %>% unique
}

