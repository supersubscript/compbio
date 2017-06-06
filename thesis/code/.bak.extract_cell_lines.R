#############################################################################
# Extract cell lines
#############################################################################
### Read in first generation and construct the cell line list
# first.track.data = all.track.data[[1]]
# cell.lines = lapply(lapply(first.track.data, function(x) x[1]), function(x) list(x)) # These are all the original cell lines
# for (ii in 1:length(first.track.data)) {
#   for (jj in 1:length(cell.lines)) {
#     is.in.cell.line = first.track.data[[ii]][1] %in% cell.lines[[jj]][[1]]
#     if (is.in.cell.line) {
#       cell.lines[[jj]] = c(cell.lines[[jj]], list(first.track.data[[ii]][-1])) # Add daughter cells
#       break
#     }
#   }
# }

# ### Now do all but the first
# for (time in 2:no.track.files) {
#   this.track.data = all.track.data[[time]] 
#   
#   # Now sort out the cell lines
#   for (ii in 1:length(this.track.data)) {
#     added = FALSE
#     for (jj in 1:length(cell.lines)) {
#       # Is the mother cell in the last time point?
#       is.in.cell.line = this.track.data[[ii]][1] %in% cell.lines[[jj]][[time]]
#       if (is.in.cell.line) {
#         next.timepoint.exists = length(cell.lines[[jj]]) == time + 1
#         if (next.timepoint.exists)
#           cell.lines[[jj]][[time + 1]] = add.to.list(cell.lines[[jj]][[time + 1]], this.track.data[[ii]][-1]) # Add daughter cells to next time frame
#         else
#           cell.lines[[jj]] = add.to.list(cell.lines[[jj]], this.track.data[[ii]][-1]) # Add daughter cells to next time frame
#         added = TRUE
#         break
#       }
#     }
#     # If we didn't find it
#     if (!added) {
#       new.cell.line = list(c(sapply(1:(time - 1), 
#                                     function(x) list(NA)), list(this.track.data[[ii]][1])))
#       if (length(this.track.data[[ii]]) > 1)
#         new.cell.line[[1]] = c(new.cell.line[[1]], 
#                                list(this.track.data[[ii]][-1]))
#       cell.lines = c(cell.lines, new.cell.line)
#     }
#   }
#   # Fill up space for the ones that didn't have daughters
#   cell.lines[sapply(cell.lines, length) == time] = lapply(cell.lines[sapply(cell.lines, length) == time], function(x) c(x, list(NA)))
# }
# 
### Now do all but the first
# for (time in 2:no.track.files) {
#   track.data = all.track.data[[time]] 
#   
#   # Now sort out the cell lines
#   for (ii in 1:length(track.data)) {
#     added = FALSE
#     for (jj in 1:length(cell.lines)) {
#       # If mother is in last time point
#       if (track.data[[ii]][1] %in% cell.lines[[jj]][[time]]) {
#         # If next timeframe doesn't already exist, add it
#         if (length(cell.lines[[jj]]) != time + 1)
#           cell.lines[[jj]] = c(cell.lines[[jj]], list(track.data[[ii]][-1])) # Add daughter cells to next time frame
#         else
#           # If it exists, add to it
#           cell.lines[[jj]][[time + 1]] = c(cell.lines[[jj]][[time + 1]], track.data[[ii]][-1]) # Add daughter cells to next time frame
#         added = TRUE
#         break
#       }
#     }
#     # If we didn't find it
#     if (!added) {
#       new.cell.line = list(c(sapply(1:(time - 1), 
#                                     function(x) list(NA)), list(track.data[[ii]][1])))
#       if (length(track.data[[ii]]) > 1)
#         new.cell.line[[1]] = c(new.cell.line[[1]], 
#                                list(track.data[[ii]][-1]))
#       cell.lines = c(cell.lines, new.cell.line)
#     }
#   }
#   # Fill up space for the ones that didn't have daughters
#   cell.lines[sapply(cell.lines, length) == time] = lapply(cell.lines[sapply(cell.lines, length) == time], function(x) c(x, list(NA)))
# }

# Make it a data frame of lists instead
cell.lines = t(sapply(cell.lines, function(x) x))
cell.lines = cell.lines[order(-apply(cell.lines, 1, function(x) sum(!is.na(x)))), ]

# Make all integer(0) NAs
for(ii in 1:ncol(cell.lines)) cell.lines[which(!(sapply(cell.lines[,ii], length))), ii] = list(NA)
