
all.data = output.data
cell.lines
cell.lines.data
interesting = (1:ncol(all.data))[-c(1:5, 8:10)] # Remove xyz and indices

# Lineages ordered by nMembers
par(mfrow=c(4,3))
for(ii in 1:12) plot(all.data[all.data[,"lineage.id"] == ii, "age"], ylim=c(0,20), main = "age")
for(ii in 1:12) plot(cell.lines.data[[ii]][,"Time"], cell.lines.data[[ii]][,"Mean.cell.intensity"], ylim=c(0,160), main = "expr")


for(ii in interesting) plot(all.data[all.data[,"lineage.id"] == ii, "age"], ylim=c(0,20), main = colnames(all.data)[ii])
