sink("curved_squares_out.txt")

### Take off axes
pdf("curved_squares.pdf", width = 3, height = 3)
par(bty = 'n',xaxt = 'n',yaxt = 'n')
par(mfrow = c(2,2), mar = c(0,0,0,0) + 0.0)

### Define the functions, upper-to-lower, lower-to-upper
x = seq(0, 1, 10  ^ -4)
lowerUTL = function(x){-sqrt(1 - (x-1)^2) + 1}
upperUTL = function(x){sqrt(1 - x^2)}
lowerLTU = function(x){-sqrt(1 - x^2) + 1}
upperLTU = function(x){sqrt(1 - (x-1)^2)}

####### FIGURE 1
### Draw the arcs
plot(x, lowerUTL(x), type='l', ylab='', xlab='')
lines(x, lowerLTU(x))
lines(x, upperLTU(x))
lines(x, upperUTL(x))
rect(0,0,1,1)

### Indices we need
intersectionIndex = which(abs(lowerUTL(x) - upperLTU(x)) < 10 ^ -4)
midpointIndex = length(x)/2

####### FIGURE 2
### Draw the arcs
plot(x, lowerUTL(x), type='l', ylab='', xlab='')
lines(x, lowerLTU(x))
lines(x, upperLTU(x))
lines(x, upperUTL(x))
rect(0,0,1,1)

### Le polygon
polygon(c(x[seq(intersectionIndex, length(x) - intersectionIndex)], 
  rev(x[seq(intersectionIndex,length(x) - intersectionIndex)])),
  c(upperLTU(x[seq(intersectionIndex, length(x) - midpointIndex)]), 
  upperUTL(x[seq(midpointIndex,length(x) - intersectionIndex)]), 
  rev(lowerLTU(x[seq(midpointIndex,length(x) - intersectionIndex)])),
  rev(lowerUTL(x[seq(intersectionIndex, midpointIndex)]))), col='purple')

####### FIGURE 3
### Draw the arcs
plot(x, lowerUTL(x), type = 'l', ylab = '', xlab = '')
lines(x, lowerLTU(x))
lines(x, upperLTU(x))
lines(x, upperUTL(x))
rect(0,0,1,1)

### Le polygon
polygon(c(x[seq(intersectionIndex, length(x) - intersectionIndex)], 
  rev(x[seq(intersectionIndex,length(x) - intersectionIndex)])),
  c(upperLTU(x[seq(intersectionIndex, length(x) - midpointIndex)]), 
  upperUTL(x[seq(midpointIndex,length(x) - intersectionIndex)]), 
  rev(lowerLTU(x[seq(midpointIndex,length(x) - intersectionIndex)])),
  rev(lowerUTL(x[seq(intersectionIndex, midpointIndex)]))), col='purple')

### Le cube
y1 =  x + upperLTU(x[intersectionIndex]) - x[intersectionIndex]
y2 = -x + upperLTU(x[intersectionIndex]) - x[intersectionIndex] + 1
y3 =  x + lowerLTU(.5) - .5
y4 = -x + lowerLTU(.5) + .5

area.square = (.5 - x[intersectionIndex])^2 + 
  (upperLTU(.5) - upperLTU(x[intersectionIndex]))^2
cat("The area of the square contained is:\t", area.square, "\n")

polygon(c(x[seq(intersectionIndex, length(x) - intersectionIndex)], 
  rev(x[seq(intersectionIndex,length(x) - intersectionIndex)])),
  c(y1[seq(intersectionIndex, length(x) - midpointIndex)], 
  y2[seq(midpointIndex,length(x) - intersectionIndex)], 
  rev(y3[seq(midpointIndex,length(x) - intersectionIndex)]),
  rev(y4[seq(intersectionIndex,midpointIndex)])), col='gray')

###### FIGURE 4
### Draw the arcs
plot(x, lowerUTL(x), type='l', ylab='', xlab='')
lines(x, lowerLTU(x))
lines(x, upperLTU(x))
lines(x, upperUTL(x))
rect(0,0,1,1)

### Le polygon
polygon(c(x[seq(intersectionIndex, length(x) - intersectionIndex)], 
  rev(x[seq(intersectionIndex,length(x) - intersectionIndex)])),
  c(upperLTU(x[seq(intersectionIndex, length(x) - midpointIndex)]), 
  upperUTL(x[seq(midpointIndex,length(x) - intersectionIndex)]), 
  rev(lowerLTU(x[seq(midpointIndex,length(x) - intersectionIndex)])),
  rev(lowerUTL(x[seq(intersectionIndex, midpointIndex)]))), col='purple')

### Naive Monte-Carlo (box around purple shape)
### Note: Our problem is symmetric, so we can limit ourselves to a quartile.
nrTries = 1e7
xHit = runif(nrTries, .5, x[length(x) - intersectionIndex])
yHit = runif(nrTries, x[intersectionIndex], .5)
hits = which(yHit > lowerLTU(xHit)) 
misses = which(yHit <= lowerLTU(xHit))

### Plot every 100th point.
points(xHit[hits[seq(1,length(hits),10000)]], yHit[hits[seq(1,length(hits),10000)]], 
  type = 'p', cex = .01, col = 'forestgreen')
points(xHit[misses[seq(1,length(misses),10000)]], yHit[misses[seq(1,length(misses), 10000)]], 
  type = 'p', cex = .01, col = 'tomato')
rect(.5,x[intersectionIndex], x[length(x) - intersectionIndex], .5)

area.estimate = length(hits)/nrTries * 
  (x[length(x) - intersectionIndex] - x[intersectionIndex])^2
cat("The estimated area of the shape is:\t", area.estimate,"\n")

### Kill off graphical output
dev.off()
