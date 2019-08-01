sink("cryptarithms_out.txt")

#####################
### FUNCTIONS
#####################

### Returns possible permutations of input
permutations <- function(n) {
  if (n == 1)
  {return(matrix(1))} else
  {
    sp <- permutations(n - 1)
    p <- nrow(sp)
    A <- matrix(nrow = n * p, ncol = n)
    for (i in 1:n)
    {
      A[(i - 1) * p + 1:p,] <- cbind(i, sp + (sp >= i))
    }
    return(A)
  }
}
getUnique = function(x)
{
  if (x == 0) {
    return(1)
  }
  x * getUnique(x - 1)
}

### Which letters can't be 0?
getFirstCharacters = function(string)
{
  string = gsub("[[:punct:]]", " ", string) # remove some symbols
  string = gsub("\\s+"," ", string) # remove extra whitespaces
  string = strsplit(string," ")[[1]]
  unique(sapply(string, function(x)
    substring(x,1,1)))
}

### Which are the unique letters?
uniqueLetters = function(inString)
{
  inString = strsplit(inString, "")[[1]]
  inString = inString[!inString %in%  c(" ", "+", "*", "-", "=", "/", "&")]
  inString = unique(inString)
}

#####################
### MAIN
#####################

### Takes in string and evaluates the corresponding cryptarithm.
### Note: Several statements are separated by '&'.
### Note: Only handles base 10 at the moment.
### Note: Very inefficient; struggles with larger strings (=> larger matrices)
crypta = function(string.in) {

  ### Get the number of row duplicates we will have
  m = permutations(10) - 1
  inString = string.in
  letters = uniqueLetters(inString)
  firstCharacters = getFirstCharacters(inString)
  inString = strsplit(inString, "")[[1]]
  colnames(m) = c(letters, rep("NA", 10 - length(letters)))
  if(length(letters) != 10)
  {
    m = m[, -grep("NA", colnames(m))] # Take off edge.
  }
  m = m[seq(1,nrow(m),getUnique(10 - length(letters))),] # Remove duplicates.
  for (ii in firstCharacters) {
    m = m[-which(m[, ii] == 0),] # Kill off those w 0's at beginning.
  }
  inString = inString[!inString %in% c(" ")] # Take out essentials
  strings = matrix(
    inString, ncol = length(inString), nrow = nrow(m), byrow = TRUE
  )
  colnames(strings) = strings[1,]
  
  ### Adjust our columns correspondingly.
  count = 1
  for (ii in colnames(strings)) {
    if (ii %in% colnames(m)) {
      strings[, count] = m[,ii]
    }
    else if (ii == "=")
    {
      strings[, count] = rep("==", nrow(m))
    }
    count = count + 1
  }
  
  ### Add strings back together
  strings = matrix(do.call(paste0, as.data.frame(strings)))
  output = sapply(strings, function(x)
    eval(parse(text = x)))
  print(m[which(output == TRUE),])
  cat("#################################","\n")
}

### Run program with our different setups.
crypta("AB * C = DE & DE + FG = HI")
crypta("send + more = money")
crypta("snow + rain = sleet")
crypta("one + two + two + three + three = eleven")