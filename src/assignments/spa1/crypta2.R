### Returns possible permutations of input
permutations <- function(n) {
   if (n == 1) 
   {
      return(matrix(1))
   } else 
   {
      sp <- permutations(n - 1)
      p <- nrow(sp)
      A <- matrix(nrow = n * p, ncol = n)
      for (i in 1:n) 
      {
         A[(i - 1) * p + 1:p, ] <- cbind(i, sp + (sp >= i))
      }
      return(A)
   }
}

### Get the number of row duplicates we will have
getUnique = function(x){
   if(x==0){
      return(1)
   }
   x * getUnique(x-1)
}

getFirstCharacters = function(string)
{
   string = gsub("[[:punct:]]", " ", string) # remove symbols
   string = gsub("\\s+"," ", string) # remove extra whitespaces
   #string = paste(string,collapse=" ") # 
   string = strsplit(string," ")[[1]]
   unique(sapply(string, function(x) substring(x,1,1)))
}

m = permutations(10) - 1
uniqueLetters = function(inString)
{
   inString = strsplit(inString, "")[[1]]
   inString = inString[!inString %in%  c(" ", "+", "*", "-", "=", "/", "&")]
   inString = unique(inString)     
}
inString = "AB * C = DE & DE + FG = HI"
#inString = "one + two + two + three + three = eleven" #"send + more = money"
letters = uniqueLetters(inString)
firstCharacters = getFirstCharacters(inString)
inString = strsplit(inString, "")[[1]]
colnames(m) = c(letters, rep("NA", 10 - length(letters)))
m = m[, -grep("NA", colnames(m))] # Slice the bread! (Take off edge.)
m = m[seq(1,nrow(m),getUnique(10 - length(letters))),] # Comb the desert! (Remove duplicates.)
for(ii in firstCharacters){
   m = m[-which(m[, ii] == 0),] # Kill off those w 0's at beginning.
}
inString = inString[!inString %in% c(" ")]
strings = matrix(inString, ncol = length(inString), nrow = nrow(m), byrow = TRUE)
colnames(strings) = strings[1,]

count = 1
for(ii in colnames(strings)){
   if(ii %in% colnames(m)){
      strings[, count] = m[,ii]
   }
   else if(ii == "=")
   {
      strings[, count] = rep("==", nrow(m))
   }
   count = count + 1
}

strings = matrix(do.call(paste0, as.data.frame(strings)))
output = sapply(strings, function(x) eval(parse(text=x)))
print(m[which(output == TRUE),])

#findSolution = function(str){
#  
#  
#}
