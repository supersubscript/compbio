rm(list = ls())
sink("word_processing_out.txt")
### 1 Number of words
words = read.table(
  "../../../src/sp_assignments/assignment_1/usr-share-dict-words")
words = unlist(words)
#words = tolower(words)
#words = unique(words)
cat("Number of words:\t", length(words), "\n")

### 2 Ascii
words.ascii = words[-grep("[^\x21-\x7e]", words)]
words.ascii = as.character(words.ascii)
cat("Number of ascii words:\t", length(words) - length(words.ascii), 
    "\n")

### 3 Vowels
noVowelWords = words.ascii[-grep("[aeiouy]", words.ascii, 
  ignore.case = TRUE)]
cat("Number of non-v words:\t", length(noVowelWords), "\n")
write(noVowelWords, "word_processing_noVowelWords.txt", sep="\t", ncolumns = 6)

### 4 Palindromes
reverseString = function(x) paste(substring(x, nchar(x):1, 
  nchar(x):1), collapse = "")
palindromes = words.ascii[which(sapply(words.ascii, 
  function(x) tolower(x) == tolower(reverseString(x))) == TRUE)]
write(palindromes, "word_processing_palindromes.txt", sep="\t", ncolumns = 6)

### 5 Apostrophes
no.apostrophes = grep("[']", words.ascii, value = TRUE, 
  invert = TRUE)
cat("Number of words with apostrophes:\t", length(words.ascii) 
    - length(no.apostrophes), "\n")

### 6 Probability distribution
letter.frequency = paste(no.apostrophes, collapse = "")
letter.frequency = summary(factor(tolower(strsplit(
  letter.frequency, "")[[1]])))
letter.frequency = letter.frequency / sum(letter.frequency)
pdf("word_processing_letter_frequency.pdf", width = 5, height = 5)
par(mfrow = c(1, 1), las = 2)
barplot(letter.frequency * 100, xlab = 'Letter', 
  names.arg = letters[1:26], ylab="Frequency (%)")
x = dev.off()

### 7 Bigram
bigram.matrix = matrix(0, length(letters) + 1, length(letters) + 1) 
bigram.data   = tolower(paste(c("", no.apostrophes, ""), 
  collapse = " "))
colnames(bigram.matrix) = c(letters, " ")
rownames(bigram.matrix) = c(letters, " ")
pairs = substring(bigram.data, 1:(nchar(bigram.data) - 1), 
  2:nchar(bigram.data))
pairs.data = table(pairs)

cat("The most frequent bigram is \'", 
  names(pairs.data)[which(pairs.data == max(pairs.data))],
  "\' with a frequency of ", 
  pairs.data[which(pairs.data == max(pairs.data))] 
  / sum(pairs.data), ".", sep = "")

for(ii in names(pairs.data))
{
  chars = strsplit(ii,"")[[1]]
  bigram.matrix[chars[1], chars[2]] = 
    bigram.matrix[chars[1], chars[2]]  + pairs.data[ii]
}

pdf("word_processing_bigram.pdf", width = 10, height = 10)
par(mfrow=c(1,1))
nrs = as.numeric(factor(rownames(bigram.matrix)))
plot(0, 0, type='n', xlim = c(1,27), ylim = c(1,27), 
     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
for(ii in nrs)
{
  for(jj in nrs)
  {
    points(ii, jj, cex = 0.22 * max(0, log(bigram.matrix[ii,jj])), 
           pch = 15, col='forestgreen')
  }  
}
axis(1, at=1:27, labels=c(letters[1:26], '_'),tck= 0,cex.axis = .7)
axis(2, at=1:27, labels=c(letters[1:26], '_'),tck= 0,cex.axis = .7)
x = dev.off()

