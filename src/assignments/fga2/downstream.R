setwd("~/compbio/text/functional_genomics/assignment_2/")
load("fga2.RData")
library(survival)
#library(RColorBrewer)
cols = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
attach(Clinical)
Clinical = sapply(Clinical, function(x) gsub("-", ".", x, fixed=TRUE))
rownames(Clinical) = Clinical[,1]
Clinical = Clinical[order(Clinical[,1]), ]

########
# 1
########
par(mfrow=c(1,1))
cox.regression.model  = coxph(Surv(T, DeathBreast) ~ (grade + size  + er_status + IntClustMemb), data = Clinical)
survival.fit = survfit(cox.regression.model)
plot(survival.fit, lw = 3, bty = 'n', ylim = c(.4, 1))

cox.zph(cox.regression.model)
summary(cox.regression.model)
########
# 2     
########
par(mfrow = c(1,1))
grade.continuous = coxph(Surv(T, DeathBreast) ~ (grade         + log(size)  + er_status + IntClustMemb))
grade.discrete   = coxph(Surv(T, DeathBreast) ~ (factor(grade) + log(size)  + er_status + IntClustMemb))

# Exp(coeff) differs between the factors in the discrete model -- not a proportional hazard
#  plot(survfit(Surv(T, DeathBreast) ~ (grade)))
# lines(survfit(Surv(T, DeathBreast) ~ (factor(grade))), col = 3:3)
# plot(grade.continuous, col = cols[1:1], lw = 3, bty = 'n', ylim = c(.4, 1), conf.int=TRUE)
# lines(grade.discrete,  col = cols[2:2], lw = 3, bty = 'n', ylim = c(.4, 1))

# What are the p-values?
summary(grade.continuous)
summary(grade.discrete)

# Is the regression valid for all variables? Seems like it. (Continuous by default.)
cox.zph(grade.continuous)
cox.zph(grade.discrete)

########
# 3
########
er.fit = coxph(Surv(T, DeathBreast) ~ er_status)
plot(survfit(er.fit))
# Get the confidence interval
confidence = c(summary(er.fit)$conf.int[3:4])
# Conventionally, probabilities lower than .05 are considered significant 
# and researchers provide a 95% confidence interval for the hazard ratio, 
# e.g. derived from the standard deviation of the Cox-model regression 
# coefficient, i.e. {\displaystyle \beta } \beta. Statistically significant 
# hazard ratios cannot include unity (one) in their confidence intervals.

########
# 4
########
increases.risk = which(summary(grade.discrete)$coefficients[, 1] > 0)
decreases.risk = which(summary(grade.discrete)$coefficients[, 1] < 0)

par(mfrow = c(1, 2))
hist(size,      breaks = 30)
hist(log(size), breaks = 30)

########
# 5
########
par(mfrow = c(1,1))
plot(survfit(coxph(Surv(T, DeathBreast) ~ (factor(grade)), data = Clinical)))
lines(survfit(coxph(Surv(T, DeathBreast) ~ (size), data = Clinical)), col = 'red')

glm.fit = glm(DeathBreast~(factor(grade) + log(size) + er_status + IntClustMemb), family = binomial())

par(mfrow = c(1,2))
plot(grade.discrete$coefficients, ylim = c(-2,1))
plot(1:13, glm.fit$coefficients[1:13], ylim = c(-2,1))

###############################
########### PART 2 ############
###############################

genes = CNA[,1]
min   = min(unname(CNA[, -1]), na.rm = TRUE)
max   = max(unname(CNA[, -1]), na.rm = TRUE)

# Absolute value of 6, we have had an alteration of the genome producing extra copies of 
# our genes, with a total of 6 (as opposed to the other case, where we have a relative 
# grouping of {0, 1, 2+}). They are somewhat comparable as we expect an alteration in 
# copy number to give rise to a change in the gene expression. However, they aren't completely
# interchangable, as they both depend on the dependence between genes. For example, a change in
# expression for one gene can be driven by a copy number variation in another gene (i.e. a 
# repressor).

## ASCAT is in absolute values. ASCAT doesn't say how many copies of the two alleles we have --
## just the copies of the gene in total. Two genes with the same copy number can still be incomp.
## as they might have a significantly different setup. SNPs play an important role here.
## In Illumina, the actual products of the sum of the alleles are measured, which ought to be more
## comparable, as it is the OUTPUT of the setup. 

########
# B
########
# Answered above. 
# 
# Copy numbers should definietely regulate transcription, and therefore serve as a measure of 
# gene expression (how much TRANSCRIPT we have, not the amount of products). 

########
# C
########
sorted.CNA = CNA[order(CNA[, 1]), order(names(CNA))]
colnames(Exp)[1] = "Hugo_Symbol"
sorted.Exp = Exp[order(Exp[, 1]), order(names(Exp))]
cn.distribution = barplot(table(sapply(sorted.CNA[,-1], "[", 1)))

par(mfrow=c(1,2))
plot(sapply(sorted.Exp[, -1], "[", 1))
hist(sapply(sorted.Exp[, -1], "[", 1))

## Candidates -- change in gene expression AND copy number?
threshold = .4
candidate.genes.indices = c()
for(ii in 1:nrow(sorted.CNA))
    if(cor(as.numeric(sorted.Exp[ii,-1]), as.numeric(sorted.CNA[ii, -1]), use="na.or.complete") > threshold)
      candidate.genes.indices = c(candidate.genes.indices, ii)

# Now we have candidates
# TODO: Play around with threshold

########
# D
########

# Make expression values and copy number to be on the same scale? Might skew results otherwise.
#Create 10 equally size folds
data.exists = complete.cases(sorted.CNA[candidate.genes.indices, ]) & complete.cases(sorted.Exp[candidate.genes.indices, ])
candidate.genes.indices = candidate.genes.indices[data.exists]
preObj     = preProcess(sorted.CNA[,-1], method=c("center", "scale"))
sorted.CNA = predict(preObj, sorted.CNA)
preObj     = preProcess(sorted.Exp, method=c("center", "scale"))
sorted.Exp = predict(preObj, sorted.Exp)

rownames(sorted.Exp) = paste0("Exp.", sorted.Exp[,1])
sorted.Exp[,1] = paste0("Exp.", sorted.Exp[,1])

cluster.data = cbind(t(sorted.CNA[candidate.genes.indices, ]), t(sorted.Exp[candidate.genes.indices, ]))
colnames(cluster.data) = cluster.data[1,]
cluster.data = cluster.data[-1,]
r.names = rownames(cluster.data)
cluster.data = apply(cluster.data, 2, as.numeric)
rownames(cluster.data) = r.names

#########
# Make matrix with GENES as variables * 2, i.e. 200 variables
# MB.data are ROWS
#########
#Perform 10 fold cross validation
folds = cut(seq(1, nrow(cluster.data)), breaks = 10, labels = FALSE)
targets = as.integer(Clinical[, 5])
for(i in 1:10){
  test.indices   = which(folds == i, arr.ind = TRUE)
  test.data      = cluster.data[test.indices, ]
  train.data     = cluster.data[-test.indices, ]
  test.targets   = targets[test.indices]
  train.targets  = targets[-test.indices]
}
## Size question: split into thirds and compare? 
## Which is a better predictor of death events?
## Do cox model fitted with either variable and 
## find which is the best predictor compared to
## the actual case.

library(nnet)
rownames(Clinical)   = gsub("-", ".", Clinical[,1], fixed=TRUE)
train.targets        = sapply(rownames(train.data), function(x){as.integer(Clinical[x, 5])})
test.targets         = sapply(rownames(test.data),  function(x){as.integer(Clinical[x, 5])})
trained.network      = nnet(as.factor(train.targets)~., data = train.data, size = 20, MaxNWts = 20000)
test.classification  = predict(trained.network, test.data,  type = "class")
train.classification = predict(trained.network, train.data, type = "class")

library(e1071)
tmodel = tune.nnet(as.factor(train.targets)~train.data, data = as.data.frame(train.data), size = 5:20, MaxNWts = 20000)

# Statistics
train.accuracy = sum(train.classification == train.targets) / length(train.targets)
test.accuracy  = sum(test.classification  == test.targets)  / length(test.targets)


###

# # library(neuralnet)
# # nn <- neuralnet(
#   IntClustMemb ~ (grade + size + er_status + T + DeathBreast),
#   data = train.data, hidden = 2, err.fct = "sse",
#   linear.output = FALSE)



###
### BOTH ARE COMPARABLE

## Values by ASCAT always comparable as they are absolute numbers.
## The expression data has to be normalised first, as you compare between arrays.