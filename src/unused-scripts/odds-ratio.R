# odds-ratio.R by Rohan Maddamsetti

## see code in core-contingency.R.
## This script sees if the odds ratio decreases or gets better over time.

library(ggplot2)

core.summary <- read.csv("../results/core-summary.csv")

ltee.mutations.50K <- read.csv("../results/ltee_mutations_50K.csv")
full.50K.table <- merge(ltee.mutations.50K, core.summary)
data.50K <- full.50K.table[full.50K.table$mutator==FALSE,]

ltee.mutations.40K <- read.csv("../results/ltee_mutations_40K.csv")
full.40K.table <- merge(ltee.mutations.40K, core.summary)
data.40K <- full.40K.table[full.40K.table$mutator==FALSE,]

ltee.mutations.30K <- read.csv("../results/ltee_mutations_30K.csv")
full.30K.table <- merge(ltee.mutations.30K, core.summary)
data.30K <- full.30K.table[full.30K.table$mutator==FALSE,]

ltee.mutations.20K <- read.csv("../results/ltee_mutations_20K.csv")
full.20K.table <- merge(ltee.mutations.20K, core.summary)
data.20K <- full.20K.table[full.20K.table$mutator==FALSE,]

ltee.mutations.10K <- read.csv("../results/ltee_mutations_10K.csv")
full.10K.table <- merge(ltee.mutations.10K, core.summary)
data.10K <- full.10K.table[full.10K.table$mutator==FALSE,]

ltee.mutations.5K <- read.csv("../results/ltee_mutations_5K.csv")
full.5K.table <- merge(ltee.mutations.5K, core.summary)
data.5K <- full.5K.table[full.5K.table$mutator==FALSE,]

ltee.mutations.2K <- read.csv("../results/ltee_mutations_2K.csv")
full.2K.table <- merge(ltee.mutations.2K, core.summary)
data.2K <- full.2K.table[full.2K.table$mutator==FALSE,]

odds.ratio <- function(data.table) {

core.rows <- data.table[data.table$panortholog==TRUE,]
noncore.rows <- data.table[data.table$panortholog==FALSE,]

core.length <- sum(core.summary[core.summary$panortholog==TRUE,]$length)
noncore.length <- sum(core.summary[core.summary$panortholog==FALSE,]$length)

core.dN <- sum(core.rows$dN)
noncore.dN <- sum(noncore.rows$dN)

core.dS <- sum(core.rows$dS)
noncore.dS <- sum(noncore.rows$dS)

contingency.table <- data.frame(length=c(core.length,noncore.length), dN=c(core.dN,noncore.dN))

## return the odds ratio
as.vector(unlist(fisher.test(contingency.table)[3]))
}

ratio.2K <- odds.ratio(data.2K)
ratio.5K <- odds.ratio(data.5K)
ratio.10K <- odds.ratio(data.10K)
ratio.20K <- odds.ratio(data.20K)
ratio.30K <- odds.ratio(data.30K)
ratio.40K <- odds.ratio(data.40K)
ratio.50K <- odds.ratio(data.50K)

odds.table <- data.frame(Generations=c(2000,5000,10000,20000,30000,40000,50000), Odds.Ratio=c(ratio.2K,ratio.5K,ratio.10K,ratio.20K,ratio.30K, ratio.40K, ratio.50K))

regression <- lm(data=odds.table, Odds.Ratio ~ Generations)

qplot(data=odds.table, x=Generations,y=Odds.Ratio,ylim=c(0.25,0.75))


## keep mutations in time interval independent,
## and ask if signal decreases over time

d.2K <- data.2K
d.2K$mutator <- NULL
d.2K$genome <- sapply(d.2K$genome, function (x) substr(x, 1,5))

d.5K <- data.5K
d.5K$mutator <- NULL
d.5K$genome <- sapply(d.5K$genome, function (x) substr(x, 1,5))

d.10K <- data.10K
d.10K$mutator <- NULL
d.10K$genome <- sapply(d.10K$genome, function (x) substr(x, 1,5))

d.20K <- data.20K
d.20K$mutator <- NULL
d.20K$genome <- sapply(d.20K$genome, function (x) substr(x, 1,5))

d.30K <- data.30K
d.30K$mutator <- NULL
d.30K$genome <- sapply(d.30K$genome, function (x) substr(x, 1,5))

d.40K <- data.40K
d.40K$mutator <- NULL
d.40K$genome <- sapply(d.40K$genome, function (x) substr(x, 1,5))

d.50K <- data.50K
d.50K$mutator <- NULL
d.50K$genome <- sapply(d.50K$genome, function (x) substr(x, 1,5))

## subtraction function from stack overflow.
df.subtract <- function (df1, df2) {
  return(df1[!duplicated(rbind(df2, df1))[-seq_len(nrow(df2))], ])
}

interval5K <- df.subtract(d.5K, d.2K)
interval10K <- df.subtract(d.10K, d.5K)
interval20K <- df.subtract(d.20K, d.10K)
interval30K <- df.subtract(d.30K, d.20K)
interval40K <- df.subtract(d.30K, d.40K)
interval50K <- df.subtract(d.50K, d.40K)

ratio1 <- odds.ratio(d.2K)
ratio2 <- odds.ratio(interval5K)
ratio3 <- odds.ratio(interval10K)
ratio4 <- odds.ratio(interval20K)
ratio5 <- odds.ratio(interval30K)
ratio6 <- odds.ratio(interval40K)
ratio7 <- odds.ratio(interval50K)

odds.table2 <- data.frame(Generations=c(2000,5000,10000,20000,30000,40000,50000), Odds.Ratio=c(ratio1,ratio2,ratio3,ratio4,ratio5, ratio6, ratio7))

regression2 <- lm(data=odds.table2, Odds.Ratio ~ Generations)

qplot(data=odds.table2, x=Generations,y=Odds.Ratio,ylim=c(0,2))
