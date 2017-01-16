##G-core-variance.R by Rohan Maddamsetti

## I ask whether core genes have higher variance in G-score (more positive and purifying selection)
## in mutators. I found this wasn't the case--any difference seem to be driven by outliers, the distributions
## actually look pretty similar.

library(ggplot2)
library(dplyr)
library(ggrepel)

## Analysis using G-score data (exactly same results as with 40K single clone data).

G.score.data <- read.csv("../data/Tenaillon-data/nature18959-s2.csv")
core.summary <- read.csv("../results/core-summary.csv")
KEIO.data <- read.csv("../data/KEIOdata.csv")

## turns out that start position is a good unique key for genes.
core.non.mutator.data <- merge(G.score.data, core.summary, all=TRUE)

core.proteome <- read.csv("../Data/Palsson-data/pnas.1501384112.sd01.csv") %>% rename(Gene.name=name) %>% select(Gene.name,gene)

core.proteome2 <- merge(core.summary,core.proteome)
tenaillon.mutator.data <- read.csv("../data/Tenaillon-data/nature18959-s3.csv")
core.mutator.data <- merge(tenaillon.mutator.data,core.summary,all=TRUE)

proteome.mutator.data <- merge(core.proteome,core.mutator.data) %>% select(Gene.name,Gene.order,Observed.nonsynonymous.mutation,Expected.nonsynonymous.mutation,G.score,panortholog)

#### take a look at negative G score genes in mutators. Are these generally core genes?
#### it's 753 core negatives to 596 non-core negatives.
negatives <- filter(core.mutator.data,G.score<0)
core.negatives <- filter(negatives,panortholog==TRUE)
noncore.negatives <- filter(negatives,panortholog==FALSE)

length(negatives$G.score)
core.negatives <- filter(negatives,panortholog==TRUE)
noncore.negatives <- filter(negatives,panortholog==FALSE)
length(noncore.negatives$G.score)
length(core.negatives$G.score)
mean(core.negatives$G.score)
mean(noncore.negatives$G.score)


## difference in variance in G.score in panorthologs in mutators? Apparently not.
dist.plot <- ggplot(core.mutator.data,aes(x=G.score,fill=panortholog)) + geom_histogram(bins=100)

var(filter(core.mutator.data,panortholog==TRUE)$G.score)
var(filter(core.mutator.data,panortholog==FALSE)$G.score)

filtered.core.mutator.data <- filter(core.mutator.data,G.score != 0)

dist.plot2 <- ggplot(filtered.core.mutator.data,aes(x=G.score,fill=panortholog)) + geom_histogram(bins=100)

var(filter(filtered.core.mutator.data,panortholog==TRUE)$G.score)
var(filter(filtered.core.mutator.data,panortholog==FALSE)$G.score)

############ looking at the core proteome in the Palsson dataset.

## plot distributions of G-score to check for approximate normality.
prot.dist.plot <- ggplot(proteome.mutator.data,aes(x=G.score,fill=panortholog)) + geom_histogram(bins=100)

var(filter(proteome.mutator.data,panortholog==TRUE)$G.score)

var(filter(proteome.mutator.data,panortholog==FALSE)$G.score)

## quick check for non-mutators.
proteome.non.mutator.data <- merge(core.proteome,core.non.mutator.data) %>% select(Gene.name,Gene.order,Observed.nonsynonymous.mutation,Expected.nonsynonymous.mutation,G.score,panortholog)

prot.dist.plot2 <- ggplot(proteome.non.mutator.data,aes(x=G.score,fill=panortholog)) + geom_histogram(bins=100)

var(filter(proteome.non.mutator.data,panortholog==TRUE)$G.score)

var(filter(proteome.non.mutator.data,panortholog==FALSE)$G.score)
