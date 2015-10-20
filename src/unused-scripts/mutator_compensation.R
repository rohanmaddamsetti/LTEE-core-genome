## mutator_compensation.R by Rohan Maddamsetti.
## This R script looks for the signature of compensatory mutations in mutator
## strains in the LTEE. My reasoning is as follows.
##
## Hypothesis:
## 1) beneficial mutations at a locus tend to exclude other
## beneficial mutations due to stability concerns.
## 2) Deleterious mutations at a locus tend to open up new beneficial mutations
## to fix protein stability.
##
## First, I predict that the variance in dN per locus should be
## greater than the variance of dS per locus in mutators.
##
## Second, I predict that the variance in dN per locus in mutators
## should be greater than the variance in dN per locus in non-mutators.

##TODO: 1) do rigorous statistics. 2) make the code cleaner--i.e. more legible!

data <- read.csv("ltee_mutations.csv", header=T)
data$locus_tag<- as.character(data$locus_tag)
data$gene <- as.character(data$gene)

##calculate the variance in dN and dS per locus among mutators.
mutator.data <- data[data$mutator==TRUE,]

mutator.variances <- data.frame(locus=unique(mutator.data$locus_tag),dN.variance=rep(0,length(unique(mutator.data$locus_tag))), dS.variance=rep(0,length(unique(mutator.data$locus_tag))))

for (locus in unique(mutator.data$locus_tag)) {
    locus.data <- mutator.data[mutator.data$locus_tag==locus,]
    dN.locus.data <- locus.data[locus.data$dN > 0,]
    dS.locus.data <- locus.data[locus.data$dS > 0,]
    ## calculate the mean number of mutations per locus in mutators.
    mean.dN <- sum(locus.data$dN)/7
    mean.dS <- sum(locus.data$dS)/7
    ## now calculate the variance.
    genomes.hit.by.dN <- length(dN.locus.data$genome)
    genomes.hit.by.dS <- length(dS.locus.data$genome)
    dN.hit.terms <- ifelse(length(dN.locus.data$dN)>0, sum(sapply(dN.locus.data$dN, function (x) (x - mean.dN)^2)), 0)
    dN.no.hit.terms <- (7 - genomes.hit.by.dN)*((0-mean.dN)^2)
    dN.locus.variance = (dN.hit.terms+dN.no.hit.terms)/7

    dS.hit.terms <- ifelse(length(dS.locus.data$dS)>0, sum(sapply(dS.locus.data$dS, function (x) (x - mean.dS)^2)), 0)
    dS.no.hit.terms <- (7 - genomes.hit.by.dS)*((0-mean.dS)^2)
    dS.locus.variance = (dS.hit.terms+dS.no.hit.terms)/7
    mutator.variances[mutator.variances$locus==locus,]$dN.variance = dN.locus.variance
    mutator.variances[mutator.variances$locus==locus,]$dS.variance = dS.locus.variance
 }


##now calculate the variance in dN for non-mutators.
non.mutator.data <- data[data$mutator==FALSE,]

non.mutator.variances <- data.frame(locus=unique(non.mutator.data$locus_tag),dN.variance=rep(0,length(unique(non.mutator.data$locus_tag))))

for (locus in unique(non.mutator.data$locus_tag)) {
    locus.data <- non.mutator.data[non.mutator.data$locus_tag==locus,]
    dN.locus.data <- locus.data[locus.data$dN > 0,]

    ## calculate the mean number of mutations per locus in mutators.
    mean.dN <- sum(locus.data$dN)/6

    ## now calculate the variance.
    genomes.hit.by.dN <- length(dN.locus.data$genome)

    dN.hit.terms <- ifelse(length(dN.locus.data$dN)>0, sum(sapply(dN.locus.data$dN, function (x) (x - mean.dN)^2)), 0)
    dN.no.hit.terms <- (6 - genomes.hit.by.dN)*((0-mean.dN)^2)
    dN.locus.variance = (dN.hit.terms+dN.no.hit.terms)/6

    non.mutator.variances[non.mutator.variances$locus==locus,]$dN.variance = dN.locus.variance
   }


mean.mutator.variance <- sum(mutator.variances$dN.variance)/length(mutator.variances$dN.variance) 

mean.non.mutator.variance <- sum(non.mutator.variances$dN.variance)/length(non.mutator.variances$dN.variance)

## now only consider loci in common in the non.mutator and mutator variances.
common.loci <- non.mutator.variances[non.mutator.variances$locus %in% mutator.variances$locus,]$locus

common.mutator.variances <- mutator.variances[mutator.variances$locus %in% common.loci,]

common.non.mutator.variances <- non.mutator.variances[non.mutator.variances$locus %in% common.loci,]

sum(common.non.mutator.variances$dN)
sum(common.mutator.variances$dN)
sum(common.mutator.variances$dS)

## It appears that there's support for my hypothesis...
## Now, do some statistics to show that the difference in the variances is significant.

## As a first pass, assume that I can treat different loci as independent units
## with regard to my variance measurements across long-term lines.

## TODO: think of a better way to analyze these data, since mutations within a line
## are not necessarily independent. Rich suggested doing some jackknifing resampling,
## although I don't yet understand why this is an appropriate thing to do.


## Is the variance in dN different from variance in dS in mutators?
t.test(mutator.variances$dN.variance, mutator.variances$dS.variance,paired=TRUE)

## Is the variance in dN different from variance in dS in mutators,
## now only considering loci hit in common with non-mutators?
t.test(common.mutator.variances$dN,common.mutator.variances$dS, paired=TRUE)

## Is the variance in dN in mutators different from variance in dN in non-mutators?
t.test(common.mutator.variances$dN,common.non.mutator.variances$dN, paired=TRUE)

## Is the variance in dS in mutators different from variance in dN in non-mutators?
t.test(common.mutator.variances$dS,common.non.mutator.variances$dN, paired=TRUE)

library(ggplot2) ## for graphing these t-tests.

##plot for the first test:
