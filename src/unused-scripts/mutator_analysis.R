##mutator_analysis.R by Rohan Maddamsetti.

library(ggplot2)

##calculate the empirical distribution of genes with n=k nonsynonymous mutations.
dN.distribution <- function (mutator.hit.table) {
  dNs <- mutator.hit.table$nonsynonymous
  counts <- rep(0,5)
  for (i in dNs)
    counts[i] = counts[i] + 1
  return(counts)
}

##subset the empirical distribution of genes with n=k nonsynonymous mutations on
##genes which have nonsynonymous mutations in the nonmutators.
subset.dN.distribution <- function(mutator.hit.table, nonmutator.hit.table) {
  relevant.mutator.hits <- subset(mutator.hit.table, locus_tag %in% nonmutator.hit.table$locus_tag)
  return(dN.distribution(relevant.mutator.hits))
}

##takes a vector of probabilities of falling in a particular bin,
##the total number of mutations to drop, and the number of runs.
null.distribution <- function(the.bins, total.mut, runs=10) {

  counts <- rep(0,5)

  for (run in 1:runs) {
    ##make a empty vector of bins.
    output.bins <- rep(0, length(the.bins))
    ##Convert the bins into  a proper probability distribution.
    bin.pdf = the.bins
    for (i in 1:length(the.bins))
      bin.pdf[i] = sum(the.bins[1:i])

    ##generate a vector of random numbers between 0 and 1.
    rands <- runif(total.mut)
    ## for each deviate, add 1 to the proper slot in output.bins.
    for (rand in rands) {
      i = 1
      while (rand > bin.pdf[i])
        i = i + 1
      
      output.bins[i] = output.bins[i] + 1
    }
    ## now make a distribution of bins with n=k mutations.
    for (i in output.bins)
      counts[i] <- counts[i] + 1
    
  }
  return (sapply(counts, function (x) x/sum(counts)))
}
  
a1.hits <- read.table("../data/mutator_40K_diffs/A-1hits.txt", header=T)
a3.hits <- read.table("../data/mutator_40K_diffs/A-3hits.txt", header=T)
non.mutator.hits <- read.table("../data/non-mutator_40K_diffs/genehits.tab", header=T)

a1.dist <- dN.distribution(a1.hits)
normalized.a1.dist <- a1.dist/sum(a1.dist)
a1.dist.subset <- subset.dN.distribution(a1.hits, non.mutator.hits)
normalized.a1.dist.subset <- a1.dist.subset/sum(a1.dist.subset)
a3.dist <- dN.distribution(a3.hits)
normalized.a3.dist <- a3.dist/sum(a3.dist)
a3.dist.subset <- subset.dN.distribution(a3.hits, non.mutator.hits)
normalized.a3.dist.subset <- a3.dist.subset/sum(a3.dist.subset)

##The probabilities of getting a mutation in any particular protein-coding region.
##This datafile was produced by simgenome.py, and it allows mutations to be dropped among genes,
##weighted by the length of the locus.
gene.bins <- as.vector(unlist(read.table("REL606_genebin_probabilities.txt", header=F)))
a1.null.dist <- null.distribution(gene.bins, sum(a1.hits$nonsynonymous), runs=10)
a3.null.dist <- null.distribution(gene.bins, sum(a3.hits$nonsynonymous), runs=10)

hit.names <- 1:5
a1.plot.data <- data.frame(mutations=hit.names, null.dist=a1.null.dist, dist=normalized.a1.dist, dist.subset=normalized.a1.dist.subset)
a3.plot.data <- data.frame(mutations=hit.names, null.dist=a3.null.dist, dist=normalized.a3.dist, dist.subset=normalized.a3.dist.subset)
##plot histograms on top of each other.
quartz()
p1 <- ggplot(a1.plot.data) +
  geom_area(aes(mutations, null.dist), fill="red", alpha = 0.2) + 
  geom_area(aes(mutations, dist), fill="blue", alpha = 0.2) +
  geom_area(aes(mutations, dist.subset), fill= "green", alpha = 0.2)
p3 <- ggplot(a3.plot.data) +
  geom_area(aes(mutations, null.dist), fill="red", alpha = 0.2) + 
  geom_area(aes(mutations, dist), fill="blue", alpha = 0.2) +
  geom_area(aes(mutations, dist.subset), fill= "green", alpha = 0.2)

p1
