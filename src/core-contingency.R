# core-contingency.R by Rohan Maddamsetti

## Tabulates dN in core and non-core genes, 
## in order to do a contingency test.
## The expected proportion of dN in the core is the sum of the length 
## of core genes/length of all genes.
## The observed proportion of dN in the core the the sum of dN 
## in the core/all dN.

ltee.mutations <- read.csv("../results/ltee_mutations_40K.csv")
core.summary <- read.csv("../results/core-summary.csv")

full.table <- merge(ltee.mutations, core.summary)

non.mutator.table <- full.table[full.table$mutator==FALSE,]

core.point.mutation.analysis <- function(data.table, do.dS=FALSE) {

  core.rows <- data.table[data.table$panortholog==TRUE,]
  noncore.rows <- data.table[data.table$panortholog==FALSE,]

  core.length <- sum(core.summary[core.summary$panortholog==TRUE,]$length)
  noncore.length <- sum(core.summary[core.summary$panortholog==FALSE,]$length)

  core.dN <- sum(core.rows$dN)
  noncore.dN <- sum(noncore.rows$dN)

  core.dS <- sum(core.rows$dS)
  noncore.dS <- sum(noncore.rows$dS)

  if (do.dS) {
    contingency.table <- data.frame(length=c(core.length,noncore.length), dS=c(core.dS,noncore.dS))
  } else {
    contingency.table <- data.frame(length=c(core.length,noncore.length), dN=c(core.dN,noncore.dN))
  }
  fisher.test(contingency.table)
}

## Test is highly significant for non.mutators.
core.point.mutation.analysis(non.mutator.table)

## Not enought dS data for non.mutators; only 25 mutations total.
core.point.mutation.analysis(non.mutator.table,do.dS=TRUE)

### Do same test on dS in mutators.

mutator.table <- full.table[full.table$mutator==TRUE,]

## Neither test is significant for mutators.
core.point.mutation.analysis(mutator.table, do.dS=FALSE)
core.point.mutation.analysis(mutator.table, do.dS=TRUE)

