# core-binomial-test.R by Rohan Maddamsetti

## Tabulates dN in core and non-core genes,
## The binomial test compares observed counts of
## dN in core and in non-core genes to expected proportions
## based on the total length of core and non-core genes.

library(dplyr)

core.point.mutation.analysis <- function(data.table, do.dS=FALSE) {

    core.rows <- data.table[data.table$panortholog==TRUE,]
    noncore.rows <- data.table[data.table$panortholog==FALSE,]

    core.length <- sum(core.summary[core.summary$panortholog==TRUE,]$length)
    noncore.length <- sum(core.summary[core.summary$panortholog==FALSE,]$length)
    expected <- c(core.length,noncore.length)/(core.length+noncore.length)
    expected.p <- core.length/(core.length+noncore.length)

    core.dN <- sum(core.rows$dN)
    noncore.dN <- sum(noncore.rows$dN)

    core.dS <- sum(core.rows$dS)
    noncore.dS <- sum(noncore.rows$dS)

    if (do.dS) {
        print("core dS:")
        print(core.dS)
        print("noncore.dS")
        print(noncore.dS)
        binom.test(c(core.dS,noncore.dS),p=expected.p,alternative="two.sided")
    } else {
        print("core dN:")
        print(core.dN)
        print("noncore.dN")
        print(noncore.dN)
        binom.test(c(core.dN,noncore.dN),p=expected.p,alternative="two.sided")
       }
}

ltee.mutations <- read.csv("../results/ltee_mutations.csv")
core.summary <- read.csv("../results/core-summary.csv")

full.table <- merge(ltee.mutations, core.summary)

non.mutator.table <- full.table[full.table$mutator==FALSE,]
## Test is highly significant for non.mutators.
core.point.mutation.analysis(non.mutator.table,do.dS=FALSE)
## Not much dS data for non.mutators; only 25 mutations total.
core.point.mutation.analysis(non.mutator.table,do.dS=TRUE)

### Do same test on dS in mutators.
mutator.table <- full.table[full.table$mutator==TRUE,]
## Neither test is significant for mutators.
core.point.mutation.analysis(mutator.table, do.dS=FALSE)
core.point.mutation.analysis(mutator.table, do.dS=TRUE)

# Additional non-parametric t-test on G-scores in core genes and non-core genes.
G.score.data <- read.csv("../data/Tenaillon-data/nature18959-s2.csv")
test <- merge(core.summary,G.score.data,all=TRUE)

test <- merge(core.summary,G.score.data,all=TRUE)
wilcox.test(test[test$panortholog==TRUE,]$G.score,test[test$panortholog==FALSE,]$G.score)


## double-check synonymous mutation numbers.
#non.mutator.table %>% arrange(genome) %>% filter(dS > 0)
