# metabolic-analysis.R by Rohan Maddamsetti

## None of these tests show anything interesting or significant.

## Condition analysis on metabolic genes in Nam database.
## Do Fisher’s exact test to test these following hypotheses:
## a)	Do LTEE substitutions fall in excess in “superessential” metabolic genes, per the Andreas Wagner paper?
## b)	Do hits fall primarily in “specialist” genes? Do Fisher’s exact test on specialist

ltee.mutations <- read.csv("../results/ltee_mutations_50K.csv")
core.summary <- read.csv("../results/core-summary.csv")

full.table <- merge(ltee.mutations, core.summary)
non.mutator.table <- full.table[full.table$mutator==FALSE,]

nam.table <- read.csv("../data/metabolic-data/Nam_Database.csv")
superessentials <- read.csv("../data/metabolic-data/Barve-Ecoli-Glucose-SuperessentialRxns.csv")

metabolic.data <- merge(superessentials, nam.table, all=TRUE)
metabolic.data$Superessential <- sapply(metabolic.data$Superessential, function(x) ifelse(is.na(x),0,1))
metabolic.data$Abbreviation <- NULL
metabolic.data <- unique(metabolic.data)
metabolic.data <-  metabolic.data[is.na(metabolic.data$Class)==FALSE,]


core.metabolic.data <- merge(core.summary, metabolic.data)

metabolic.substitution.data <- merge(core.metabolic.data, ltee.mutations)
#write.csv(core.metabolic.data,file="/Users/Rohandinho/Desktop/test.csv")

## Test hypothesis a). Compare superessential genes again the rest of the core.
superessential.rows <- metabolic.substitution.data[metabolic.substitution.data$Superessential==1 & metabolic.substitution.data$panortholog==TRUE,]
non.superessential.rows <- metabolic.substitution.data[metabolic.substitution.data$Superessential==0 & metabolic.substitution.data$panortholog==TRUE,]
##########################
superessential.length <- sum(core.metabolic.data[core.metabolic.data$Superessential==1,]$length)
non.superessential.length <- sum(core.metabolic.data[core.metabolic.data$Superessential==0,]$length)

superessential.dN <- sum(superessential.rows$dN)
non.superessential.dN <- sum(non.superessential.rows$dN)

superessential.contingency.table <- data.frame(length=c(superessential.length, non.superessential.length), dN=c(superessential.dN, non.superessential.dN))

fisher.test(superessential.contingency.table)

## p-value = 0.734, insignificant.

## Test hypothesis b)

spec.rows <- metabolic.substitution.data[metabolic.substitution.data$Class=="Specialist",]
gen.rows <-  metabolic.substitution.data[metabolic.substitution.data$Class=="Generalist",]

spec.length <- sum(core.metabolic.data[core.metabolic.data$Class=="Specialist",]$length)
gen.length <- sum(core.metabolic.data[core.metabolic.data$Class=="Generalist",]$length)

spec.dN <- sum(spec.rows$dN)
gen.dN <- sum(gen.rows$dN)

spec.gen.contingency.table <- data.frame(length=c(spec.length,gen.length),dN=c(spec.dN,gen.dN))

fisher.test(spec.gen.contingency.table)

## p-value = 0.5169, not significant.

## Test hypothesis c)? That specialists tend to be panorthologs?

## Did not do a formal test, seems irrelevant. Nothing striking about the data.

spec.rows <- core.metabolic.data[core.metabolic.data$Class=="Specialist",]
gen.rows <-  core.metabolic.data[core.metabolic.data$Class=="Generalist",]

spec.panorthologs <- spec.rows[spec.rows$panortholog==TRUE,]
spec.non.panorthologs <- spec.rows[spec.rows$panortholog==FALSE,]
gen.panorthologs <- gen.rows[gen.rows$panortholog==TRUE,]
gen.non.panorthologs <- gen.rows[gen.rows$panortholog==FALSE,]

spec.panortholog.count <- nrow(spec.panorthologs)
  spec.non.panortholog.count <- nrow(spec.non.panorthologs)
  gen.panortholog.count <- nrow(gen.panorthologs)
  gen.non.panortholog.count <- nrow(gen.non.panorthologs)


