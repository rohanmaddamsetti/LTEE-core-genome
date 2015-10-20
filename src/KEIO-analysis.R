# KEIO-analysis.R by Rohan Maddamsetti

# Do a t-test on KEIO essentiality and growth defect scores comparing
# panorthologs and non-panorthologs.

library(ggplot2)

core.summary <- read.csv("../results/core-summary.csv")
KEIO.data <- read.csv("../data/KEIOdata.csv")

functional.data <- merge(core.summary, KEIO.data)

panortholog.data <- functional.data[functional.data$panortholog==TRUE,]

ltee.40K.mutations <- read.csv("../results/ltee_mutations_40K.csv")
non.mutator.data <- subset(ltee.40K.mutations,ltee.40K.mutations$mutator==FALSE)

## count the number of lineages with at least one dN per locus.  
non.mutator.data$lineage <- sapply(non.mutator.data$genome,function(x) substr(x,1,5))
hits.data <- non.mutator.data
hits.data$genome <- NULL
hits.data$mutator <- NULL
hits.data$dS <- NULL
hits.data <- hits.data[hits.data$dN>0,]
  
locus <- unique(hits.data$locus_tag)
hits <- rep(0,length(locus))
for (i in 1:length(locus)) {
  l = locus[i]
  cur.hits <- subset(hits.data, locus_tag==l)
  hit.lineages <- length(unique(cur.hits$lineage))
  hits[i] <- hit.lineages
}
    
final.hits <- data.frame(locus_tag=locus,lineages=hits)
final.data <- merge(final.hits, panortholog.data,all=TRUE)
final.data$lineages <- sapply(final.data$lineages, function (x) ifelse(is.na(x), 0, x))

## p-value = 0.6343. No difference in essentiality between panorthologs that mutated in the LTEE and those that did not.
t.test(final.data[final.data$lineages>0,]$EssentialityScore, final.data[final.data$lineages==0,]$EssentialityScore)

## p-value = 0.17. No difference in LB growth between panorthologs that mutated in the LTEE and those that did not.
t.test(final.data[final.data$lineages>0,]$LB_22hr_OD600, final.data[final.data$lineages==0,]$LB_22hr_OD600)

## p-value = 0.89 No difference in MOPS 24 hour growth between panorthologs that mutated in the LTEE and those that did not.
t.test(final.data[final.data$lineages>0,]$MOPS_24hr_OD600, final.data[final.data$lineages==0,]$MOPS_24hr_OD600)

## p-value = 0.4247. No difference in MOPS 48 hour growth between panorthologs that mutated in the LTEE and those that did not.
t.test(final.data[final.data$lineages>0,]$MOPS_48hr_OD600, final.data[final.data$lineages==0,]$MOPS_48hr_OD600)

nonpanortholog.data <- functional.data[functional.data$panortholog==FALSE,]

## panorthologs are more essential:  p-value = 2.441e-11
t.test(panortholog.data$EssentialityScore, nonpanortholog.data$EssentialityScore,alternative=c("greater"))
## Essentiality plot
quartz()
essentiality.plot <- ggplot(functional.data, aes(x=EssentialityScore, fill=panortholog)) + scale_fill_manual(values=c("purple", "green")) + geom_histogram(alpha=0.5, binwidth=0.5)
essentiality.plot

## panortholog KOs have worse growth in LB: p-value = 7.534e-05
t.test(panortholog.data$LB_22hr_OD600, nonpanortholog.data$LB_22hr_OD600,alternative=c("less"))
## LB histogram plot
LB.plot <- ggplot(functional.data, aes(x=LB_22hr_OD600, fill=panortholog)) + scale_fill_manual(values=c("purple", "green")) + geom_histogram(alpha=0.5, binwidth=0.01)
LB.plot

## panortholog KOs have worse growth in MOPS after 24 hours: p-value = 5.721e-09 
t.test(panortholog.data$MOPS_24hr_OD600, nonpanortholog.data$MOPS_24hr_OD600,alternative=c("less"))
## MOPS 24 hr histogram plot.
MOPS24hr.plot <- ggplot(functional.data, aes(x=MOPS_24hr_OD600, fill=panortholog)) + scale_fill_manual(values=c("purple", "green")) + geom_histogram(alpha=0.5, binwidth=0.01)
MOPS24hr.plot


## panortholog KOs have worse growth in MOPS after 48 hours: p-value = 3.943e-07
t.test(panortholog.data$MOPS_48hr_OD600, nonpanortholog.data$MOPS_48hr_OD600, alternative=c("less"))


## Ben Kerr's question: for genes with mutations in the LTEE,
#what the KEIO fitness effects, for genes with one hit versus 6 hits?
#Can I distinguish between effects where the gene is completely broken,
#versus the function of the gene being tuned?

## At the moment, I am not sure how to answer this question using the KEIO
## collection data.

#ltee.40K.mutations <- read.csv("../results/ltee_mutations_40K.csv")
# non.mutator.data <- subset(ltee.40K.mutations,ltee.40K.mutations$mutator==FALSE)


