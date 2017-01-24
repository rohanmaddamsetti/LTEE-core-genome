## KEIO-analysis.R by Rohan Maddamsetti

## Do a t-test on KEIO essentiality and growth defect scores comparing
## panorthologs and non-panorthologs.

## NEW: compare essentiality and G-score data.

library(ggplot2)
library(dplyr)
library(ggrepel)

## Analysis using G-score data (exactly same results as with 40K single clone data).

G.score.data <- read.csv("../data/Tenaillon-data/nature18959-s2.csv")
core.summary <- read.csv("../results/core-summary.csv")
KEIO.data <- read.csv("../data/KEIOdata.csv")

## turns out that start position is a good unique key for genes.
first.merge <- merge(G.score.data, core.summary, all=TRUE)

## for other analyses, write out the top 57 genes with 2 or more dN.
first.merge.57 <- filter(first.merge,Observed.nonsynonymous.mutation>1) %>%
    arrange(Gene.name)
write.csv(x=first.merge.57,file="../results/57-core-status.csv")


## drop all irrelevant columns. "note" is key for merge with KEIO data.
second.merge <- first.merge[c("Gene.name", "Gene.order", "Start.position", "G.score", "locus_tag", "panortholog", "note","Observed.nonsynonymous.mutation")]

final.data <- merge(second.merge, KEIO.data)

panortholog.data <- final.data[final.data$panortholog==TRUE,]
positive.panorthologs <- panortholog.data[panortholog.data$G.score>0,]
zero.panorthologs <- panortholog.data[panortholog.data$G.score==0,]

nonpanortholog.data <- final.data[final.data$panortholog==FALSE,]

top.57.panorthologs <- panortholog.data[panortholog.data$Observed.nonsynonymous.mutation>=2,]
single.hit.panorthologs <- panortholog.data[panortholog.data$Observed.nonsynonymous.mutation==1,]

wilcox.test(top.57.panorthologs$EssentialityScore,single.hit.panorthologs$EssentialityScore)
wilcox.test(top.57.panorthologs$EssentialityScore,zero.panorthologs$EssentialityScore)

panortholog.data2 <- panortholog.data
panortholog.data2$testcol <- as.factor(sapply(panortholog.data$Observed.nonsynonymous.mutation,function(x) if (x==0) {1} else if (x==1) {2} else if (x>1){3}))

## as a sanity check, plot Essentiality Score distributions.
quartz()
test.plot <- ggplot(panortholog.data2,aes(x=EssentialityScore,fill=testcol)) + geom_histogram()
test.plot

panortholog.set1 <- filter(panortholog.data2,testcol==1)
panortholog.set1.c <- filter(panortholog.data2,testcol!=1)

## p-value = 0.12. No difference in essentiality between panorthologs that mutated in the LTEE and those that did not.
t.test(panortholog.set1.c$EssentialityScore,panortholog.set1$EssentialityScore)
t.test(positive.panorthologs$EssentialityScore, zero.panorthologs$EssentialityScore)

## p-value = 0.43. No difference in LB growth between panorthologs that mutated in the LTEE and those that did not.
t.test(positive.panorthologs$LB_22hr_OD600, zero.panorthologs$LB_22hr_OD600)

## p-value = 0.72 No difference in MOPS 24 hour growth between panorthologs that mutated in the LTEE and those that did not.
t.test(positive.panorthologs$MOPS_24hr_OD600, zero.panorthologs$MOPS_24hr_OD600)

## p-value = 0.76. No difference in MOPS 48 hour growth between panorthologs that mutated in the LTEE and those that did not.
t.test(positive.panorthologs$MOPS_48hr_OD600, zero.panorthologs$MOPS_48hr_OD600)

## panorthologs are more essential:  p-value = 2.441e-11
t.test(panortholog.data$EssentialityScore, nonpanortholog.data$EssentialityScore,alternative=c("greater"))

## Essentiality plot
quartz()
essentiality.plot <- ggplot(final.data, aes(x=EssentialityScore, fill=panortholog)) + scale_fill_manual(values=c("purple", "green")) + geom_histogram(alpha=0.5, binwidth=0.5)
essentiality.plot

## panortholog KOs have worse growth in LB: p-value = 7.534e-05
t.test(panortholog.data$LB_22hr_OD600, nonpanortholog.data$LB_22hr_OD600,alternative=c("less"))
## LB histogram plot
LB.plot <- ggplot(final.data, aes(x=LB_22hr_OD600, fill=panortholog)) + scale_fill_manual(values=c("purple", "green")) + geom_histogram(alpha=0.5, binwidth=0.01)
LB.plot

## panortholog KOs have worse growth in MOPS after 24 hours: p-value = 5.721e-09
t.test(panortholog.data$MOPS_24hr_OD600, nonpanortholog.data$MOPS_24hr_OD600,alternative=c("less"))
## MOPS 24 hr histogram plot.
MOPS24hr.plot <- ggplot(final.data, aes(x=MOPS_24hr_OD600, fill=panortholog)) + scale_fill_manual(values=c("purple", "green")) + geom_histogram(alpha=0.5, binwidth=0.01)
MOPS24hr.plot


## panortholog KOs have worse growth in MOPS after 48 hours: p-value = 3.943e-07
t.test(panortholog.data$MOPS_48hr_OD600, nonpanortholog.data$MOPS_48hr_OD600, alternative=c("less"))

## compare essentiality against G-scores.
## Add an asterisk to the name of genes that that parallel mutations at the AA level:
asterisk.these <- factor(c('atoC','hflB','rpsD','infC','nadR','pykF','yijC','rplF','infB','spoT'))
final.data2 <- merge(first.merge, KEIO.data) %>%
    filter(Observed.nonsynonymous.mutation>1) %>%
    mutate(possible.KO=(IS.insertion+Short.indel+Large.deletion)>0) %>%
    arrange(EssentialityScore) %>% ## next line adds the asterisk.
    mutate(Gene.name = as.character(Gene.name)) %>%
    mutate(Gene.name=ifelse(Gene.name %in% asterisk.these,paste('***',Gene.name),Gene.name))

my.lab <- expression(paste('log(',italic('G'),' score + 1)'))

essentiality.G.score.plot <- ggplot(data=final.data2,aes(x=EssentialityScore,y=log10(G.score+1),shape=panortholog)) + geom_point(size=1) + scale_y_continuous(breaks=seq(0,2,by=0.5)) + scale_x_continuous(breaks=seq(-4,4,by=1)) + geom_text_repel(aes(label=Gene.name,colour=possible.KO,fontface='italic'), segment.size=0.1,size=2.5,angle=45) + scale_colour_manual(values=c("purple", "green4")) + theme_classic() + ylab(my.lab) + xlab("KEIO Essentiality Score") + theme(legend.position="none")

ggsave(filename="~/Desktop/Essentiality_G-score.pdf",plot=essentiality.G.score.plot, width=2*3.25,height=2*3.25,units='in')

## Do Fisher's exact test on # KO mutation containing genes among genes with positive Essentiality scores,
## versus those with negative Essentiality scores.
ES.positives <- filter(final.data2,EssentialityScore > 0)
ES.negatives <- filter(final.data2,EssentialityScore < 0)
ES.zeros <- filter(final.data2,EssentialityScore == 0)

## Do a binomial test on number of possible KO mutation genes compared to all genes in top 17 versus bottom.
filter(ES.positives,possible.KO==TRUE) ## 2 positive genes with possible KO mutations out of 16 positives.
filter(ES.negatives,possible.KO==TRUE) ## 17 negative genes with possible KO mutations out of 40 negatives.
## 2  | 14
## 17 | 23
fisher.table <- matrix(c(2, 14, 17, 23),nrow = 2)
fisher.test(fisher.table,alternative="less")
