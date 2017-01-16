## proteome-sector-analysis.R by Rohan Maddamsetti

## I found nothing; Not enough statistical power.

library(ggplot2)
library(dplyr)

core.summary <- read.csv("../results/core-summary.csv")

## There are 1053 genes in the Hui dataset.
proteome.sector.data <- read.csv("../data/Hui-data/proteome-sector-assignments.csv")

## Non-mutator data.
nonmut.data <- read.csv("../data/Tenaillon-data/nature18959-s2.csv")
nonmut.data2 <- merge(nonmut.G.score.data,core.summary,all=TRUE)
nonmut.data3 <- merge(nonmut.data2,proteome.sector.data,all=TRUE)

## Mutator data.
mut.data <- read.csv("../data/Tenaillon-data/nature18959-s3.csv")
mut.data2 <- merge(mut.G.score.data,core.summary,all=TRUE)
mut.data3 <- merge(mut.data2,proteome.sector.data,all=TRUE)

## how do core gene assignments relate to proteome sectors?
## are core genes over or under-represented in the Hui dataset?

nonmut.proteome.outsiders <- nonmut.data3 %>% filter(is.na(Sector.assigned))
mut.proteome.outsiders <- mut.data3 %>% filter(is.na(Sector.assigned))

nonmut.proteome.insiders <- nonmut.data3 %>% filter(!is.na(Sector.assigned))
mut.proteome.insiders <- mut.data3 %>% filter(!is.na(Sector.assigned))



## these should be the same
count(nonmut.proteome.outsiders)
count(mut.proteome.outsiders)

count(nonmut.proteome.insiders)
count(mut.proteome.insiders)


nonmut.proteome.outsiders %>% filter(panortholog==TRUE) %>% count() ## 1342 genes
nonmut.proteome.outsiders %>% filter(panortholog==FALSE) %>% count() ## 1939 genes

nonmut.proteome.insiders %>% filter(panortholog==TRUE) %>% count() ## 626 genes
nonmut.proteome.insiders %>% filter(panortholog==FALSE) %>% count() ## 292 genes

## total is 4199 genes (total should be related to difference b/t K-12 and REL606.

## Do LTEE mutations tend to fall in or out of the Hui dataset?
nonmut.proteome.insiders %>% filter(G.score > 0) %>% count() ## 87 genes
nonmut.proteome.outsiders %>% filter(G.score > 0) %>% count() ## 198 genes

mut.proteome.insiders %>% filter(G.score > 0) %>% count() ## 327 genes
mut.proteome.outsiders %>% filter(G.score > 0) %>% count() ## 1507 genes

## count the number of genes assigned to each proteome sector.
e1 <- nonmut.proteome.insiders %>% group_by(Sector.assigned) %>% summarise(total=n())
x1 <- nonmut.proteome.insiders %>% filter(G.score > 0) %>% group_by(Sector.assigned) %>% summarize(total=n())

e2 <- mut.proteome.insiders %>% group_by(Sector.assigned) %>% summarise(total=n())
x2 <- mut.proteome.insiders %>% filter(G.score > 0) %>% group_by(Sector.assigned) %>% summarize(total=n())

chisq.test(x1$total,e1$total,correct=FALSE)
chisq.test(x2$total,e2$total,correct=FALSE)
