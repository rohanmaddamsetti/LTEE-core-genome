## plot-KO-parallelism.R by Rohan Maddamsetti.

## make a scatterplot of:
##  number of KO mutations for each gene across mutators and non-mutators.
##  number of lineages with mutation in this gene.

library(ggplot2)
library(dplyr)
library(stringr)

## count number of possible KO mutations across mutators and non-mutators.
tenaillon.non.mutator.data <- read.csv("../data/Tenaillon-data/nature18959-s2.csv")
tenaillon.mutator.data <- read.csv("../data/Tenaillon-data/nature18959-s3.csv")

KO.count.non.mutators <- tenaillon.non.mutator.data %>% mutate(non.mutator.possible.KOs=IS.insertion+Short.indel+Large.deletion)
KO.count.mutators <- tenaillon.mutator.data %>% mutate(mutator.possible.KOs=IS.insertion+Short.indel+Large.deletion)
all.KO.counts <- merge(KO.count.non.mutators,KO.count.mutators,all=TRUE)
all.KO.counts[is.na(all.KO.counts)] <- 0
all.KO.counts <- all.KO.counts %>% mutate(all.possible.KOs=non.mutator.possible.KOs+mutator.possible.KOs)

## count number of lineages with mutations in gene X (excluding intergenic and synonymous mutations).
processed.shiny.data <- read.csv("../results/processed_shiny_data.csv")
processed.shiny.data2 <- processed.shiny.data %>% mutate(lineage.count=str_count(lineages,":")+1) %>% rename(Gene.name = gene_name)

final.data <- merge(all.KO.counts,processed.shiny.data2)

quartz()
my.plot <- ggplot(final.data,aes(x=all.possible.KOs,y=lineage.count)) + geom_point() + theme_classic()

all.KO.counts %>% filter(all.possible.KOs>=20)
