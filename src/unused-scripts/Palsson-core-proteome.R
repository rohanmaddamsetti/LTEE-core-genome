## Palsson-core-proteome.R by Rohan Maddamsetti

library(ggplot2)
library(dplyr)

core.summary <- read.csv("../results/core-summary.csv") %>% rename(gene=note)

core.proteome <- read.csv("../Data/Palsson-data/pnas.1501384112.sd01.csv") %>% rename(Gene.name=name) %>% select(Gene.name,gene)

core.proteome2 <- merge(core.summary,core.proteome)

tenaillon.non.mutator.data <- read.csv("../data/Tenaillon-data/nature18959-s2.csv")
tenaillon.mutator.data <- read.csv("../data/Tenaillon-data/nature18959-s3.csv")

proteome.mutator.data <- merge(core.proteome,tenaillon.mutator.data) %>% select(Gene.name,Gene.order,Observed.nonsynonymous.mutation,Expected.nonsynonymous.mutation,G.score)
proteome.non.mutator.data <- merge(core.proteome,tenaillon.non.mutator.data) %>% select(Gene.name,Gene.order,Observed.nonsynonymous.mutation,Expected.nonsynonymous.mutation,G.score)

proteome.mutator.data2 <- proteome.mutator.data %>% filter(G.score !=0)
proteome.non.mutator.data2 <- proteome.non.mutator.data %>% filter(G.score !=0)

