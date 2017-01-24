## make-tree-orthology-summary.R by Rohan Maddamsetti

library(dplyr)

salmonella.summary <- read.csv("/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/results/salmonella-orthology-summary.csv")

ecoli.core <- read.csv("/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/results/core-summary.csv") %>% filter(panortholog == TRUE)

tree.genes <- inner_join(salmonella.summary,ecoli.core)
write.csv(tree.genes,file="/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/results/tree/tree-orthology-summary.csv",quote=FALSE,row.names=FALSE)
