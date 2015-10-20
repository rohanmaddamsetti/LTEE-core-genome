## Chattopadhyay-analysis.R by Rohan Maddamsetti

core.data <- read.csv("../results/core-summary.csv")

chattopadhyay <- read.table("../data/CFT073-NC_012967-RBB--evalue.rbb-panorthologs.txt")

chattopadhyay$Protein_id <- sapply(chattopadhyay$V2, function (x) substring(x, 23))

chattopadhyay$V2 <- NULL
chattopadhyay$CFT073_gene <- chattopadhyay$V1
chattopadhyay$V1 <- NULL
chattopadhyay$chattopadhyay_gene_hit <- rep(TRUE,length(chattopadhyay$Protein_id))

core.data2 <- merge(core.data,chattopadhyay, all=TRUE)
core.data2[c("chattopadhyay_gene_hit")][is.na(core.data2[c("chattopadhyay_gene_hit")])] <- FALSE

diversity.data <- read.csv("../results/ecoli-protein-diversity.csv")

data <- merge(core.data2,diversity.data,all=TRUE)
data <- data[data$panortholog==TRUE,]

## Do Mann Whitney U test to see if Chattopadhyay genes are more variable.
mutated.in.chattopadhyay <- data[data$chattopadhyay_gene_hit==TRUE,]$diversity
absent.in.chattopadhyay <- data[data$chattopadhyay_gene_hit==FALSE,]$diversity

## Chattopadhyay genes are not more or less variable than other panorthologs.
wilcox.test(mutated.in.chattopadhyay,absent.in.chattopadhyay)

library(ggplot2)

A.plot <- ggplot(data, aes(x=chattopadhyay_gene_hit,diversity)) + geom_boxplot() + scale_x_discrete("mutated in chattopadhyay") + theme_classic()

## now do Jaccard index test.

## first calculate Jaccard index between genes that mutated in LTEE and in Chattopadhyay paper.

ltee.hit.data <- read.csv("../results/protein-site-diversity.csv")
ltee.hit.data$mutated_site_diversity <- NULL
ltee.hit.data$rest_protein_diversity <- NULL

ltee.hit.genes <- ltee.hit.data$locus_tag
chattopadhyay.hit.genes <- data[data$chattopadhyay_gene_hit==TRUE,]$locus_tag

jaccard.index <- function(g1,g2) {length(intersect(g1,g2))/length(union(g1,g2))}

jaccard.value <- jaccard.index(ltee.hit.genes,chattopadhyay.hit.genes)
## value is 0.02. how extreme is this value when you calculate the jaccard index of LTEE hits to
## random sets of panorthologs?

N <- 50000 #number of random jaccard values.
g2.size <- length(chattopadhyay.hit.genes)

panortholog.locus.tags <- data[data$panortholog==TRUE,]$locus_tag

randos <- sapply(1:N,function (x) jaccard.index(ltee.hit.genes,sample(panortholog.locus.tags,g2.size)))

p.value <- sum(sapply(randos,function(x) ifelse(x<jaccard.value | x> (1-jaccard.value),1,0)))/N
## p.value is 0.25-- no significant convergence between LTEE and Chattopadhyay data sets.
