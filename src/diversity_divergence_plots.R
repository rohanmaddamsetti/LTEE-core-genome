## diversity_divergence_plots.R by Rohan Maddamsetti

## This script makes boxplots of the mean pairwise difference per site in amino acid
## sequence for each gene shared between REL606 and Salmonella typhimurium, and
## also for those shared among 60 E. coli genomes.
## I calculate a correlation between protein diversity
## and the number of independent hits per nonmutator lineage in the LTEE
## at 40K generations.

library(ggplot2)

## This function does the analysis; the last line of the script
## calls this function.

plot.diversity <- function(diverge=TRUE) {
  if (diverge == TRUE) {
    fname <- "ecoli-salmonella-protein-diversity"
    differences <- "../results/ecoli-salmonella-protein-diversity.csv"
  } else {
    fname <- "ecoli-protein-diversity"
    differences <- "../results/ecoli-protein-diversity.csv"
  }
    
  diversity.data <- read.csv(differences)

    ## count the number of lineages with at least one dN per locus.
  ltee.40K.mutations <- read.csv("../results/ltee_mutations_40K.csv")
  non.mutator.data <- subset(ltee.40K.mutations,ltee.40K.mutations$mutator==FALSE)
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
  
  final.data <- merge(final.hits, diversity.data, all=TRUE)
  final.data$lineages <- sapply(final.data$lineages, function(x) ifelse(is.na(x),0,x))

### some data is missing (should be mutations outside of panorthologs):
  problems <- final.data[is.na(final.data$diversity),]

###### Figure part A. Mann-Whitney test between 0s and 1-6s.
  A.data <- final.data
  A.data$mutated <- sapply(A.data$lineages, function (x) ifelse(x>0,"1-6","0"))

## add numerical labels to the lineage categories.
  ## first count number of points in each category.
  x1 <- levels(as.factor(A.data$mutated))
  x2 <- sapply(x1, function (y) sum(sapply(A.data$mutated, function (x) ifelse(x==y,1,0))))
  ## get the max diversity value in each category to put the label position
  x3 <- sapply(x1, function (y) max(A.data[A.data$mutated==y,]$diversity,na.rm=TRUE))
  labelz <- data.frame(mutated=x1, counts=x2, ordinate=x3)

  y.max <- max(A.data$diversity, na.rm=TRUE)
  y.max <- y.max + 0.1*y.max  ## so the numerical labels don't fall off the panel.

  A.plot <- ggplot(A.data, aes(x=mutated,diversity)) + geom_boxplot() + scale_x_discrete("Number of lineages with non-synonymous mutations") + geom_text(data = labelz, aes(x=factor(mutated), y=ordinate, label=paste("N =", counts)), hjust=-0.3, size=3) + ylim(0,y.max) + theme_classic()

  A.plot + coord_flip()

  ggsave(paste0("/Users/Rohandinho/Desktop/", fname,"-A" ,".pdf"))

  ## Do Mann Whitney U test to see if mutated genes are more conserved.
  A.mutated.data <- A.data[A.data$mutated=="1-6",]$diversity
  A.non.mutated.data <- A.data[A.data$mutated=="0",]$diversity
  wilcox.test(A.mutated.data,A.non.mutated.data)

############################# Figure part B. plot only genes that mutated.
  B.data <- final.data
  B.data <- B.data[B.data$lineages>0,]

## add numerical labels to the lineage categories.
  x1 <- as.integer(levels(as.factor(B.data$lineages)))
  ## count number of points in each category.
  x2 <- sapply(x1, function (y) sum(sapply(B.data$lineages, function (x) ifelse(x==y,1,0))))
  ## get the max diversity value in each category to put the label position
  x3 <- sapply(x1, function (y) max(B.data[B.data$lineages==y,]$diversity,na.rm=TRUE))
  labelz <- data.frame(lineages=x1, counts=x2, ordinate=x3)

  y.max <- max(final.data$diversity, na.rm=TRUE)
  y.max <- y.max + 0.1*y.max  ## so the numerical labels don't fall off the panel.
  
  ## graph with box plots.
  
  B.plot <- ggplot(B.data, aes(x=factor(lineages), diversity)) + geom_boxplot() + scale_x_discrete("Number of lineages with non-synonymous mutations") + geom_text(data = labelz, aes(x=factor(lineages), y=ordinate, label=paste("N =", counts)), hjust=-0.3, size=3) + ylim(0,y.max) + theme_classic()

  B.plot + coord_flip()
  ggsave(paste0("/Users/Rohandinho/Desktop/", fname,"-B",".pdf"))
  cor.test(x=B.data$lineages, y=B.data$diversity,method="spearman",exact=FALSE)

  
############################## Figure part C. make diversity plot.
  C.data <- final.data
## add numerical labels to the lineage categories.
  x1 <- as.integer(levels(as.factor(C.data$lineages)))
  ## count number of points in each category.
  x2 <- sapply(x1, function (y) sum(sapply(C.data$lineages, function (x) ifelse(x==y,1,0))))
  ## get the max diversity value in each category to put the label position
  x3 <- sapply(x1, function (y) max(C.data[C.data$lineages==y,]$diversity,na.rm=TRUE))
  labelz <- data.frame(lineages=x1, counts=x2, ordinate=x3)

  y.max <- max(final.data$diversity, na.rm=TRUE)
  y.max <- y.max + 0.1*y.max  ## so the numerical labels don't fall off the panel.
  
  ## graph with box plots.
  
  C.plot <- ggplot(C.data, aes(x=factor(lineages), diversity)) + geom_boxplot() + scale_x_discrete("Number of lineages with non-synonymous mutations") + geom_text(data = labelz, aes(x=factor(lineages), y=ordinate, label=paste("N =", counts)), hjust=-0.3, size=3) + ylim(0,y.max) + theme_classic()
  
  C.plot + coord_flip()

  ggsave(paste0("/Users/Rohandinho/Desktop/", fname,"-C",".pdf"))

  cor.test(x=C.data$lineages, y=C.data$diversity,method="spearman",exact=FALSE)
}

##################################
   

## analyze E. coli-Salmonella divergence.
plot.diversity(diverge=TRUE)
## analyze E. coli polymorphism.
plot.diversity(diverge=FALSE)

## The result is FAR weaker for the polymorphism data.
## My interpretation is that since these data are already subset
## on the core, there is less signal that more conserved/less polymorphic
##core genes are under positive selection in the LTEE.
