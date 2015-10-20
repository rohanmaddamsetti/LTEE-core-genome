## divergence_scatterplot.R by Rohan Maddamsetti

## This script makes a simple scatterplot of E. coli-Salmonella DNA divergence for
## each gene, and the number of independent hits per nonmutator lineage in the LTEE
##at 50K generations.

library(ggplot2)

## This function does the analysis; the last line of the script
## calls this function.

plot.differences <- function(diverge=TRUE, dna=TRUE) {
  if (diverge == TRUE) {
    differences <- "../results/ecoli-salmonella-divergence.csv"
    if (dna == TRUE)
      the.title <- "Sequence divergence"
    else
      the.title <- "Protein divergence"
  } else {
    differences <- "../results/ecoli-polymorphism.csv"
    if (dna == TRUE)
      the.title <- "Sequence polymorphism"
    else
      the.title <- "Protein polymorphism"
  }
  
  divergence.data <- read.csv(differences)
  divergence.data$divergence <- divergence.data$differences/divergence.data$sites
  
  ltee.50K.mutations <- read.csv("../results/ltee_mutations_50K.csv")
  non.mutator.data <- subset(ltee.50K.mutations,ltee.50K.mutations$mutator==FALSE)
  
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
    
  final.data <- merge(final.hits, divergence.data, all=TRUE)
  final.data$lineages <- sapply(final.data$lineages, function(x) ifelse(is.na(x),0,x))
  final.data <- final.data[!is.na(final.data$divergence),]

## make divergence plot.
  
  quartz()
  divergence.plot <- ggplot(final.data, aes(x=lineages, y=divergence)) + geom_point() +
    theme_bw() +
                                        #eliminates background, gridlines, and chart border
      theme(
        plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        )
  
  divergence.plot + scale_x_continuous("Number of lineages containing a novel allele") + scale_y_continuous(the.title) + ggtitle("Conserved genes evolve in parallel")


 ## calculate rank-correlation on divergence/polymorphism scatterplot,
## and do permutation test for significance.

## Results of test on E.coli-Salmonella divergence:
##	Spearman's rank correlation rho

##data:  final.data$lineages and final.data$divergence
##S = 4708178572, p-value = 6.417e-10
##alternative hypothesis: true rho is less than 0
##sample estimates:
##       rho 
##-0.1116351 

## Results of test on E.coli polymorphism: 
##  	Spearman's rank correlation rho

##data:  final.data$lineages and final.data$divergence
##S = 1341251545, p-value = 0.005429
##alternative hypothesis: true rho is less than 0
##sample estimates:
##        rho 
##-0.05742243 

  cor.test(x=final.data$lineages, y=final.data$divergence,method="spearman",exact=FALSE, alternative="less")
}

## analyze E. coli-Salmonella divergence.
plot.differences(diverge=TRUE, dna=TRUE)
plot.differences(diverge=TRUE, dna=FALSE)
## analyze E. coli polymorphism.
plot.differences(diverge=FALSE, dna=TRUE)
plot.differences(diverge=FALSE, dna=FALSE)

## The result is FAR weaker for the polymorphism data.
## My interpretation is that since these data are already subset
## on the core, there is less signal that more conserved/less polymorphic
##core genes are under positive selection in the LTEE.
