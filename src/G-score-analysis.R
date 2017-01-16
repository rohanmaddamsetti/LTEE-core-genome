## G-score-analysis.R by Rohan Maddamsetti.

## This script makes scatterplots of the G-scores calculated in Tenaillon et al. 2016
## to protein diversity in panorthologs among E. coli and between REL606 and Salmonella.
## Protein diversity is the mean pairwise difference per site in proteins.

library(ggplot2)
library(dplyr)
library(asbio) ## This is to calculate confidence interval on median
## for figure part C.

## This function does the analysis and plots figures.
## The last line of the script calls this function.

G.score.analysis <- function(G.score.data,mutator=FALSE,diverge=TRUE) {
    ## EDIT THIS TO MAKE DIFFERENT FIG FOR MUTATORS.
    if (diverge == TRUE) {
        fname <- "Fig-2"
        differences <- "../results/ecoli-salmonella-protein-diversity.csv"
        print("ANALYSIS FOR E. COLI -- SALMONELLA DIVERGENCE")
    } else {
        fname <- "Fig-1"
        differences <- "../results/ecoli-protein-diversity.csv"
        print("ANALYSIS FOR E. COLI DIVERSITY")
    }
    if (mutator == TRUE) {
        fname <- paste("supp",fname,sep='-')
        print("HYPER-MUTATOR DATA ANALYSIS")
    } else {
        print("NON-MUTATOR DATA ANALYSIS")
    }


    diversity.data <- read.csv(differences)

    core.summary <- read.csv("../results/core-summary.csv")

    ## Turns out that start position is a good unique key for genes.
    first.merge <- merge(G.score.data, core.summary, all=TRUE)
    ## drop all irrelevant columns.
    first.merge <- first.merge[c("Gene.name", "Gene.order", "Start.position", "G.score", "locus_tag", "panortholog")]

    final.data <- merge(first.merge,diversity.data,all=TRUE)

    ## we want to only merge genes found in diversity.data (i.e. panorthologs).
    ## so since we have merged all genes,
    ## remove locus_tags that are not in diversity.data.
    final.data <- final.data[final.data$locus_tag %in% diversity.data$locus_tag,]

    my.lab <- expression(paste('log(',italic('G'),' score + 1)'))

###### Figure part A.
    A.data <- final.data
    ## log-transform G-score. For mutator data, deal with negative G-scores symmetrically.
    if (mutator == FALSE) {
        A.data$logG <- sapply(A.data$G.score, function(x) log10(x+1))
    } else {
        A.data$logG <- sapply(A.data$G.score, function(x) ifelse(x>0,log10(x+1),-log10(-x+1)))
    }

    ## sqrt-transform diversity.
    A.data$sqrt.diversity <- sapply(A.data$diversity, function(x) sqrt(x))
    A.plot <- ggplot(A.data,aes(logG,sqrt.diversity)) + geom_point(size=0.3) + scale_x_continuous(my.lab) + theme_classic(base_size=9)
    if (diverge == TRUE){
        A.plot <- A.plot + scale_y_continuous("sqrt(Divergence)")
    } else {
        A.plot <- A.plot + scale_y_continuous("sqrt(Diversity)")
    }
    quartz()
    A.plot + coord_flip()
    ggsave(paste0("/Users/Rohandinho/Desktop/", fname,"-A" ,".pdf"), scale=0.32)
    dev.off()
    print(cor.test(x=A.data$G.score, y=A.data$diversity,method="spearman",exact=FALSE))

############################# Figure part B. plot only genes that mutated.
    B.data <- final.data
    B.data <- B.data[B.data$G.score!=0,]
    ## log-transform G-score. For mutator data, deal with negative G-scores symmetrically.
    if (mutator == FALSE) {
        B.data$logG <- sapply(B.data$G.score, function(x) log10(x+1))
    } else {
        B.data$logG <- sapply(B.data$G.score, function(x) ifelse(x>0,log10(x+1),-log10(-x+1)))
    }
    ## sqrt-transform diversity.
    B.data$sqrt.diversity <- sapply(B.data$diversity, function(x) sqrt(x))
    B.plot <- ggplot(B.data,aes(logG, sqrt.diversity)) + geom_point(size=0.3) + scale_x_continuous(my.lab) + theme_classic(base_size=9)
    if (diverge == TRUE){
        B.plot <- B.plot + scale_y_continuous("sqrt(Divergence)")
    } else {
        B.plot <- B.plot + scale_y_continuous("sqrt(Diversity)")
    }
    quartz()
    B.plot + coord_flip()
    ggsave(paste0("/Users/Rohandinho/Desktop/", fname,"-B",".pdf"), scale=0.32)
    dev.off()
    print(cor.test(x=B.data$G.score, y=B.data$diversity,method="spearman",exact=FALSE))

######### FOR NON-MUTATOR DATA:
######### Figure part C. Mann-Whitney test between G-score = 0 and G-score > 0 for non-mutators.
    if (mutator == FALSE) {
    C.data <- final.data
    C.data$mutated <- as.factor(sapply(C.data$G.score, function(x) ifelse(x>0,"Positive", "Zero")))

    ## sqrt-transform diversity.
    C.data$sqrt.diversity <- sapply(C.data$diversity, function(x) sqrt(x))

    ## To calculate 95% CI for the median of the distributions and for M-W U-test.
    ## These use sqrt-transformed data!!!
    G.mutated.diversity <- C.data[C.data$mutated=='Positive',]$sqrt.diversity
    G.non.mutated.diversity <- C.data[C.data$mutated=='Zero',]$sqrt.diversity

    ## Do Mann Whitney U test to see if mutated genes are more conserved.
    print("Mann-Whitney U-test for Positive G-scores against zero G-scores")
    print(wilcox.test(G.mutated.diversity,G.non.mutated.diversity,correct=FALSE))

    ## make a data frame to plot confidence interval of median for the two dists.
    mutated.conf.int <- ci.median(G.mutated.diversity)$ci
    non.mutated.conf.int <- ci.median(G.non.mutated.diversity)$ci

    summary.C.data <- data.frame(category=c('Zero', 'Positive'),
                                 median=c(non.mutated.conf.int[1],
                                          mutated.conf.int[1]),
                                 left.CI=c(non.mutated.conf.int[2],
                                           mutated.conf.int[2]),
                                 right.CI=c(non.mutated.conf.int[3],
                                            mutated.conf.int[3]))

    ## reverse levels 'Non-zero', 'Zero' for plotting.
    summary.C.data$category <- factor(summary.C.data$category,levels(summary.C.data$category)[c(2,1)])

    ## add numerical labels to the lineage categories.
    ## first count number of points in each category.

    x1 <- levels(C.data$mutated)
    x2 <- sapply(x1, function (y) sum(sapply(C.data$mutated, function (x) ifelse(x==y,1,0))))
    ## get the right.CI value in the summary to put the label position
    x3 <- sapply(x1, function (y) summary.C.data[summary.C.data$category==y,]$median)
    labelz <- data.frame(mutated=x1, counts=x2, ordinate=x3)

g.score.lab <- expression(paste(italic('G'),' score'))

    C.plot <- ggplot(data=summary.C.data, aes(x=category,y=median)) + geom_errorbar(aes(ymin=left.CI,ymax=right.CI), width=0.25) + geom_point() + scale_x_discrete(g.score.lab) + geom_text(data = labelz, aes(x=factor(mutated), y=ordinate, label=paste("N =", counts)), vjust=-3, size=2.6) + theme_classic(base_size=9)


    ## the 'limits' option is to keep the number of genes label on the figure.
    if (diverge == TRUE){
        C.plot <- C.plot + scale_y_continuous("sqrt(Divergence)", breaks=scales::pretty_breaks(n=3))
    } else {
        C.plot <- C.plot + scale_y_continuous("sqrt(Diversity)", breaks=scales::pretty_breaks(n=3))
    }

    quartz()
    C.plot + coord_flip()
    ggsave(paste0("/Users/Rohandinho/Desktop/", fname,"-C" ,".pdf"), scale=0.32)
    dev.off()
    ## redo the test by resampling 163 random genes from the data table, and asking
    ## how often the random sample contains less diverged genes than those with G-scores.
    ##  observed.mean.conservation <- mean(G.mutated.diversity)
    ##  test.p.val <- sum(sapply(replicate(10000,mean(sample_n(final.data,163)$diversity)), function(x) ifelse(x<observed.mean.conservation,1,0)))/10000
    ## both tests are significant,
    ## second test is significant: p < 0.015.
    ##  print(test.p.val)
    }
    ###############################################################
    ## FOR MUTATOR DATA:
    if (mutator == TRUE) {
        C.data <- final.data
        C.data$mutated <- as.factor(sapply(C.data$G.score, function(x) if(x>0) {"Positive"} else if (x<0) {"Negative"} else {"Zero"}))

        ## sqrt-transform diversity.
        C.data$sqrt.diversity <- sapply(C.data$diversity, function(x) sqrt(x))

        ## To calculate 95% CI for the median of the distributions and for M-W U-test.
        ## These use sqrt-transformed data!!!
        G.positive.diversity <- C.data[C.data$mutated=='Positive',]$sqrt.diversity
        G.zeros.diversity <- C.data[C.data$mutated=='Zero',]$sqrt.diversity
        G.negative.diversity <- C.data[C.data$mutated=='Negative',]$sqrt.diversity
        G.mutated.diversity <- c(G.positive.diversity,G.negative.diversity)

        ## Do Mann Whitney U test to see if mutated genes are more conserved.
        print("Mann-Whitney U-test for Positive G-scores against Zero G-scores")
        print(wilcox.test(G.positive.diversity,G.zeros.diversity,correct=FALSE))

        print("Mann-Whitney U-test for Negative G-scores against Zero G-scores")
        print(wilcox.test(G.negative.diversity,G.zeros.diversity,correct=FALSE))


        ## make a data frame to plot confidence interval of median for the three dists.
        positive.conf.int <- ci.median(G.positive.diversity)$ci
        zeros.conf.int <- ci.median(G.zeros.diversity)$ci
        negative.conf.int <- ci.median(G.negative.diversity)$ci

    summary.C.data <- data.frame(category=c('Negative','Zero', 'Positive'),
                                 median=c(negative.conf.int[1],
                                          zeros.conf.int[1],
                                          positive.conf.int[1]),
                                 left.CI=c(negative.conf.int[2],
                                           zeros.conf.int[2],
                                           positive.conf.int[2]),
                                 right.CI=c(negative.conf.int[3],
                                            zeros.conf.int[3],
                                            positive.conf.int[3]))

    ## reverse levels 'Positive', 'Zero', 'Negative' for plotting. Note that levels are ordered alphabetically.
    summary.C.data$category <- factor(summary.C.data$category,levels(summary.C.data$category)[c(1,3,2)])

    ## add numerical labels to the lineage categories.
    ## first count number of points in each category.

    x1 <- levels(C.data$mutated)
    x2 <- sapply(x1, function (y) sum(sapply(C.data$mutated, function (x) ifelse(x==y,1,0))))
    ## get the right.CI value in the summary to put the label position
    x3 <- sapply(x1, function (y) summary.C.data[summary.C.data$category==y,]$median)
    labelz <- data.frame(mutated=x1, counts=x2, ordinate=x3)

    C.plot <- ggplot(data=summary.C.data, aes(x=category,y=median)) + geom_errorbar(aes(ymin=left.CI,ymax=right.CI), width=0.25) + geom_point() + scale_x_discrete(g.score.lab) + geom_text(data = labelz, aes(x=factor(mutated), y=ordinate, label=paste("N =", counts)), vjust=-3, size=2.6) + theme_classic(base_size=9)


    ## the 'limits' option is to keep the number of genes label on the figure.
    if (diverge == TRUE){
        C.plot <- C.plot + scale_y_continuous("sqrt(Divergence)", breaks=scales::pretty_breaks(n=3))
    } else {
        C.plot <- C.plot + scale_y_continuous("sqrt(Diversity)", breaks=scales::pretty_breaks(n=3))
    }

    quartz()
    C.plot + coord_flip()
    ggsave(paste0("/Users/Rohandinho/Desktop/", fname,"-C" ,".pdf"), scale=0.32)
    dev.off()
    }

    return()
}


non.mut.G.score.data <- read.csv("../data/Tenaillon-data/nature18959-s2.csv")

## analyze E. coli polymorphism.
G.score.analysis(non.mut.G.score.data,mutator=FALSE,diverge=FALSE)

## analyze E. coli-Salmonella divergence.
G.score.analysis(non.mut.G.score.data,mutator=FALSE,diverge=TRUE)

## The result is FAR weaker for the polymorphism data.
## My interpretation is that since these data are already subset
## on the core, there is less signal that more conserved/less polymorphic
##core genes are under positive selection in the LTEE.
