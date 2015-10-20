## site-diversity-analysis.R by Rohan Maddamsetti.

data <- read.csv("/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/results/protein-site-diversity.csv")

wilcox.test(data$mutated_site_diversity, data$rest_protein_diversity, paired=TRUE)

## p = 2.139 * 10^-5. Slow-evolving codons are under selection in the LTEE.
## protein sites under selection in the LTEE are clearly not fast-evolving in nature.

wilcox.test(data$mutated_site_diversity, data$rest_protein_diversity, paired=TRUE,alternative=c("less"))

## Only consider positions with any diversity.
variable.site.data <- data[data$mutated_site_diversity>0,]

wilcox.test(variable.site.data$mutated_site_diversity, variable.site.data$rest_protein_diversity, paired=TRUE, alternative=c("greater"))

sum(variable.site.data$LTEE_dN)
sum(data$LTEE_dN)
