## assess-runs.R by Rohan Maddamsetti.

## look at the output to judge burn-in. 

## NOTE!: need replicate runs to check for convergence! I haven't done this.

## import omegaMap R functions.
source("/Users/Rohandinho/Desktop/Projects/LTEE-evolutionary-rates/src/omegamap/Documentation/R-functions.txt")

error <- stop ## Because the function error should be the function stop in Daniel
              ## Wilson's R code. 

run1 <- open.omegaMap("/Users/Rohandinho/Desktop/Projects/LTEE-evolutionary-rates/results/omegaMap-output/polymorphism/ECB_03156.txt")

run2 <- open.omegaMap("/Users/Rohandinho/Desktop/Projects/LTEE-evolutionary-rates/results/omegaMap-output/polymorphism/ECB_02110.txt")

run3 <- open.omegaMap("/Users/Rohandinho/Desktop/Projects/LTEE-evolutionary-rates/results/omegaMap-output/polymorphism/ECB_00110.txt")

runs <- list(run3=run3)

## assess burn-in. 15000 seems like an adequate burn-in period.
quartz()
trace.omegaMap(runs,"mu")
trace.omegaMap(runs,"mu",xlim=c(0,15000))
trace.omegaMap(runs,"mu",xlim=c(15000,150000)) 


### Assess convergence of omega. This doesn't make sense unless examining replicate
### runs.
mycols <- c("red","green") # Change the default colours 
plot.omega.converged(runs,cols=mycols) 


