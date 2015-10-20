##transpositionGrapher.R by Rohan Maddamsetti

##this quick script makes a stacked bar graph using the data in LTEE_40K_IS_events.csv.

library(ggplot2)
library(gridExtra)
data <- read.csv("LTEE_40K_IS_events.csv")
attach(data)
quartz()
plot1 <- ggplot(data,aes(x=Strain)) + geom_bar(aes(fill = Event))
plot2 <- ggplot(data,aes(x=Strain)) + geom_bar(aes(fill = Event), position = 'fill')
##now put the two plots side by side.
grid.arrange(plot1,plot2,ncol=1)
