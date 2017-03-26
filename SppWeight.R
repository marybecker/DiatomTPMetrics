setwd ("/Users/tbecker/Documents/Projects/GitHubProjects/DiatomTPMetrics") 

library(ggplot2)

spp <- read.csv ("data/SPP.csv",header=TRUE,row.names=1)

taxa.names <- names(spp)

sppSCnt<- as.data.frame(colSums(spp))
colnames(sppSCnt)[1]<-"Cnt"

ggplot(sppSCnt,aes(Cnt))+
  geom_histogram(binwidth=10)+
  geom_vline(xintercept=c(31.5,69))

quantile(sppSCnt$Cnt,0.269)

