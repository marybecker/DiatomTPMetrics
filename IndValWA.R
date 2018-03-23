#######Calculate Weight Averages,Tolerance Values and#################
######Indicator Value Species with occurence greater than 5###########
######replicating Potapova & Charles (2007) with CT specific data#####

setwd ("")##Set WD

library(indicspecies)
library(rioja)

#####Read in Data####################################################
spp <- read.csv ("SPP_RelAbund.csv",header=TRUE,row.names=1)
spp[is.na(spp)] <- 0
TP <- read.csv ("TP.csv",header=TRUE,row.names=1)
colnames(TP)[1]<- "TP"

taxa.names <- names(spp)

#####Calculate Weighted Average Optima and Tolerance (Stdev)#########

WA<- WA(spp,TP$TP,tolDW=TRUE)
WA<- as.data.frame(coef(WA))
WA$TORatio<- WA$Tolerances/WA$Optima
Q25<-quantile(WA$Optima,0.25)
Q75<-quantile(WA$Optima,0.75)
WA$WAGRP <- ifelse(WA$Optima <= Q25,"Decreasing",ifelse(WA$Optima >= Q75, "Increasing","I"))
WA$WAInclude<- ifelse(WA$WAGRP=="Decreasing" & WA$TORatio <3 
                      | WA$WAGRP=="Increasing" & WA$TORatio <3,1,0)

#####Indicator Species Analysis#####################################

TP$GRP <- ifelse(TP$TP <= 0.02,"L",ifelse(TP$TP > 0.065, "H","I"))

indval<- multipatt(spp,TP$GRP,control=how(nperm=999),duleg=TRUE)

indval.summary <- indval$sign
indval.summary$IndValGRP<- ifelse(indval.summary$index==1,"Increasing",
                                  ifelse(indval.summary$index=="2","I","Decreasing"))
indval.summary$Sig<- ifelse(indval.summary$p.value<0.05,1,0)

####Merge Results#################################################

IndValWASpp<- merge(indval.summary,WA,by=0)
IndValWASpp$tolcl<- ifelse(IndValWASpp$Sig==1,IndValWASpp$IndValGRP,
                            ifelse(IndValWASpp$WAInclude==1,IndValWASpp$WAGRP,"NA"))
IndValWASigSpp<- IndValWASpp[which(IndValWASpp$tolcl != "NA" & IndValWASpp$tolcl != "I" ),]
colnames(IndValWASigSpp)[1]<- "Taxa"

write.csv(IndValWASpp,"IndValWA_SppResults.csv")
write.csv(IndValWASigSpp,"IndValWA_SignificantSpp.csv",row.names=FALSE)




