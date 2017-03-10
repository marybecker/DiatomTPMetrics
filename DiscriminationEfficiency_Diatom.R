setwd ()#SETWD

INDICES <-  read.csv ("data/DiatomMetrics_112316.csv",sep=",",header=TRUE)
pHINDICES<- read.csv ("data/DiatomMetricspH.csv",sep=",",header=TRUE)
TINDICES<-  read.csv ("data/DiatomMetricsJulyTemp.csv",sep=",",header=TRUE)


library(dunn.test)
library(ggplot2)
library(grid)
library(plyr)
library(reshape2)

#####MULTI-Plot Function##############
#######################################

#call this with p1,p,2,... cols=4
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#Calculate the discrimination efficiency of metrics#

######Calculates that Discrimination Efficiency for Tolerant Metric for different combinations of TP Groups#
#############################################
 LGrp <- INDICES[which(INDICES$GRP=='L'),]
 L75 <- quantile(LGrp$H,0.75)
 MGrp <- INDICES[which(INDICES$GRP=='M'),]
 M75 <- quantile(MGrp$H,0.75)

#LH DE#
GrpTotal <- INDICES[which(INDICES$GRP=='H'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$GRP=='H'& INDICES$H >L75),]
n_subset <- length(GrpSubset$ID)
n_subset/n

#LM DE#
GrpTotal <- INDICES[which(INDICES$GRP=='M'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$GRP=='M'& INDICES$H >L75),]
n_subset <- length(GrpSubset$ID)
n_subset/n

#MH#
GrpTotal <- INDICES[which(INDICES$GRP=='H'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$GRP=='H'& INDICES$H >M75),]
n_subset <- length(GrpSubset$ID)
n_subset/n

#####Calculates that Discrimination Efficiency for Sensitive metric for different combinations of TP Groups##
################################################
 LGrp <- INDICES[which(INDICES$GRP=='L'),]
 L25 <- quantile(LGrp$L,0.25)
 MGrp <- INDICES[which(INDICES$GRP=='M'),]
 M25 <- quantile(MGrp$L,0.25)

#LH DE#
GrpTotal <- INDICES[which(INDICES$GRP=='H'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$GRP=='H'& INDICES$L <L25),]
n_subset <- length(GrpSubset$ID)
n_subset/n

#LM DE#
GrpTotal <- INDICES[which(INDICES$GRP=='M'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$GRP=='M'& INDICES$L <L25),]
n_subset <- length(GrpSubset$ID)
n_subset/n

#MH#
GrpTotal <- INDICES[which(INDICES$GRP=='H'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$GRP=='H'& INDICES$L <M25),]
n_subset <- length(GrpSubset$ID)
n_subset/n

#####Calculates that Discrimination Efficiency for T to S Ratio metric for different combinations of TP Groups##
################################################
LGrp <- INDICES[which(INDICES$GRP=='L'),]
L75 <- quantile(LGrp$R,0.75)
MGrp <- INDICES[which(INDICES$GRP=='M'),]
M75 <- quantile(MGrp$R,0.75)

#LH DE#
GrpTotal <- INDICES[which(INDICES$GRP=='H'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$GRP=='H'& INDICES$R >L75),]
n_subset <- length(GrpSubset$ID)
n_subset/n

#LM DE#
GrpTotal <- INDICES[which(INDICES$GRP=='M'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$GRP=='M'& INDICES$R >L75),]
n_subset <- length(GrpSubset$ID)
n_subset/n

#MH#
GrpTotal <- INDICES[which(INDICES$GRP=='H'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$GRP=='H'& INDICES$R >M75),]
n_subset <- length(GrpSubset$ID)
n_subset/n


#####Boxplots##
################################################
INDICES$GRP<- factor(INDICES$GRP,levels=c("L","M","H"))
pHINDICES$GRP<- factor(pHINDICES$GRP,levels=c("L","M","H"))
TINDICES$GRP<- factor(TINDICES$GRP,levels=c("L","M","H"))

bp1<- ggplot(INDICES,aes(x=GRP,y=H,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  labs(y="Relative Abundance TP Tolerant")+
  theme(legend.position="none",axis.title.x=element_blank(),plot.title=element_text(hjust=0),
        plot.title=element_text(size=14),axis.text=element_text(size=12),
        axis.text.x=element_blank(),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  scale_fill_manual(values=c("black","darkgray","white"))+
  annotate("text",x=1,y=1,label="Total Phosphorus",size=5)

bp2<- ggplot(INDICES,aes(x=GRP,y=L,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  labs(y="Relative Abundance TP Sensitive")+
  theme(legend.position="none",axis.title.x=element_blank(),plot.title=element_text(hjust=0),
        plot.title=element_text(size=14),axis.text=element_text(size=12),
        axis.text.x=element_blank(),axis.title.y=element_text(size=14))+
  scale_fill_manual(values=c("black","darkgray","white"))+
  annotate("text",x=1,y=1,label="Total Phosphorus",size=5)

bp3<- ggplot(INDICES,aes(x=GRP,y=R,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,10)+
  labs(y="Tolerant to Sensitive TP Index")+
  theme(legend.position="none",
        plot.title=element_text(hjust=0),
        plot.title=element_text(size=14),axis.text=element_text(size=12),
        axis.title.y=element_text(size=14),axis.title.x=element_blank())+
  scale_fill_manual(values=c("black","darkgray","white"))+
  annotate("text",x=1,y=10,label="Total Phosphorus",size=5)

bp4<- ggplot(pHINDICES,aes(x=GRP,y=H,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title=element_text(hjust=0),
        plot.title=element_text(size=14),axis.text=element_text(size=12),axis.text.y=element_blank(),
        axis.text.x=element_blank())+
  scale_fill_manual(values=c("black","darkgray","white"))+
  annotate("text",x=1,y=1,label="pH",size=5)

bp5<- ggplot(pHINDICES,aes(x=GRP,y=L,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title=element_text(hjust=0),
        plot.title=element_text(size=14),axis.text=element_text(size=12),axis.text.y=element_blank(),
        axis.text.x=element_blank())+
  scale_fill_manual(values=c("black","darkgray","white"))+
  annotate("text",x=1,y=1,label="pH",size=5)

bp6<- ggplot(pHINDICES,aes(x=GRP,y=R,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,10)+
  theme(legend.position="none",axis.title.y=element_blank(),
        plot.title=element_text(hjust=0),
        plot.title=element_text(size=14),axis.text=element_text(size=12),
        axis.text.y=element_blank(),axis.title.x=element_blank())+
  scale_fill_manual(values=c("black","darkgray","white"))+
  annotate("text",x=1,y=10,label="pH",size=5)

bp7<- ggplot(TINDICES,aes(x=GRP,y=H,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title=element_text(hjust=0),
        plot.title=element_text(size=14),axis.text=element_text(size=12),axis.text.y=element_blank(),
        axis.text.x=element_blank(),axis.title.x=element_text(size=14))+
  scale_fill_manual(values=c("black","darkgray","white"))+
  annotate("text",x=1,y=1,label="Temperature",size=5)

bp8<- ggplot(TINDICES,aes(x=GRP,y=L,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title=element_text(hjust=0),
        plot.title=element_text(size=14),axis.text=element_text(size=12),axis.text.y=element_blank(),
        axis.text.x=element_blank())+
  scale_fill_manual(values=c("black","darkgray","white"))+
  annotate("text",x=1,y=1,label="Temperature",size=5)

bp9<- ggplot(TINDICES,aes(x=GRP,y=R,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,10)+
  theme(legend.position="none",axis.title.y=element_blank(),
        plot.title=element_text(hjust=0),
        plot.title=element_text(size=14),axis.text=element_text(size=12),
        axis.text.y=element_blank(),axis.title.x=element_blank())+
  scale_fill_manual(values=c("black","darkgray","white"))+
  annotate("text",x=1,y=10,label="Temperature",size=5)

tiff(filename="DiatomRelAbundPlots.tiff",width=1100,height=1100)
multiplot(bp1,bp2,bp3,bp4,bp5,bp6,bp7,bp8,bp9,cols=3)
dev.off()



####Kruskal Wallis and Dunn Test for CT DEEP and NAWQA TP Groups############
###########################################################################

dunn.test(INDICES$R,INDICES$GRP,kw=TRUE)
dunn.test(INDICES$H,INDICES$GRP,kw=TRUE)
dunn.test(INDICES$L,INDICES$GRP,kw=TRUE)

#########Overall Percentiles of Metrics#######################
#############################################################

H_quant <- quantile(INDICES$H,c(0.1,0.25,0.5,0.75,0.9))
L_quant <- quantile(INDICES$L,c(0.1,0.25,0.5,0.75,0.9))
R_quant <- quantile(INDICES$R,c(0.1,0.25,0.5,0.75,0.9))

quant <- rbind(H_quant,L_quant,R_quant)
write.csv(quant,"DiatomMetricsPercentiles_101916.csv")


