setwd ("/Users/tbecker/Documents/Projects/GitHubProjects/DiatomTPMetrics")#SET WD#

INDICES <-  read.csv ("data/DiatomMetrics_040517.csv",sep=",",header=TRUE)
pHINDICES<- read.csv ("data/DiatomMetricspH.csv",sep=",",header=TRUE)
TINDICES<-  read.csv ("data/DiatomMetricsJulyTemp.csv",sep=",",header=TRUE)

TPQuant<- quantile(INDICES$TP_MGL,0.5)
INDICES$GRP <- ifelse(INDICES$TP_MGL > TPQuant,"H", "L")
pHQuant<- quantile(pHINDICES$pH,0.5)
pHINDICES$GRP <- ifelse(pHINDICES$pH > pHQuant,"H", "L")
TQuant<- quantile(TINDICES$Jtemp,0.5)
TINDICES$GRP <- ifelse(TINDICES$Jtemp > TQuant,"H", "L")


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


#####Calculate the discrimination efficiency of metrics#
VARDF<- INDICES  ##Change DF input to calculate for different dataset
VAR<- which(colnames(VARDF)=="TP_MGL")  ##Change column input to associate with correct dataset variable
QUANT<-quantile(VARDF[,VAR],0.5)
LGrp<- VARDF[which(VARDF[,VAR]<=QUANT),]
L75 <- quantile(LGrp$H,0.75)
L25 <- quantile(LGrp$L,0.25)
LR75 <- quantile(LGrp$R,0.75)
GrpTotal <- VARDF[which(VARDF[,VAR]>QUANT),]
n<- length(GrpTotal$ID)



#####Tolerant Metric DE
GrpSubset <- VARDF[which(VARDF[,VAR]>QUANT& VARDF$H >L75),]
n_subset <- length(GrpSubset$ID)
TolMetDE<- n_subset/n

#####Sensitive Metric DE
GrpSubset <- VARDF[which(VARDF[,VAR]>QUANT& VARDF$L <L25),]
n_subset <- length(GrpSubset$ID)
SenMetDE<-n_subset/n

#####Diatom Index DE
GrpSubset <- VARDF[which(VARDF[,VAR]>QUANT& VARDF$R >LR75),]
n_subset <- length(GrpSubset$ID)
DIDE<-n_subset/n


#####Run appropriate variable to save and compare after run above#####
TPDiatom<-rbind(TolMetDE,SenMetDE,DIDE)
pHDiatom<-rbind(TolMetDE,SenMetDE,DIDE)
TempDiatom<-rbind(TolMetDE,SenMetDE,DIDE)



#####Boxplots##
################################################

INDICES$GRP<- factor(INDICES$GRP,levels=c("L","H"))
pHINDICES$GRP<- factor(pHINDICES$GRP,levels=c("L","H"))
TINDICES$GRP<- factor(TINDICES$GRP,levels=c("L","H"))

bp1<- ggplot(INDICES,aes(x=GRP,y=H,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  labs(y="Relative Abundance TP Tolerant")+
  theme(legend.position="none",axis.title.x=element_blank(),
       # plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=30),
        axis.title.y=element_text(size=30))+
  scale_fill_manual(values=c("black","white"))+
  annotate("text",x=1,y=1,label="Total Phosphorus",size=10)#+
  #geom_hline(yintercept = quantile(INDICES$H,0.5),size=2,col="gray46")


bp2<- ggplot(INDICES,aes(x=GRP,y=L,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  labs(y="Relative Abundance TP Sensitive")+
  theme(legend.position="none",axis.title.x=element_blank(),
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.text=element_text(size=30),
        axis.text.x=element_blank(),axis.title.y=element_text(size=30))+
  scale_fill_manual(values=c("black","white"))+
  annotate("text",x=1,y=1,label="Total Phosphorus",size=10)#+
  #geom_hline(yintercept = quantile(INDICES$L,0.5),size=2,col="gray46")

bp3<- ggplot(INDICES,aes(x=GRP,y=R,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,10)+
  labs(y="Tolerant to Sensitive TP Index")+
  theme(legend.position="none",
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.text=element_text(size=30),
        axis.title.y=element_text(size=30),axis.title.x=element_blank())+
  scale_fill_manual(values=c("black","white"))+
  annotate("text",x=1,y=10,label="Total Phosphorus",size=10)#+
  #geom_hline(yintercept = quantile(INDICES$R,0.5),size=2,col="gray46")

bp4<- ggplot(pHINDICES,aes(x=GRP,y=H,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text=element_text(size=30),axis.text.y=element_blank(),
        axis.text.x=element_blank())+
  scale_fill_manual(values=c("black","white"))+
  annotate("text",x=1,y=1,label="pH",size=10)

bp5<- ggplot(pHINDICES,aes(x=GRP,y=L,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text=element_text(size=30),axis.text.y=element_blank(),
        axis.text.x=element_blank())+
  scale_fill_manual(values=c("black","white"))+
  annotate("text",x=1,y=1,label="pH",size=10)

bp6<- ggplot(pHINDICES,aes(x=GRP,y=R,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,10)+
  theme(legend.position="none",axis.title.y=element_blank(),
        axis.text=element_text(size=30),
        axis.text.y=element_blank(),axis.title.x=element_blank())+
  scale_fill_manual(values=c("black","white"))+
  annotate("text",x=1,y=10,label="pH",size=10)

bp7<- ggplot(TINDICES,aes(x=GRP,y=H,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.y=element_blank())+
  scale_fill_manual(values=c("black","white"))+
  annotate("text",x=1,y=1,label="Temperature",size=10)

bp8<- ggplot(TINDICES,aes(x=GRP,y=L,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text=element_text(size=30),axis.text.y=element_blank(),
        axis.text.x=element_blank())+
  scale_fill_manual(values=c("black","white"))+
  annotate("text",x=1,y=1,label="Temperature",size=10)

bp9<- ggplot(TINDICES,aes(x=GRP,y=R,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,10)+
  theme(legend.position="none",axis.title.y=element_blank(),
        axis.text=element_text(size=30),
        axis.text.y=element_blank(),axis.title.x=element_blank())+
  scale_fill_manual(values=c("black","white"))+
  annotate("text",x=1,y=10,label="Temperature",size=10)

tiff(file="DiatomRelAbundPlots.tiff",width=2400,height=1800)
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


