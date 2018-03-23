setwd ("")#SET WD#

INDICES <-  read.csv ("DiatomMetrics_TESTCHEM.csv",sep=",",header=TRUE)
INDICES[is.na(INDICES)] <- 0
INDICES$R<-(INDICES$H*10)/(INDICES$H+INDICES$L)
INDICES_GAM <-  read.csv ("DiatomMetrics_TESTGAM.csv",sep=",",header=TRUE)
INDICES_GAM[is.na(INDICES_GAM)] <- 0
INDICES_GAM$R<-(INDICES_GAM$H*10)/(INDICES_GAM$H+INDICES_GAM$L)
INDICES_ANSPR<- read.csv("DiatomMetrics_TESTANSPRegional.csv",sep=",",header=TRUE)
INDICES_ANSPR[is.na(INDICES_ANSPR)] <- 0
INDICES_ANSPR$R<-(INDICES_ANSPR$H*10)/(INDICES_ANSPR$H+INDICES_ANSPR$L)
INDICES_ANSP<- read.csv("DiatomMetrics_TESTANSP.csv",sep=",",header=TRUE)
INDICES_ANSP[is.na(INDICES_ANSP)] <- 0
INDICES_ANSP$R<-(INDICES_ANSP$H*10)/(INDICES_ANSP$H+INDICES_ANSP$L)
INDICES_JTEMP<- read.csv("DiatomMetrics_TESTJTEMP.csv",sep=",",header=TRUE)
INDICES_JTEMP[is.na(INDICES_JTEMP)] <- 0
INDICES_JTEMP$R<-(INDICES_JTEMP$H*10)/(INDICES_JTEMP$H+INDICES_JTEMP$L)

#####Calculate the discrimination efficiency of metrics#####
############################################################

VARDF<- INDICES_GAM  ##Change DF input to calculate for different dataset
VAR<- which(colnames(VARDF)=="TP_MGL")  ##Change column input to associate with correct dataset variable
QUANT<- 0.065##Comment out when variable is something other than TP
##QUANT<-quantile(VARDF[,VAR],0.75)##Use when variable is something other than TP
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
IndValMet<-rbind(TolMetDE,SenMetDE,DIDE)
GAMMet<-rbind(TolMetDE,SenMetDE,DIDE)
ANSPReg<-rbind(TolMetDE,SenMetDE,DIDE)
ANSPMet<-rbind(TolMetDE,SenMetDE,DIDE)

######Spearman correlations#########################################
chem<- INDICES_JTEMP[,c("ID","Jtemp")]
chem<- merge(chem,INDICES_ANSPR,by="ID")
VARDF <- chem
VAR<- "Jtemp"

corH<- cor(VARDF[,VAR],VARDF$H,method="spearman")
corS<- cor(VARDF[,VAR],VARDF$L,method="spearman")
corR<- cor(VARDF[,VAR],VARDF$R,method="spearman")

cor<-rbind(corH,corS,corR)


VARDF<- INDICES_ANSPR
wilcox.test(H~GRP, data=VARDF)
wilcox.test(L~GRP, data=VARDF)
wilcox.test(R~GRP, data=VARDF)


#####Boxplots##
################################################
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

##################################################################################
INDICES$GRP<- ifelse(INDICES$TP_MGL<=0.065,"L","H")
INDICES$GRP<- factor(INDICES$GRP,levels=c("L","H"))
INDICES_GAM$GRP<- ifelse(INDICES_GAM$TP_MGL<=0.065,"L","H")
INDICES_GAM$GRP<- factor(INDICES_GAM$GRP,levels=c("L","H"))
INDICES_ANSPR$GRP<- ifelse(INDICES_ANSPR$TP_MGL<=0.065,"L","H")
INDICES_ANSPR$GRP<- factor(INDICES_ANSPR$GRP,levels=c("L","H"))
INDICES_ANSP$GRP<- ifelse(INDICES_ANSP$TP_MGL<=0.065,"L","H")
INDICES_ANSP$GRP<- factor(INDICES_ANSP$GRP,levels=c("L","H"))

bp1<- ggplot(INDICES_GAM,aes(x=GRP,y=H,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  labs(y="Relative Abundance TP Tolerant")+
  theme(legend.position="none",axis.title.x=element_blank(),
       # plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=7),
        axis.title.y=element_text(size=7),
       panel.background = element_rect(fill = NA),
       panel.grid.major = element_line(colour = NA),
       panel.border=element_rect(colour = "grey85",fill=NA))+
  scale_fill_manual(values=c("black","grey50"))+
  annotate("text",x=1,y=1,label="CT GAM",size=3)#+
  #geom_hline(yintercept = quantile(INDICES$H,0.5),size=2,col="gray46")


bp2<- ggplot(INDICES_GAM,aes(x=GRP,y=L,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  labs(y="Relative Abundance TP Sensitive")+
  theme(legend.position="none",axis.title.x=element_blank(),
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.text=element_text(size=7),
        axis.text.x=element_blank(),axis.title.y=element_text(size=7),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  scale_fill_manual(values=c("black","grey50"))+
  annotate("text",x=1,y=1,label="CT GAM",size=3)#+
  #geom_hline(yintercept = quantile(INDICES$L,0.5),size=2,col="gray46")

bp3<- ggplot(INDICES_GAM,aes(x=GRP,y=R,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,10)+
  labs(y="Tolerant to Sensitive TP Index")+
  theme(legend.position="none",
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.text=element_text(size=7),
        axis.title.y=element_text(size=7),axis.title.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  scale_fill_manual(values=c("black","grey50"))+
  annotate("text",x=1,y=10,label="CT GAM",size=3)#+
  #geom_hline(yintercept = quantile(INDICES$R,0.5),size=2,col="gray46")

bp4<- ggplot(INDICES,aes(x=GRP,y=H,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        # plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  scale_fill_manual(values=c("black","grey50"))+
  annotate("text",x=1,y=1,label="CT IndVal WA",size=3)#+
#geom_hline(yintercept = quantile(INDICES$H,0.5),size=2,col="gray46")


bp5<- ggplot(INDICES,aes(x=GRP,y=L,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  scale_fill_manual(values=c("black","grey50"))+
  annotate("text",x=1,y=1,label="CT IndVal WA",size=3)#+
#geom_hline(yintercept = quantile(INDICES$L,0.5),size=2,col="gray46")

bp6<- ggplot(INDICES,aes(x=GRP,y=R,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,10)+
  theme(legend.position="none",
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  scale_fill_manual(values=c("black","grey50"))+
  annotate("text",x=1,y=10,label="CT IndVal WA",size=3)#+
#geom_hline(yintercept = quantile(INDICES$R,0.5),size=2,col="gray46")

bp7<- ggplot(INDICES_ANSP,aes(x=GRP,y=H,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),
        # plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  scale_fill_manual(values=c("black","grey50"))+
  annotate("text",x=1,y=1,label="National",size=3)#+
#geom_hline(yintercept = quantile(INDICES$H,0.5),size=2,col="gray46")


bp8<- ggplot(INDICES_ANSP,aes(x=GRP,y=L,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  scale_fill_manual(values=c("black","grey50"))+
  annotate("text",x=1,y=1,label="National",size=3)#+
#geom_hline(yintercept = quantile(INDICES$L,0.5),size=2,col="gray46")

bp9<- ggplot(INDICES_ANSP,aes(x=GRP,y=R,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,10)+
  theme(legend.position="none",
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  scale_fill_manual(values=c("black","grey50"))+
  annotate("text",x=1,y=10,label="National",size=3)#+
#geom_hline(yintercept = quantile(INDICES$R,0.5),size=2,col="gray46")

bp10<- ggplot(INDICES_ANSPR,aes(x=GRP,y=H,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),
        # plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  scale_fill_manual(values=c("black","grey50"))+
  annotate("text",x=1,y=1,label="Regional",size=3)#+
#geom_hline(yintercept = quantile(INDICES$H,0.5),size=2,col="gray46")

bp11<- ggplot(INDICES_ANSPR,aes(x=GRP,y=L,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  scale_fill_manual(values=c("black","grey50"))+
  annotate("text",x=1,y=1,label="Regional",size=3)#+
#geom_hline(yintercept = quantile(INDICES$L,0.5),size=2,col="gray46")

bp12<- ggplot(INDICES_ANSPR,aes(x=GRP,y=R,fill=GRP))+
  geom_boxplot(aes(fill=GRP))+
  ylim(0,10)+
  theme(legend.position="none",
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  scale_fill_manual(values=c("black","grey50"))+
  annotate("text",x=1,y=10,label="Regional",size=3)#+
#geom_hline(yintercept = quantile(INDICES$R,0.5),size=2,col="gray46")

tiff(file="DiatomRelAbundPlots.tiff",width=190,height=142.4,units="mm",
     pointsize=1/300,res=300)
multiplot(bp1,bp2,bp3,bp4,bp5,bp6,bp7,bp8,bp9,bp10,bp11,bp12,cols=4)
dev.off()

#########Overall Percentiles of Metrics#######################
#############################################################

H_quant <- quantile(INDICES_GAM$H,c(0.1,0.25,0.5,0.75,0.9))
L_quant <- quantile(INDICES_GAM$L,c(0.1,0.25,0.5,0.75,0.9))
R_quant <- quantile(INDICES_GAM$R,c(0.1,0.25,0.5,0.75,0.9))

quant <- rbind(H_quant,L_quant,R_quant)
write.csv(quant,"DiatomMetricsPercentiles_101916.csv")

########SCATTERPLOTS#####################################
##########################################################
sp1<- ggplot(INDICES_GAM,aes(x=TP_MGL,y=H))+
  geom_point()+
  geom_smooth(se=FALSE,colour="black")+
  ylim(0,1)+
  scale_x_log10()+
  labs(y="Relative Abundance TP Tolerant")+
  theme(legend.position="none",axis.title.x=element_blank(),
        # plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=7),
        axis.title.y=element_text(size=7),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  scale_fill_manual(values=c("black","white"))+
  annotate("text",x=0.025,y=0.875,label="TP (mg/L)",size=3)#+



sp2<- ggplot(INDICES_GAM,aes(x=TP_MGL,y=L))+
  geom_point()+
  geom_smooth(se=FALSE,colour="black")+
  scale_x_log10()+
  ylim(0,1)+
  labs(y="Relative Abundance TP Sensitive")+
  theme(legend.position="none",axis.title.x=element_blank(),
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.text=element_text(size=7),
        axis.text.x=element_blank(),axis.title.y=element_text(size=7),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  scale_fill_manual(values=c("black","white"))+
  annotate("text",x=0.025,y=0.875,label="TP (mg/L)",size=3)#+


sp3<- ggplot(INDICES_GAM,aes(x=TP_MGL,y=R))+
  geom_point()+
  geom_smooth(se=FALSE,colour="black")+
  scale_x_log10()+
  ylim(0,10)+
  labs(y="Tolerant to Sensitive TP Index")+
  theme(legend.position="none",
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.text=element_text(size=7),
        axis.title.y=element_text(size=7),axis.title.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  scale_fill_manual(values=c("black","white"))+
  annotate("text",x=0.025,y=8.75,label="TP (mg/L)",size=3)#+


sp4<- ggplot(INDICES_GAM,aes(x=Chloride,y=H))+
  geom_point()+
  geom_smooth(se=FALSE,colour="black")+
  scale_x_log10()+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        # plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  annotate("text",x=5,y=0.875,label="Chloride (mg/L)",size=3)#+



sp5<- ggplot(INDICES_GAM,aes(x=Chloride,y=L))+
  geom_point()+
  geom_smooth(se=FALSE,colour="black")+
  scale_x_log10()+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  annotate("text",x=5,y=0.875,label="Chloride (mg/L)",size=3)#+

sp6<- ggplot(INDICES_GAM,aes(x=Chloride,y=R))+
  geom_point()+
  geom_smooth(se=FALSE,colour="black")+
  scale_x_log10()+
  ylim(0,10)+
  theme(legend.position="none",
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  annotate("text",x=5,y=8.75,label="Chloride (mg/L)",size=3)#+


sp7<- ggplot(INDICES_GAM,aes(x=pH,y=H))+
  geom_point()+
  geom_smooth(se=FALSE,colour="black")+
  #scale_x_log10()+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),
        # plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  annotate("text",x=6,y=0.875,label="pH",size=3)#+



sp8<- ggplot(INDICES_GAM,aes(x=pH,y=L))+
  geom_point()+
  geom_smooth(se=FALSE,colour="black")+
  #scale_x_log10()+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  annotate("text",x=6,y=0.875,label="pH",size=3)#+


sp9<- ggplot(INDICES_GAM,aes(x=pH,y=R))+
  geom_point()+
  geom_smooth(se=FALSE,colour="black")+
  ylim(0,10)+
  theme(legend.position="none",
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  annotate("text",x=6,y=8.75,label="pH",size=3)#+


sp10<- ggplot(INDICES_JTEMP,aes(x=Jtemp,y=H))+
  geom_point()+
  geom_smooth(se=FALSE,colour="black")+
  #scale_x_log10()+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  annotate("text",x=18.5,y=0.875, parse=TRUE,
           label=as.character(expression(paste("July Temp(",degree ~ C,")"))),size=3)


sp11<- ggplot(INDICES_JTEMP,aes(x=Jtemp,y=L))+
  geom_point()+
  geom_smooth(se=FALSE,colour="black")+
  #scale_x_log10()+
  ylim(0,1)+
  theme(legend.position="none",axis.title.x=element_blank(),
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  annotate("text",x=18.5,y=0.875, parse=TRUE,
           label=as.character(expression(paste("July Temp(",degree ~ C,")"))),size=3)


sp12<- ggplot(INDICES_JTEMP,aes(x=Jtemp,y=R))+
  geom_point()+
  geom_smooth(se=FALSE,colour="black")+
  ylim(0,10)+
  theme(legend.position="none",
        #plot.title=element_text(hjust=0),plot.title=element_text(size=10),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.border=element_rect(colour = "grey85",fill=NA))+
  annotate("text",x=18.5,y=8.75, parse=TRUE,
           label=as.character(expression(paste("July Temp(",degree ~ C,")"))),size=3)

tiff(file="DiatomMetricScatterPlots.tiff",width=190,height=142.4,units="mm",
     pointsize=1/300,res=300)
multiplot(sp1,sp2,sp3,sp4,sp5,sp6,sp7,sp8,sp9,sp10,sp11,sp12,cols=4)
dev.off()