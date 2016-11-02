setwd (#SET WD)

INDICES <- read.csv ("data/DiatomMetricsChloride.csv",sep=",",header=TRUE)


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

######Calculates that Discrimination Efficiency for Tolerant Metric for different combinations of CL Groups#
#############################################
 LGrp <- INDICES[which(INDICES$CL_GRP=='L'),]
 L75 <- quantile(LGrp$H,0.75)
 MGrp <- INDICES[which(INDICES$CL_GRP=='M'),]
 M75 <- quantile(MGrp$H,0.75)

#LH DE#
GrpTotal <- INDICES[which(INDICES$CL_GRP=='H'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$CL_GRP=='H'& INDICES$H >L75),]
n_subset <- length(GrpSubset$ID)
n_subset/n

#LM DE#
GrpTotal <- INDICES[which(INDICES$CL_GRP=='M'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$CL_GRP=='M'& INDICES$H >L75),]
n_subset <- length(GrpSubset$ID)
n_subset/n

#MH#
GrpTotal <- INDICES[which(INDICES$CL_GRP=='H'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$CL_GRP=='H'& INDICES$H >M75),]
n_subset <- length(GrpSubset$ID)
n_subset/n

#####Calculates that Discrimination Efficiency for Sensitive metric for different combinations of CL Groups##
################################################
 LGrp <- INDICES[which(INDICES$CL_GRP=='L'),]
 L25 <- quantile(LGrp$L,0.25)
 MGrp <- INDICES[which(INDICES$CL_GRP=='M'),]
 M25 <- quantile(MGrp$L,0.25)

#LH DE#
GrpTotal <- INDICES[which(INDICES$CL_GRP=='H'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$CL_GRP=='H'& INDICES$L <L25),]
n_subset <- length(GrpSubset$ID)
n_subset/n

#LM DE#
GrpTotal <- INDICES[which(INDICES$CL_GRP=='M'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$CL_GRP=='M'& INDICES$L <L25),]
n_subset <- length(GrpSubset$ID)
n_subset/n

#MH#
GrpTotal <- INDICES[which(INDICES$CL_GRP=='H'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$CL_GRP=='H'& INDICES$L <M25),]
n_subset <- length(GrpSubset$ID)
n_subset/n

#####Calculates that Discrimination Efficiency for T to S Ratio metric for different combinations of CL Groups##
################################################
LGrp <- INDICES[which(INDICES$CL_GRP=='L'),]
L75 <- quantile(LGrp$R,0.75)
MGrp <- INDICES[which(INDICES$CL_GRP=='M'),]
M75 <- quantile(MGrp$R,0.75)

#LH DE#
GrpTotal <- INDICES[which(INDICES$CL_GRP=='H'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$CL_GRP=='H'& INDICES$R >L75),]
n_subset <- length(GrpSubset$ID)
n_subset/n

#LM DE#
GrpTotal <- INDICES[which(INDICES$CL_GRP=='M'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$CL_GRP=='M'& INDICES$R >L75),]
n_subset <- length(GrpSubset$ID)
n_subset/n

#MH#
GrpTotal <- INDICES[which(INDICES$CL_GRP=='H'),]
n<- length(GrpTotal$ID)
GrpSubset <- INDICES[which(INDICES$CL_GRP=='H'& INDICES$R >M75),]
n_subset <- length(GrpSubset$ID)
n_subset/n


#####Boxplots and Scatterplots##
################################################
INDICES$CL_GRP<- factor(INDICES$CL_GRP,levels=c("L","M","H"))

p1<-ggplot(INDICES,aes(x=CL,y=H))+
  geom_point()+
  geom_smooth(se=FALSE,colour="black")+
  ylim(0,1)+
  labs(title="A",y="Relative Abundance of Tolerant Diatoms")+
  theme(legend.position="none",axis.title.x=element_blank(),plot.title=element_text(hjust=0),
        plot.title=element_text(size=14),axis.text=element_text(size=12),
        axis.text.x=element_blank(),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))

p2<-ggplot(INDICES,aes(x=CL,y=L))+
  geom_point()+
  geom_smooth(se=FALSE,colour="black")+
  ylim(0,1)+
  labs(title="B",y="Relative Abundance of Sensitive Diatoms")+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_text(size=14),plot.title=element_text(hjust=0),
        plot.title=element_text(size=14),axis.text=element_text(size=12),
        axis.text.x=element_blank())

p3<-ggplot(INDICES,aes(x=CL,y=R))+
  geom_point()+
  geom_smooth(se=FALSE,colour="black")+
  ylim(0,10)+
  labs(title="C",y="Tolerant to Sensitive Diatom Index",x="Total Chloride mg/L")+
  theme(legend.position="none",axis.title.y=element_text(size=14),axis.title.x=element_text(size=14),
        plot.title=element_text(hjust=0),plot.title=element_text(size=14),axis.text=element_text(size=12))

bp4<- ggplot(INDICES,aes(x=CL_GRP,y=H,fill=CL_GRP))+
  geom_boxplot(aes(fill=CL_GRP))+
  ylim(0,1)+
  labs(title="D")+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title=element_text(hjust=0),
        plot.title=element_text(size=14),axis.text=element_text(size=12),axis.text.y=element_blank(),
        axis.text.x=element_blank(),axis.title.x=element_text(size=14))+
  scale_fill_manual(values=c("black","darkgray","white"))

bp5<- ggplot(INDICES,aes(x=CL_GRP,y=L,fill=CL_GRP))+
  geom_boxplot(aes(fill=CL_GRP))+
  ylim(0,1)+
  labs(title="E")+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title=element_text(hjust=0),
        plot.title=element_text(size=14),axis.text=element_text(size=12),axis.text.y=element_blank(),
        axis.text.x=element_blank())+
  scale_fill_manual(values=c("black","darkgray","white"))

bp6<- ggplot(INDICES,aes(x=CL_GRP,y=R,fill=CL_GRP))+
  geom_boxplot(aes(fill=CL_GRP))+
  ylim(0,10)+
  labs(title="F",x="Total Chloride Group")+
  theme(legend.position="none",axis.title.y=element_blank(),axis.title.x=element_text(size=14),
        plot.title=element_text(hjust=0),
        plot.title=element_text(size=14),axis.text=element_text(size=12),
        axis.text.y=element_blank())+
  scale_fill_manual(values=c("black","darkgray","white"))

tiff(filename="DiatomChloridePlots101816.tiff",width=1100,height=1100)
multiplot(p1,p2,p3,bp4,bp5,bp6,cols=2)
dev.off()

####Kruskal Wallis and Dunn Test for CT DEEP and NAWQA CL Groups############
###########################################################################

dunn.test(INDICES$R,INDICES$CL_GRP,kw=TRUE)
dunn.test(INDICES$H,INDICES$CL_GRP,kw=TRUE)
dunn.test(INDICES$L,INDICES$CL_GRP,kw=TRUE)

