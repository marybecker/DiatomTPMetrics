#######Calculate Tolerance Values for Individual Species#####
setwd ("")#SET WD

library(gam)

spp <- read.csv ("SPP.csv",header=TRUE,row.names=1)
spp[is.na(spp)] <- 0
TP <- read.csv ("TP.csv",header=TRUE,row.names=1)
TP$TP <- log(TP$TP_MGL+1)

taxa.names <- names(spp)

#############Define Curve Shape Function (Yuan, 2006)###################
################################################

unimod.test <- function(mnr, ubnd, lbnd) {
  
  # Find the maximum and minimum predicted mean probabilities
  lmax <- max(mnr)
  lmin <- min(mnr)
  
  # Find index locations for these probabilities
  imax <- match(lmax, mnr)
  imin <- match(lmin, mnr)
  
  x.out <- F
  y.out <- F
  
  # Compare mean predicted probability to the left of maximum point 
  # with upper confidence bound.  Store a T in x.out if
  # any point in the mean response deviates from the 
  # upper confidence limit
  if (imax > 1) {
    x.out <- sum(lmax == pmax(lmax, ubnd[1:(imax-1)])) > 0
  }
  
  # Store a T in y.out if any point in the mean probability
  # to the right of the maximum point deviates from the upper
  # confidence limit
  if (imax < length(ubnd)) {
    y.out <- sum(lmax == pmax(lmax, 
                              ubnd[(imax+1):length(ubnd)])) > 0
  }
  
  # Perform same set of tests for lower confidence limit
  a.out <- F
  b.out <- F
  
  if (imin > 1) {
    a.out <- sum(lmin == pmin(lmin, lbnd[1:(imin-1)])) > 0
  }
  if (imin > length(lbnd)) {
    b.out <- sum(lmin == pmin(lmin, 
                              lbnd[(imin+1):length(lbnd)])) > 0
  }
  
  # The information on where the mean curve deviates from the
  # confidence limits tells us its curve shape...
  if (x.out & y.out) {
    return("Unimodal")
  }
  if (a.out & b.out) {
    return("Concave up")
  }
  if (x.out | b.out) {
    return("Increasing")
  }
  if (y.out | a.out) {
    return("Decreasing")
  }
  if (! (x.out | y.out | a.out | b.out)) return(NA)
}


#######Whole Model Calc With Plots#############
######################################################

#Merge ENV & SPP Datasets#

SPP.TP <- merge (spp,TP,by=0)

# Create storage list for models

modlist.gam <- as.list(rep(NA, times = length(taxa.names)))

for (i in 1:length(taxa.names)) {
  # Create a logical vector is true if taxon is
  #   present and false if taxon is absent.
  resp <- SPP.TP[, taxa.names[i]] > 0
  
  # Fit the regression model, specifying two degrees of freedom
  # to the curve fit.
  modlist.gam[[i]] <- gam(resp ~ s(SPP.TP$TP, df = 2),
                          family = "binomial")
  
  print(summary(modlist.gam[[i]]))
}

for (i in 1:length(taxa.names)) {
  
  # Compute mean predicted probability of occurrence
  # and standard errors about this predicted probability.
  predres <- predict(modlist.gam[[i]], type= "link",se.fit=T)
  
  # Compute approximate upper and lower 90% confidence limits
  up.bound.link <- predres$fit + 1.65*predres$se.fit
  low.bound.link <- predres$fit - 1.65*predres$se.fit
  mean.resp.link <- predres$fit
  
  # Convert from logit transformed values to probability.
  up.bound <- exp(up.bound.link)/(1+exp(up.bound.link))
  low.bound <- exp(low.bound.link)/(1+exp(low.bound.link))
  mean.resp <- exp(mean.resp.link)/(1+exp(mean.resp.link))
  
  # Sort the environmental variable.
  iord <- order(SPP.TP$TP)
  
  # Define bins to summarize observational data as
  # probabilities of occurrence
  nbin <- 5
  
  # Define bin boundaries so each bin has approximately the same
  # number of observations.
  cutp <- quantile(SPP.TP$TP,
                   probs = seq(from = 0, to = 1, length = 20))
  
  # Compute the midpoint of each bin
  cutm <- 0.5*(cutp[-1] + cutp[-nbin])
  
  # Assign a factor to each bin
  cutf <- cut(SPP.TP$TP, cutp, include.lowest = T)
  
  # Compute the mean of the presence/absence data within each bin.
  vals <- tapply(SPP.TP[, taxa.names[i]] > 0, cutf, mean)
  
  # Now generate the plot
  # Plot binned observational data as symbols.
  
  file_name<- file.path("",##ADD FILE PATH HERE
                        paste("GAM","_",taxa.names[i],".tiff"))
  tiff(file=file_name,width=600,height=600,pointsize=20)
  
  plot(x=cutm, y=vals, xlab = "Log + 1 Total Phosphorus (mg/L)", 
       ylab = "Probability of occurrence", ylim = c(0,1), log="x",
       main = taxa.names[i],pch=16,col="black")		
  # Plot mean fit as a solid line.		
  lines(x= SPP.TP$TP[iord], y= mean.resp[iord],col="red",lwd=2)
  # Plot confidence limits as dotted lines.				
  lines (SPP.TP$TP[iord], up.bound[iord], lty = "dotted")
  lines (x=SPP.TP$TP[iord], y=low.bound[iord], lty = "dotted")
  
  
  dev.off()
}

chi <- as.vector(matrix(0.0,ncol=length(taxa.names),nrow=1))
# Conduct chi-square tests on nested parametric models
for (i in 1:length(taxa.names)) {
  
  print(taxa.names[i])
  resp <- SPP.TP[,taxa.names[i]] > 0
  modcmp <- gam(resp ~ 1, family = binomial, data = SPP.TP)
  modout <- anova(modlist.gam[[i]], modcmp, test = "Chi")
  print(modout)
  chi[i] <- modout[2,'Pr(>Chi)']
  if (modout[2,'Pr(>Chi)'] < 0.05) {
    print("Model significant compared to constant")
  }
}

names(chi)<- taxa.names
chi<-as.data.frame(chi)

#########ROC##########
roc <- rep(NA, times = length(taxa.names))

for (i in 1:length(taxa.names)) {
  # Compute mean predicted probability of occurrence
  predout <- predict(modlist.gam[[i]], type = "response")
  
  # Generate logical vector corresponding to presence/absence
  resp <- SPP.TP[, taxa.names[i]] > 0
  
  # Divide predicted probabilities into sites where
  # species is present ("x") and sites where the species is
  # absent ("y").
  x <- predout[resp]
  y <- predout[! resp]
  
  # Now perform all pairwise comparisons of x vs. y
  # and store results in a matrix
  rocmat <- matrix(NA, nrow = length(x), ncol = length(y))
  for (j in 1:length(x)) {
    rocmat[j,] <- as.numeric(x[j] > y)
  }
  
  # Summarize all comparisons to compute area under ROC
  roc[i] <- sum(rocmat)/(length(x)*length(y))
}

names(roc)<- taxa.names
roc<-as.data.frame(roc)

#########CURVE SHAPE#########

tolcl <- rep("", times = length(taxa.names)) 
for (i in 1:length(taxa.names)) {
  predres <- predict(modlist.gam[[i]], type= "link", se.fit = T)
  
  # Compute upper and lower 90% confidence limits
  up.bound.link <- predres$fit + 1.65*predres$se.fit
  low.bound.link <- predres$fit - 1.65*predres$se.fit
  mean.resp.link <- predres$fit
  
  # Convert from logit transformed values to probability.
  up.bound <- exp(up.bound.link)/(1+exp(up.bound.link))
  low.bound <- exp(low.bound.link)/(1+exp(low.bound.link))
  mean.resp <- exp(mean.resp.link)/(1+exp(mean.resp.link))
  
  # unimod.test requires that the responses be sorted by 
  # the value of the environmental variable.
  iord <- order(SPP.TP$TP)
  
  tolcl[i] <- unimod.test(mean.resp[iord], up.bound[iord], 
                          low.bound[iord])
}

names(tolcl)<- taxa.names
tolcl<-as.data.frame(tolcl)

gam_results<-cbind(chi,roc,tolcl)
gam_results$SigChi<-ifelse(gam_results$chi<=0.05,1,0)
gam_results$SigRoc<-ifelse(gam_results$roc>=0.6,1,0)
gam_results$Sig<-ifelse(gam_results$SigChi+gam_results$SigRoc==2,1,0)
write.table(gam_results,"gam_results.csv",sep=",",row.names=TRUE,col.names=NA)

SigSpp<-gam_results[which(gam_results$Sig==1),]
write.table(SigSpp,"SigSppGAM.csv",sep=",",row.names=TRUE,col.names=NA)
