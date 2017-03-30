#######Calculate Tolerance Values for Individual Species#####

setwd ("/Users/tbecker/Documents/Projects/GitHubProjects/DiatomTPMetrics")

library(gam)

spp <- read.csv ("data/SPP_032717.csv",header=TRUE,row.names=1)
TP <- read.csv ("data/TP.csv",header=TRUE,row.names=1)

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

#Merge ENV & SPP Datasets#

SPP.TP1 <- merge (spp,TP,by=0)
SPP.TP1<-SPP.TP1[,c(2:128)]

#########GAM & ROC Permutations##########
##############################################

n<- 1000
CURVEoutput<- matrix(nrow=length(taxa.names),ncol=n)
ROCoutput<- matrix(nrow=length(taxa.names),ncol=n)

for(k in 1:n) {

#Randomly sample rows from SPP.TP dataset
TPsamp <- sample(nrow(SPP.TP1),size=253)
SPP.TP <- SPP.TP1[TPsamp,]

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
}

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

tolcl_mat<- as.matrix(tolcl)
CURVEoutput[,k]<- tolcl_mat[,1]

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

roc_mat<- as.matrix(roc)
ROCoutput[,k]<- roc_mat[,1]

}

#######Sum Permutations and Export Tables#############
######################################################

IncrCnt<- rowSums(CURVEoutput=="Increasing")
DecrCnt<- rowSums(CURVEoutput=="Decreasing")
UniCnt<- rowSums(CURVEoutput=="Unimodal")
ROCCnt<- rowSums(ROCoutput>=0.6)
names(ROCCnt)<-taxa.names

write.table(ROCCnt,"ROC_SppGAM_032717.csv",sep=",",row.names=TRUE,col.names=NA)
write.table(IncrCnt,"Inc_SppGAM_032717.csv",sep=",",row.names=TRUE,col.names=NA)
write.table(DecrCnt,"Dec_SppGAM_032717.csv",sep=",",row.names=TRUE,col.names=NA)
write.table(UniCnt,"Uni_SppGAM_032717.csv",sep=",",row.names=TRUE,col.names=NA)

