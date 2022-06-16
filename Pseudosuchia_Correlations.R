library(BAMMtools)
library(ggplot2)
library(adephylo)
library(pracma)
library(caper)
library(phangorn)
library(phytools)

# set files here
tree <- read.tree("FBD_Pseudosuchia_041219.tre")
bamm_data <- "PseudosuchiaNov2020_event_data.txt"
habitat="Marine" #set to whichever partition is required
df=10

#set variables for correlation
delo18 <- read.csv("science_delO18_curve.csv", header=FALSE)
haq <- read.csv("haq_sealevel_curve.csv",header=TRUE)
LTT <- read.csv('coords.csv', header=TRUE)

# Sort out temperature curves
# remove some autocorrelation
finaltemp = smooth(smooth(delo18$V2))
# make the correlation to temperature (inverse of DelO18)
finaltemp = -finaltemp

#DCCA-based LRCC test (giving the test statistic, standard error and $p$-value)
DCCA_LRCC<-function(x,y,s){
  xx<-cumsum(x)
  yy<-cumsum(y)
  t<-1:length(xx)
  F2sj_xy<-runif(floor(length(xx)/s))
  F2sj_xx<-F2sj_xy
  F2sj_yy<-F2sj_xy
  for(ss in seq(1,(floor(length(xx)/s)*s),by=s)){
    F2sj_xy[(ss-1)/s+1]<-sum((summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals))/(s-1)
    F2sj_xx[(ss-1)/s+1]<-sum((summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals))/(s-1)
    F2sj_yy[(ss-1)/s+1]<-sum((summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals))/(s-1)
  }
  rho<-mean(F2sj_xy)/sqrt(mean(F2sj_xx)*mean(F2sj_yy))
  return(c(rho,1/sqrt(length(xx)),1-pnorm(abs(rho),mean=0,sd=1/sqrt(length(xx)))))
}

# sort out your tree and load in BAMM data
edata <- getEventData(tree, eventdata = bamm_data, burnin=0.1)
lengths = distRoot(tree, tree$tip.label)
root = max(lengths)
count = 1

# load in habitat data
habitatdata <- read.csv("Marine.txt", header=F, stringsAsFactors = FALSE)

# get habitat data
subtreetaxa <- habitatdata$V1

# work out where this node lives in tree age in Ma
H<-nodeHeights(tree)
nn<-findMRCA(tree,subtreetaxa)
h<-H[tree$edge==nn][1]
clade_root = root - h

# get subtree
stree <- subtreeBAMM(edata, tips=subtreetaxa)

# get rate through time for subtree
rtt <- getRateThroughTimeMatrix(stree)
# the subtree gets shifted to an arbitrary zero, so 
# we need to turn time around (from Myr to Ma) and
# then shift to the correct node age
times = abs(rtt$times-max(rtt$times))
times = times + (clade_root - max(rtt$times))

#times = times + 6.2895 #for marine
#times = times + 0.0059 #for terrestrial

numberOfSims = length(rtt$lambda)/length(rtt$times)

## SPECIATION ##

# This is where we store the correlation numbers for each simulation in BAMM
cors_temp <- rep(NA, numberOfSims)
cors_haq <- rep(NA, numberOfSims)
cors_LTT <- rep(NA, numberOfSims)

# Speciation Correlation for Temperature
for (i in 1:numberOfSims ) {
    interpdiv = approx(times, rtt$lambda[i,], delo18$V1, method='linear', rule=1)
    end = which(is.na(interpdiv$y))
    if (length(end) == 0) {
        div_rates = interpdiv$y
        ft = finaltemp
    } else {
        div_rates = interpdiv$y[-end]
        ft = finaltemp[-end]
    }
    t = DCCA_LRCC(as.numeric(unlist(div_rates)),as.numeric(unlist(ft)),length(ft)/df)
    cors_temp[i] = t[1]
}

plot=qplot(cors_temp, geom="histogram")
plot2 = plot + geom_histogram(fill="lightgrey") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_vline(aes(xintercept=mean(cors_temp)+sd(cors_temp)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_temp)-sd(cors_temp)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_temp)), linetype="solid") + ylab("Frequency") + xlab("Correlation co-efficient")
ggsave(plot2,file=paste0("Speciation","_",habitat,"_","Temperature",".pdf"))
sink(file=paste0("Speciation","_",habitat,"_","Temperature",".txt"))
print( summary(cors_temp) )
print( quantile(cors_temp, c(0.025, 0.975)) )
print (sd(cors_temp) )
print(wilcox.test(cors_temp, mu=0.0))
sink()

# Speciation Correlation for Haq Sea Level
for (i in 1:numberOfSims ) {
  interpdiv = approx(times, rtt$lambda[i,], haq$Age, method='linear', rule=1)
  end = which(is.na(interpdiv$y))
  if (length(end) == 0) {
    div_rates = interpdiv$y
    ft = haq$SL
  } else {
    div_rates = interpdiv$y[-end]
    ft = haq$SL[-end]
  }
  t = DCCA_LRCC(as.numeric(unlist(div_rates)),as.numeric(unlist(ft)),length(ft)/df)
  cors_haq[i] = t[1]
}

plot=qplot(cors_haq, geom="histogram")
plot2 = plot + geom_histogram(fill="lightgrey") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_vline(aes(xintercept=mean(cors_haq)+sd(cors_haq)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_haq)-sd(cors_haq)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_haq)), linetype="solid") + ylab("Frequency") + xlab("Correlation co-efficient")
ggsave(plot2,file=paste0("Speciation","_",habitat,"_","HaqSeaLevel",".pdf"))
sink(file=paste0("Speciation","_",habitat,"_","HaqSeaLevel",".txt"))
print( summary(cors_haq) )
print( quantile(cors_haq, c(0.025, 0.975)) )
print (sd(cors_haq) )
print(wilcox.test(cors_haq, mu=0.0))
sink()

# Speciation correlation for LTT

for (i in 1:numberOfSims ) {
    interpdiv = approx(times, rtt$lambda[i,], LTT$time, method='linear', rule=1)
    end = which(is.na(interpdiv$y))
    if (length(end) == 0) {
        div_rates = interpdiv$y
        ft = LTT$N
    } else {
        div_rates = interpdiv$y[-end]
        ft = LTT$N[-end]
    }
    t = DCCA_LRCC(as.numeric(unlist(div_rates)),as.numeric(unlist(ft)),length(ft)/df)
    cors_LTT[i] = t[1]
}

plot=qplot(cors_LTT, geom="histogram")
plot2 = plot + geom_histogram(fill="lightgrey") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_vline(aes(xintercept=mean(cors_LTT)+sd(cors_LTT)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_LTT)-sd(cors_LTT)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_LTT)), linetype="solid") + ylab("Frequency") + xlab("Correlation co-efficient")
ggsave(plot2,file=paste0("Speciation","_",habitat,"_","LTT",".pdf"))
sink(file=paste0("Speciation","_",habitat,"_","LTT",".txt"))
print( summary(cors_LTT) )
print( quantile(cors_LTT, c(0.025, 0.975)) )
print (sd(cors_LTT) )
print(wilcox.test(cors_LTT, mu=0.0))
sink()

## EXTINCTION ##

# This is where we store the correlation numbers for each simulation in BAMM
cors_temp <- rep(NA, numberOfSims)
cors_haq <- rep(NA, numberOfSims)
cors_LTT <- rep(NA, numberOfSims)

# Extinction Correlation for Temperature
for (i in 1:numberOfSims ) {
  interpdiv = approx(times, rtt$mu[i,], delo18$V1, method='linear', rule=1)
  end = which(is.na(interpdiv$y))
  if (length(end) == 0) {
    div_rates = interpdiv$y
    ft = finaltemp
  } else {
    div_rates = interpdiv$y[-end]
    ft = finaltemp[-end]
  }
  t = DCCA_LRCC(as.numeric(unlist(div_rates)),as.numeric(unlist(ft)),length(ft)/df)
    cors_temp[i] = t[1]
}

plot=qplot(cors_temp, geom="histogram")
plot2 = plot + geom_histogram(fill="lightgrey") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_vline(aes(xintercept=mean(cors_temp)+sd(cors_temp)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_temp)-sd(cors_temp)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_temp)), linetype="solid") + ylab("Frequency") + xlab("Correlation co-efficient")
ggsave(plot2,file=paste0("Extinction","_",habitat,"_","Temperature",".pdf"))
sink(file=paste0("Extinction","_",habitat,"_","Temperature",".txt"))
print( summary(cors_temp) )
print( quantile(cors_temp, c(0.025, 0.975)) )
print (sd(cors_temp) )
print(wilcox.test(cors_temp, mu=0.0))
sink()

# Extinction Correlation for Haq Sea Level
for (i in 1:numberOfSims ) {
  interpdiv = approx(times, rtt$mu[i,], haq$Age, method='linear', rule=1)
  end = which(is.na(interpdiv$y))
  if (length(end) == 0) {
    div_rates = interpdiv$y
    ft = haq$SL
  } else {
    div_rates = interpdiv$y[-end]
    ft = haq$SL[-end]
  }
  t = DCCA_LRCC(as.numeric(unlist(div_rates)),as.numeric(unlist(ft)),length(ft)/df)
  cors_haq[i] = t[1]
}

plot=qplot(cors_haq, geom="histogram")
plot2 = plot + geom_histogram(fill="lightgrey") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_vline(aes(xintercept=mean(cors_haq)+sd(cors_haq)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_haq)-sd(cors_haq)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_haq)), linetype="solid") + ylab("Frequency") + xlab("Correlation co-efficient")
ggsave(plot2,file=paste0("Extinction","_",habitat,"_","HaqSeaLevel",".pdf"))
sink(file=paste0("Extinction","_",habitat,"_","HaqSeaLevel",".txt"))
print( summary(cors_haq) )
print( quantile(cors_haq, c(0.025, 0.975)) )
print (sd(cors_haq) )
print(wilcox.test(cors_haq, mu=0.0))
sink()

# Extinction Correlation for LTT
for (i in 1:numberOfSims ) {
  interpdiv = approx(times, rtt$mu[i,], LTT$time, method='linear', rule=1)
  end = which(is.na(interpdiv$y))
  if (length(end) == 0) {
    div_rates = interpdiv$y
    ft = LTT$N
  } else {
    div_rates = interpdiv$y[-end]
    ft = LTT$N[-end]
  }
  t = DCCA_LRCC(as.numeric(unlist(div_rates)),as.numeric(unlist(ft)),length(ft)/df)
  cors_LTT[i] = t[1]
}

plot=qplot(cors_LTT, geom="histogram")
plot2 = plot + geom_histogram(fill="lightgrey") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_vline(aes(xintercept=mean(cors_LTT)+sd(cors_LTT)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_LTT)-sd(cors_LTT)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_LTT)), linetype="solid") + ylab("Frequency") + xlab("Correlation co-efficient")
ggsave(plot2,file=paste0("Extinction","_",habitat,"_","LTT",".pdf"))
sink(file=paste0("Extinction","_",habitat,"_","LTT",".txt"))
print( summary(cors_LTT) )
print( quantile(cors_LTT, c(0.025, 0.975)) )
print (sd(cors_LTT) )
print(wilcox.test(cors_LTT, mu=0.0))
sink()

## NET DIVERSIFICATION ##

# This is where we store the correlation numbers for each simulation in BAMM
cors_temp <- rep(NA, numberOfSims)
cors_haq <- rep(NA, numberOfSims)
cors_LTT <- rep(NA, numberOfSims)

netdiv<-(rtt$lambda)-(rtt$mu)

# NetDiversification Correlation for Temperature
for (i in 1:numberOfSims ) {
  interpdiv = approx(times, netdiv[i,], delo18$V1, method='linear', rule=1)
  end = which(is.na(interpdiv$y))
  if (length(end) == 0) {
    div_rates = interpdiv$y
    ft = finaltemp
  } else {
    div_rates = interpdiv$y[-end]
    ft = finaltemp[-end]
  }
  t = DCCA_LRCC(as.numeric(unlist(div_rates)),as.numeric(unlist(ft)),length(ft)/df)
  cors_temp[i] = t[1]
}

plot=qplot(cors_temp, geom="histogram")
plot2 = plot + geom_histogram(fill="lightgrey") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_vline(aes(xintercept=mean(cors_temp)+sd(cors_temp)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_temp)-sd(cors_temp)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_temp)), linetype="solid") + ylab("Frequency") + xlab("Correlation co-efficient")
ggsave(plot2,file=paste0("NetDiversification","_",habitat,"_","Temperature",".pdf"))
sink(file=paste0("NetDiversification","_",habitat,"_","Temperature",".txt"))
print( summary(cors_temp) )
print( quantile(cors_temp, c(0.025, 0.975)) )
print (sd(cors_temp) )
print(wilcox.test(cors_temp, mu=0.0))
sink()

# NetDiversification Correlation for Haq sea level
for (i in 1:numberOfSims ) {
  interpdiv = approx(times, netdiv[i,], haq$Age, method='linear', rule=1)
  end = which(is.na(interpdiv$y))
  if (length(end) == 0) {
    div_rates = interpdiv$y
    ft = haq$SL
  } else {
    div_rates = interpdiv$y[-end]
    ft = haq$SL[-end]
  }
  t = DCCA_LRCC(as.numeric(unlist(div_rates)),as.numeric(unlist(ft)),length(ft)/df)
  cors_haq[i] = t[1]
}

plot=qplot(cors_haq, geom="histogram")
plot2 = plot + geom_histogram(fill="lightgrey") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_vline(aes(xintercept=mean(cors_haq)+sd(cors_haq)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_haq)-sd(cors_haq)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_haq)), linetype="solid") + ylab("Frequency") + xlab("Correlation co-efficient")
ggsave(plot2,file=paste0("NetDiversification","_",habitat,"_","HaqSeaLevel",".pdf"))
sink(file=paste0("NetDiversification","_",habitat,"_","HaqSeaLevel",".txt"))
print( summary(cors_haq) )
print( quantile(cors_haq, c(0.025, 0.975)) )
print (sd(cors_haq) )
print(wilcox.test(cors_haq, mu=0.0))
sink()

# Net Diversification Correlation for LTT
for (i in 1:numberOfSims ) {
  interpdiv = approx(times, netdiv[i,], LTT$time, method='linear', rule=1)
  end = which(is.na(interpdiv$y))
  if (length(end) == 0) {
    div_rates = interpdiv$y
    ft = LTT$N
  } else {
    div_rates = interpdiv$y[-end]
    ft = LTT$N[-end]
  }
  t = DCCA_LRCC(as.numeric(unlist(div_rates)),as.numeric(unlist(ft)),length(ft)/df)
  cors_LTT[i] = t[1]
}

plot=qplot(cors_LTT, geom="histogram")
plot2 = plot + geom_histogram(fill="lightgrey") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_vline(aes(xintercept=mean(cors_LTT)+sd(cors_LTT)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_LTT)-sd(cors_LTT)), linetype="dashed") + geom_vline(aes(xintercept=mean(cors_LTT)), linetype="solid") + ylab("Frequency") + xlab("Correlation co-efficient")
ggsave(plot2,file=paste0("NetDiversification","_",habitat,"_","LTT",".pdf"))
sink(file=paste0("NetDiversification","_",habitat,"_","LTT",".txt"))
print( summary(cors_LTT) )
print( quantile(cors_LTT, c(0.025, 0.975)) )
print (sd(cors_LTT) )
print(wilcox.test(cors_LTT, mu=0.0))
sink()

###SAVE DATA###

save(list = ls(all.names = TRUE), file = "BAMM_Pseudosuchia_Correlations.RData", envir = .GlobalEnv)

