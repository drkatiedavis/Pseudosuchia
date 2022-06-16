library(BAMMtools)
library(adephylo)
library(pracma)
library(caper)
library(phangorn)
library(phytools)
library(depmixS4)
library(RTransferEntropy)
library(ggplot2)

# set files

tree <- read.tree("FBD_Pseudosuchia_041219.tre")
bamm_data <- "PseudosuchiaNov2020_event_data.txt"
numberOfSims <- 100

#set variables for correlation

delo18 <- read.csv("science_delO18_curve.csv", header=FALSE)
haq <- read.csv("haq_sealevel_curve.csv",header=TRUE)
LTT <- read.csv('coords.csv', header=TRUE)

# Sort out temperature curves
# remove some autocorrelation

finaltemp = smooth(smooth(delo18$V2))

# make the correlation to temperature (inverse of DelO18)
finaltemp = -finaltemp

# sort out your tree and load in BAMM data
edata <- getEventData(tree, eventdata = bamm_data, burnin=0.1)
lengths = distRoot(tree, tree$tip.label)
root = max(lengths)
count = 1

# load in habitat data
habitatdata <- read.csv("Pseudosuchia_Habitat.txt", header=F, stringsAsFactors = FALSE)

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

set.seed(123)


##################################################
#################TRANSFER ENTROPY#################
##################################################


# new interpolation function to simplify code
interpOntoVariable<-function(times1, ts1, times2, ts2) {
    interpTS1 <- approx(times1, ts1, times2, method='linear', rule=1)
    end <- which(is.na(interpTS1$y))
    if (length(end) == 0) {
        finalTS1 <- interpTS1$y
        finalTS2 <- ts2
    } else {
        finalTS1 <- interpTS1$y[-end]
        finalTS2 <- ts2[-end]
    }
    return(data.frame("ts1" = finalTS1, "ts2" = finalTS2))
}

# function to find correct number of bins
findNoBins<-function(ts) {
    for (i in 3:(0.1*length(ts))) {
        bins <- seq(min(ts),max(ts),length.out=i)
        h<-hist(ts, breaks=bins, plot=F)
        if (min(h$counts)==0){
            break
        }
    }
    # want the last i, then bins is breaks-1 
    return(i-2)
}

# function to find correct number of states
findNstates<-function(ts) {
    minAIC = 1000000.0
    dat <- data.frame(as.numeric(ts))
    for (i in 2:20) {
        mod = depmix(ts ~ 1, nstates = i, data=dat)
        fitmod = fit(mod)
        aic_val <- AIC(fitmod)
        if (aic_val < minAIC) {
            minAIC = aic_val
        } else {
            break
        }
    }
    return(i-1)
}

###SPECIATION###

#LTT

transfer_entropy_values  <- rep(NA, numberOfSims)
significance <- rep(NA, numberOfSims)
for (i in 1:numberOfSims ) {
    # interpolate onto same times
    timeseries <- interpOntoVariable(times, rtt$lambda[i,], LTT$time, LTT$N)
    lambda <- timeseries$ts1
    variableOI <- timeseries$ts2
    # now find bins
    bins <- findNoBins(variableOI)
    # now find the HMM states; max of 20 states
    ##########################################
    # these are wrapped in try blocks so errors are ignored
    #
    # The output might there contain NAs, so account for that later
    ###########################################
    tryCatch( {
        states <- findNstates(variableOI);
        # finally, you can do the TE and store
        results <-transfer_entropy(variableOI, lambda, lx=states, ly=states, type=c('bins'), bins=bins);
        # check the $ variables here - I don't know what the output contains
        transfer_entropy_values[i] <- results$coef[1];
        significance[i] <- results$coef[7];
        }, error = function(e) {significance[i] <- 2.0 }) # so we ignore this output
}

# the 0.05 is the p-value of interest
sigp <- 0.05
indices <- which(significance < sigp)
sig_te <- transfer_entropy_values[indices]

plot=qplot(sig_te, geom="histogram")
plot2 = plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(plot2,file="transfer_entropy_values_VARIABLE_LTT_Speciation.pdf")
sink(file="TE_VARIABLE_LTT_Speciation.txt")
print( paste0(length(sig_te), " significant values to p<", sigp)) 
print( summary(sig_te) )
print (sd(sig_te) )
#print(wilcox.test(sig_te, mu=0.0))
sink()

#Temperature

transfer_entropy_values  <- rep(NA, numberOfSims)
significance <- rep(NA, numberOfSims)
for (i in 1:numberOfSims ) {
    # interpolate onto same times
    timeseries <- interpOntoVariable(times, rtt$lambda[i,], delo18$V1, finaltemp)
    lambda <- timeseries$ts1
    variableOI <- timeseries$ts2
    # now find bins
    bins <- findNoBins(variableOI)
    # now find the HMM states; max of 20 states
    ##########################################
    # these are wrapped in try blocks so errors are ignored
    #
    # The output might there contain NAs, so account for that later
    ###########################################
    tryCatch( {
        states <- findNstates(variableOI);
        # finally, you can do the TE and store
        results <-transfer_entropy(variableOI, lambda, lx=states, ly=states, type=c('bins'), bins=bins);
        # check the $ variables here - I don't know what the output contains
        transfer_entropy_values[i] <- results$coef[1];
        significance[i] <- results$coef[7];
        }, error = function(e) {significance[i] <- 2.0 }) # so we ignore this output
}

# the 0.05 is the p-value of interest
sigp <- 0.05
indices <- which(significance < sigp)
sig_te <- transfer_entropy_values[indices]

plot=qplot(sig_te, geom="histogram")
plot2 = plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(plot2,file="transfer_entropy_values_VARIABLE_Temp_Speciation.pdf")
sink(file="TE_VARIABLE_Temp_Speciation.txt")
print( paste0(length(sig_te), " significant values to p<", sigp)) 
print( summary(sig_te) )
print (sd(sig_te) )
#print(wilcox.test(sig_te, mu=0.0))
sink()

#Sea level

transfer_entropy_values  <- rep(NA, numberOfSims)
significance <- rep(NA, numberOfSims)
for (i in 1:numberOfSims ) {
    # interpolate onto same times
    timeseries <- interpOntoVariable(times, rtt$lambda[i,], haq$Age, haq$SL)
    lambda <- timeseries$ts1
    variableOI <- timeseries$ts2
    # now find bins
    bins <- findNoBins(variableOI)
    # now find the HMM states; max of 20 states
    ##########################################
    # these are wrapped in try blocks so errors are ignored
    #
    # The output might there contain NAs, so account for that later
    ###########################################
    tryCatch( {
        states <- findNstates(variableOI);
        # finally, you can do the TE and store
        results <-transfer_entropy(variableOI, lambda, lx=states, ly=states, type=c('bins'), bins=bins);
        # check the $ variables here - I don't know what the output contains
        transfer_entropy_values[i] <- results$coef[1];
        significance[i] <- results$coef[7];
        }, error = function(e) {significance[i] <- 2.0 }) # so we ignore this output
}

# the 0.05 is the p-value of interest
sigp <- 0.05
indices <- which(significance < sigp)
sig_te <- transfer_entropy_values[indices]

plot=qplot(sig_te, geom="histogram")
plot2 = plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(plot2,file="transfer_entropy_values_VARIABLE_Haq_Speciation.pdf")
sink(file="TE_VARIABLE_Haq_Speciation.txt")
print( paste0(length(sig_te), " significant values to p<", sigp)) 
print( summary(sig_te) )
print (sd(sig_te) )
#print(wilcox.test(sig_te, mu=0.0))
sink()

###EXTINCTION###

#LTT

transfer_entropy_values  <- rep(NA, numberOfSims)
significance <- rep(NA, numberOfSims)
for (i in 1:numberOfSims ) {
    # interpolate onto same times
    timeseries <- interpOntoVariable(times, rtt$mu[i,], LTT$time, LTT$N)
    mu <- timeseries$ts1
    variableOI <- timeseries$ts2
    # now find bins
    bins <- findNoBins(variableOI)
    # now find the HMM states; max of 20 states
    ##########################################
    # these are wrapped in try blocks so errors are ignored
    #
    # The output might there contain NAs, so account for that later
    ###########################################
    tryCatch( {
        states <- findNstates(variableOI);
        # finally, you can do the TE and store
        results <-transfer_entropy(variableOI, mu, lx=states, ly=states, type=c('bins'), bins=bins);
        # check the $ variables here - I don't know what the output contains
        transfer_entropy_values[i] <- results$coef[1];
        significance[i] <- results$coef[7];
        }, error = function(e) {significance[i] <- 2.0 }) # so we ignore this output
}

# the 0.05 is the p-value of interest
sigp <- 0.05
indices <- which(significance < sigp)
sig_te <- transfer_entropy_values[indices]

plot=qplot(sig_te, geom="histogram")
plot2 = plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(plot2,file="transfer_entropy_values_VARIABLE_LTT_Extinction.pdf")
sink(file="TE_VARIABLE_LTT_Extinction.txt")
print( paste0(length(sig_te), " significant values to p<", sigp)) 
print( summary(sig_te) )
print (sd(sig_te) )
#print(wilcox.test(sig_te, mu=0.0))
sink()

#Temperature

transfer_entropy_values  <- rep(NA, numberOfSims)
significance <- rep(NA, numberOfSims)
for (i in 1:numberOfSims ) {
    # interpolate onto same times
    timeseries <- interpOntoVariable(times, rtt$mu[i,], delo18$V1, finaltemp)
    mu <- timeseries$ts1
    variableOI <- timeseries$ts2
    # now find bins
    bins <- findNoBins(variableOI)
    # now find the HMM states; max of 20 states
    ##########################################
    # these are wrapped in try blocks so errors are ignored
    #
    # The output might there contain NAs, so account for that later
    ###########################################
    tryCatch( {
        states <- findNstates(variableOI);
        # finally, you can do the TE and store
        results <-transfer_entropy(variableOI, mu, lx=states, ly=states, type=c('bins'), bins=bins);
        # check the $ variables here - I don't know what the output contains
        transfer_entropy_values[i] <- results$coef[1];
        significance[i] <- results$coef[7];
        }, error = function(e) {significance[i] <- 2.0 }) # so we ignore this output
}

# the 0.05 is the p-value of interest
sigp <- 0.05
indices <- which(significance < sigp)
sig_te <- transfer_entropy_values[indices]

plot=qplot(sig_te, geom="histogram")
plot2 = plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(plot2,file="transfer_entropy_values_VARIABLE_Temp_Extinction.pdf")
sink(file="TE_VARIABLE_Temp_Extinction.txt")
print( paste0(length(sig_te), " significant values to p<", sigp)) 
print( summary(sig_te) )
print (sd(sig_te) )
#print(wilcox.test(sig_te, mu=0.0))
sink()

#Sea level

transfer_entropy_values  <- rep(NA, numberOfSims)
significance <- rep(NA, numberOfSims)
for (i in 1:numberOfSims ) {
    # interpolate onto same times
    timeseries <- interpOntoVariable(times, rtt$mu[i,], haq$Age, haq$SL)
    mu <- timeseries$ts1
    variableOI <- timeseries$ts2
    # now find bins
    bins <- findNoBins(variableOI)
    # now find the HMM states; max of 20 states
    ##########################################
    # these are wrapped in try blocks so errors are ignored
    #
    # The output might there contain NAs, so account for that later
    ###########################################
    tryCatch( {
        states <- findNstates(variableOI);
        # finally, you can do the TE and store
        results <-transfer_entropy(variableOI, mu, lx=states, ly=states, type=c('bins'), bins=bins);
        # check the $ variables here - I don't know what the output contains
        transfer_entropy_values[i] <- results$coef[1];
        significance[i] <- results$coef[7];
        }, error = function(e) {significance[i] <- 2.0 }) # so we ignore this output
}

# the 0.05 is the p-value of interest
sigp <- 0.05
indices <- which(significance < sigp)
sig_te <- transfer_entropy_values[indices]

plot=qplot(sig_te, geom="histogram")
plot2 = plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(plot2,file="transfer_entropy_values_VARIABLE_Haq_Extinction.pdf")
sink(file="TE_VARIABLE_Haq_Extinction.txt")
print( paste0(length(sig_te), " significant values to p<", sigp)) 
print( summary(sig_te) )
print (sd(sig_te) )
#print(wilcox.test(sig_te, mu=0.0))
sink()

###NET DIVERSIFICATION###

netdiv<-(rtt$lambda)-(rtt$mu)

#LTT

transfer_entropy_values  <- rep(NA, numberOfSims)
significance <- rep(NA, numberOfSims)
for (i in 1:numberOfSims ) {
    # interpolate onto same times
    timeseries <- interpOntoVariable(times, netdiv[i,], LTT$time, LTT$N)
    netdiv <- timeseries$ts1
    variableOI <- timeseries$ts2
    # now find bins
    bins <- findNoBins(variableOI)
    # now find the HMM states; max of 20 states
    ##########################################
    # these are wrapped in try blocks so errors are ignored
    #
    # The output might there contain NAs, so account for that later
    ###########################################
    tryCatch( {
        states <- findNstates(variableOI);
        # finally, you can do the TE and store
        results <-transfer_entropy(variableOI, netdiv, lx=states, ly=states, type=c('bins'), bins=bins);
        # check the $ variables here - I don't know what the output contains
        transfer_entropy_values[i] <- results$coef[1];
        significance[i] <- results$coef[7];
        }, error = function(e) {significance[i] <- 2.0 }) # so we ignore this output
}

# the 0.05 is the p-value of interest
sigp <- 0.05
indices <- which(significance < sigp)
sig_te <- transfer_entropy_values[indices]

plot=qplot(sig_te, geom="histogram")
plot2 = plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(plot2,file="transfer_entropy_values_VARIABLE_LTT_Extinction.pdf")
sink(file="TE_VARIABLE_LTT_Extinction.txt")
print( paste0(length(sig_te), " significant values to p<", sigp)) 
print( summary(sig_te) )
print (sd(sig_te) )
#print(wilcox.test(sig_te, mu=0.0))
sink()

#Temperature

transfer_entropy_values  <- rep(NA, numberOfSims)
significance <- rep(NA, numberOfSims)
for (i in 1:numberOfSims ) {
    # interpolate onto same times
    timeseries <- interpOntoVariable(times, netdiv[i,], delo18$V1, finaltemp)
    netdiv <- timeseries$ts1
    variableOI <- timeseries$ts2
    # now find bins
    bins <- findNoBins(variableOI)
    # now find the HMM states; max of 20 states
    ##########################################
    # these are wrapped in try blocks so errors are ignored
    #
    # The output might there contain NAs, so account for that later
    ###########################################
    tryCatch( {
        states <- findNstates(variableOI);
        # finally, you can do the TE and store
        results <-transfer_entropy(variableOI, netdiv, lx=states, ly=states, type=c('bins'), bins=bins);
        # check the $ variables here - I don't know what the output contains
        transfer_entropy_values[i] <- results$coef[1];
        significance[i] <- results$coef[7];
        }, error = function(e) {significance[i] <- 2.0 }) # so we ignore this output
}

# the 0.05 is the p-value of interest
sigp <- 0.05
indices <- which(significance < sigp)
sig_te <- transfer_entropy_values[indices]

plot=qplot(sig_te, geom="histogram")
plot2 = plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(plot2,file="transfer_entropy_values_VARIABLE_Temp_Extinction.pdf")
sink(file="TE_VARIABLE_Temp_Extinction.txt")
print( paste0(length(sig_te), " significant values to p<", sigp)) 
print( summary(sig_te) )
print (sd(sig_te) )
#print(wilcox.test(sig_te, mu=0.0))
sink()

#Sea level

transfer_entropy_values  <- rep(NA, numberOfSims)
significance <- rep(NA, numberOfSims)
for (i in 1:numberOfSims ) {
    # interpolate onto same times
    timeseries <- interpOntoVariable(times, netdiv[i,], haq$Age, haq$SL)
    netdiv <- timeseries$ts1
    variableOI <- timeseries$ts2
    # now find bins
    bins <- findNoBins(variableOI)
    # now find the HMM states; max of 20 states
    ##########################################
    # these are wrapped in try blocks so errors are ignored
    #
    # The output might there contain NAs, so account for that later
    ###########################################
    tryCatch( {
        states <- findNstates(variableOI);
        # finally, you can do the TE and store
        results <-transfer_entropy(variableOI, netdiv, lx=states, ly=states, type=c('bins'), bins=bins);
        # check the $ variables here - I don't know what the output contains
        transfer_entropy_values[i] <- results$coef[1];
        significance[i] <- results$coef[7];
        }, error = function(e) {significance[i] <- 2.0 }) # so we ignore this output
}

# the 0.05 is the p-value of interest
sigp <- 0.05
indices <- which(significance < sigp)
sig_te <- transfer_entropy_values[indices]

plot=qplot(sig_te, geom="histogram")
plot2 = plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(plot2,file="transfer_entropy_values_VARIABLE_Haq_Extinction.pdf")
sink(file="TE_VARIABLE_Haq_Extinction.txt")
print( paste0(length(sig_te), " significant values to p<", sigp)) 
print( summary(sig_te) )
print (sd(sig_te) )
#print(wilcox.test(sig_te, mu=0.0))
sink()

