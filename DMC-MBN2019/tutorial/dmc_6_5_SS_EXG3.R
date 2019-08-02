##################  DMC Lesson 6: More Models

### Lesson 6.5: ExGaussian stop-signal model with 1 stop & 2 go accumulators 


### An EX-Gaussian stop-signal example with TRIGGER and GO FAILURE and 
# CONTEXT INDEPENDENT parametrization (i.e., same go parameters on go and stop trials).
# This example uses a 3-accumulator race suitable for high error rates on the go task.
# It also uses a probit scale for the tf and gf parameters, and truncated normal
# priors (this could also be done with the 2 accumulator model), which work best 
# with hierarchical models (although these are not fit here).

rm(list=ls())
# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("EXG-SS","exgSSprobit.R")
# load_data("dmc_6_5.RData")

# We will look at the most realistic scenario:
is.tf <- TRUE
is.gf <- TRUE
use.staircase <- TRUE

# Note that in cases where we set trigger failure and go failure to zero on the
# probit scale we set it to -12 (pnorm(-12)=1.776482e-33)
if (!is.tf & !is.gf) {
  
  model <- model.dmc(
    # SS stands for trial type (GO or Stop-signal [SS]):
    factors=list(S=c("s1","s2"),SS=c("GO","SS")),
    # NR stands for "No response", i.e., go omission & successful inhibitions:
    responses=c("NR","r1","r2"),
    # Match scores correct responses for each GO stimulus as usual, and also scores 
    # the "correct" stimulus corresponding to an NR response, but the latter has
    # no effect (it just avoids a standard check making sure that each response is scored): 
    match.map=list(M=list(s1="r1",s2="r2",s1="NR")),
    p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1"),
    # No trigger failures and no go failures:            
    constants=c(tf=-12,gf=-12),
    type="exgss")
  
  # This parameter setting will produce about 25% (relatively slow) errors:
  p.vector  <- c(mu.true=.5,mu.false=.60,muS=.2,
                 sigma.true=.05,sigma.false=.030,sigmaS=.03,
                 tau.true=.08,tau.false=.04,tauS=.05)
}

if (is.tf & !is.gf) {
  
  model <- model.dmc(
    # SS stands for trial type (GO or Stop-signal [SS]):
    factors=list(S=c("s1","s2"),SS=c("GO","SS")),
    # NR stands for "No response", i.e., go omission & successful inhibitions:
    responses=c("NR","r1","r2"),
    # Match scores correct responses for each GO stimulus as usual, and also scores 
    # the "correct" stimulus corresponding to an NR response, but the latter has
    # no effect (it just avoids a standard check making sure that each response is scored): 
    match.map=list(M=list(s1="r1",s2="r2",s1="NR")),
    p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1"),
    # No go failures: 
    constants=c(gf=-12),
    type="exgss")
  
  # This parameter setting will produce about 25% (relatively slow) errors:
  p.vector  <- c(mu.true=.5,mu.false=.60,muS=.2,
                 sigma.true=.05,sigma.false=.030,sigmaS=.03,
                 tau.true=.08,tau.false=.04,tauS=.05,
                 tf=qnorm(.1))
}

if (!is.tf & is.gf) {
  
  model <- model.dmc(
    # SS stands for trial type (GO or Stop-signal [SS]):
    factors=list(S=c("s1","s2"),SS=c("GO","SS")),
    # NR stands for "No response", i.e., go omission & successful inhibitions:
    responses=c("NR","r1","r2"),
    # Match scores correct responses for each GO stimulus as usual, and also scores 
    # the "correct" stimulus corresponding to an NR response, but the latter has
    # no effect (it just avoids a standard check making sure that each response is scored): 
    match.map=list(M=list(s1="r1",s2="r2",s1="NR")),
    p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1"),
    # No trigger failures: 
    constants=c(tf=-12),
    type="exgss")
  
  # This parameter setting will produce about 25% (relatively slow) errors:
  p.vector  <- c(mu.true=.5,mu.false=.60,muS=.2,
                 sigma.true=.05,sigma.false=.030,sigmaS=.03,
                 tau.true=.08,tau.false=.04,tauS=.05,
                 gf=qnorm(.1))
}

if (is.tf & is.gf) {
  
  model <- model.dmc(
    # SS stands for trial type (GO or Stop-signal [SS]):
    factors=list(S=c("s1","s2"),SS=c("GO","SS")),
    # NR stands for "No response", i.e., go omission & successful inhibitions:
    responses=c("NR","r1","r2"),
    # Match scores correct responses for each GO stimulus as usual, and also scores 
    # the "correct" stimulus corresponding to an NR response, but the latter has
    # no effect (it just avoids a standard check making sure that each response is scored): 
    match.map=list(M=list(s1="r1",s2="r2",s1="NR")),
    p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1"),
    type="exgss")
  
  # This parameter setting will produce about 25% (relatively slow) errors:
  p.vector  <- c(mu.true=.5,mu.false=.60,muS=.2,
                 sigma.true=.05,sigma.false=.030,sigmaS=.03,
                 tau.true=.08,tau.false=.04,tauS=.05,
                 tf=qnorm(.1),gf=qnorm(.1))
}

check.p.vector(p.vector,model)

# This's how (toy) stop-signal data look like with fixed and staircase SSDs;
# Now we have errors:
n <- 6
data_fixed <- data.model.dmc(simulate.dmc(p.vector,model,n=n,SSD=c(Inf,Inf,.32,.32)),model)
data_fixed

data_stair <- data.model.dmc(simulate.dmc(p.vector,model,n=n,staircase=.05, SSD=c(Inf,Inf,.25,.25)),model)
data_stair

# Make more realistic data:
n <- c(375,375,125,125)
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,n=n,SSD=c(Inf,Inf,.25,.25)),model)

# Plot go RT distribution;
# P(NA) = Probability of go omission;
# Accuracy is computed as correct go/(all go - go omissions):
correct <- as.numeric(data$S)==(as.numeric(data$R)-1)
layout(1)
plot.cell.density(data[data$SS=="GO",],
                  C=correct[data$SS=="GO"],
                  xlim=c(0,5),ymax=5,main="Go RTs")

# Overall accuracy (i.e., all trials-(errors + omission)
tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]
#  GO 
# 0.676  

# Show the different SSDs:
sort(tapply(as.character(data$SSD),data[,c("SS")],unique)$SS)
# "0.15" "0.2"  "0.25" "0.3"  "0.35" "0.4"  "0.45"

# Show the number of trials for each SSD:
Ns = tapply(data$RT,data$SSD,length)
Ns
# 0.15  0.2 0.25  0.3 0.35  0.4 0.45  Inf 
#    1   10   43   79   77   35    5  750 

# Show response rate:
tapply(!is.na(data$RT),data[,c("SS")],mean)
#   GO    SS 
# 0.896 0.500 

# Response rate broken down by SSD & corresponding inhibition function:
tapply(!is.na(data$RT),data[,c("SS","SSD")],mean)
layout(1)
plot_SS_if.dmc(data)  # P(Respond) increases as a function of SSD, as it should
Ns

# Median signal-respond RT per SSD:
tapply(data$RT,data$SSD,median,na.rm=TRUE)
layout(1)
plot_SS_srrt.dmc(data) # Median SRRT increases as a function of SSD, as it should
Nr = tapply(!is.na(data$RT),data[,c("SS","SSD")],sum)[2,]
Nr

# Signal-respond RTs should be faster than go RTs:
hist(data$RT[data$SS=="GO"],breaks="fd",main="SRRT vs. Go RT",freq=F,ylim=c(0,8))
lines(density(data$RT[data$SS=="SS"],na.rm=T),col="red",lwd=2)

###----------------------------------------- Let's start fitting

# Truncated normal priors:
p1 <- p.vector; p1[1:length(p1)] <- 1; p1[10:11] <- 3
p.prior <- prior.p.dmc(
  dists = rep("tnorm",length(p1)),p1=p.vector,p2=p1, 
  lower=c(rep(0,length(p1)-2),-12,-12),upper=c(rep(2,3),rep(.5,6),rep(12,2)) 
)
par(mfcol=c(2,6)); for (i in names(p.prior)) plot.prior(i,p.prior)

# Start sampling
samples <- samples.dmc(nmc=100,p.prior,data)
samples <- run.unstuck.dmc(samples,report=10,cores=4,p.migrate=.05,verbose=TRUE)
layout(1)
plot.dmc(samples,pll.chain=TRUE)

samples2 <- run.converge.dmc(samples.dmc(samples=samples,nmc=50,thin=5),
                             report=10,cores=4,cut=1.1,verbose=TRUE,nmc=50,minN=500,max.try=20)

# Posterior loglikelihood looks nice:
layout(1)
plot.dmc(samples2,pll.chain=TRUE,density=FALSE,smooth=FALSE)
# Parameter chains look like fat hairy catapillars:
plot.dmc(samples2,layout=c(2,6))

# Rhat looks good:
gelman.diag.dmc(samples2)
#Potential scale reduction factors:
#  
#  Point est. Upper C.I.
#mu.true           1.03       1.04
#mu.false          1.07       1.09
#sigma.true        1.01       1.02
#sigma.false       1.04       1.06
#tau.true          1.02       1.03
#tau.false         1.09       1.14
#tf                1.03       1.04
#muS               1.03       1.05
#sigmaS            1.02       1.03
#tauS              1.03       1.04
#gf                1.03       1.04
#
#Multivariate psrf
#
#1.09

# Good parameter recovery:
summary.dmc(samples2)
p.vector
#1. Empirical mean and standard deviation for each variable,
#plus standard error of the mean:
#  
#  Mean       SD  Naive SE Time-series SE
#mu.true      0.50292 0.005443 6.700e-05      0.0002148
#mu.false     0.60534 0.009991 1.230e-04      0.0004402
#sigma.true   0.04753 0.003389 4.171e-05      0.0001379
#sigma.false  0.03363 0.004245 5.225e-05      0.0001726
#tau.true     0.07751 0.007566 9.313e-05      0.0002945
#tau.false    0.03811 0.012191 1.501e-04      0.0005468
#tf          -3.45827 1.921862 2.366e-02      0.0858126
#muS          0.18987 0.024473 3.012e-04      0.0009652
#sigmaS       0.04200 0.025486 3.137e-04      0.0010060
#tauS         0.06058 0.030896 3.803e-04      0.0013116
#gf          -1.37043 0.063004 7.755e-04      0.0025367
#
#2. Quantiles for each variable:
#  
#  2.5%      25%      50%      75%    97.5%
#mu.true      0.492729  0.49911  0.50287  0.50650  0.51399
#mu.false     0.589104  0.59853  0.60421  0.61052  0.63179
#sigma.true   0.041395  0.04524  0.04743  0.04971  0.05465
#sigma.false  0.025766  0.03076  0.03349  0.03628  0.04284
#tau.true     0.063062  0.07239  0.07754  0.08255  0.09218
#tau.false    0.007963  0.03139  0.03847  0.04609  0.06055
#tf          -8.330897 -4.53135 -2.99293 -1.94615 -1.08751
#muS          0.140582  0.17440  0.18969  0.20540  0.23878
#sigmaS       0.002867  0.02143  0.03959  0.05928  0.09518
#tauS         0.005740  0.03936  0.06077  0.07881  0.12399
#gf          -1.495243 -1.41347 -1.36885 -1.32663 -1.25073

# Check parameter recovery:
check.recovery.dmc(samples2,p.vector)
#mu.true mu.false sigma.true sigma.false tau.true tau.false    tf   muS sigmaS tauS    gf
#True              0.50     0.60       0.05        0.03     0.08      0.04 -1.28  0.20   0.03 0.05 -1.28
#2.5% Estimate     0.49     0.59       0.04        0.03     0.06      0.01 -8.33  0.14   0.00 0.01 -1.50
#50% Estimate      0.50     0.60       0.05        0.03     0.08      0.04 -2.99  0.19   0.04 0.06 -1.37
#97.5% Estimate    0.51     0.63       0.05        0.04     0.09      0.06 -1.09  0.24   0.10 0.12 -1.25
#Median-True       0.00     0.00       0.00        0.00     0.00      0.00 -1.71 -0.01   0.01 0.01 -0.09
#attr(,"ci50")
#
#25%         75%
#  mu.true      0.49911311  0.50649508
#mu.false     0.59852586  0.61051502
#sigma.true   0.04523772  0.04971335
#sigma.false  0.03075932  0.03627898
#tau.true     0.07238736  0.08254630
#tau.false    0.03138788  0.04609305
#tf          -4.53134586 -1.94614847
#muS          0.17439977  0.20540290
#sigmaS       0.02143479  0.05928286
#tauS         0.03936175  0.07881407
#gf          -1.41346591 -1.32663429

# Plot the posteriors:
# Note that for given n, EXG3 typically results in wider posteriors than EXG2, 
# because EXG3 estimates three extra parameters:
plot.dmc(samples2,layout=c(2,6),p.prior=p.prior)

# Good fits, as expected for simulated data: 
pp <- post.predict.dmc(samples2,n.post=1000)
pp2 <- post.predict.dmc(samples2,n.post=1000,save.simulation=TRUE)

# Compare cdf of observed and predicted go RTs:
plot.pp.dmc(pp,style="cdf",layout=c(2,2),ylim=c(0,1),fits.pcol="grey",
            model.legend=FALSE,dname="",aname="")

# Compare observed (black dots) and predicted (gray violins) response rates on the different SSDs
# using posterior-predictive simulations. 
layout(1)
plot_SS_if.dmc(data=samples2$data,sim=pp2)
#Inhibition function 
#SSD  n     p
#1 0.15  1    NA
#2 0.20 10 0.745
#3 0.25 40 0.623
#4 0.30 73 0.512
#5 0.35 73 0.419
#6 0.40 40 0.716
#7 0.45 12 0.442
#8 0.50  1    NA

# Compare observed (black dots) and predicted (gray violins) median signal-respond RTs using
# posterior-predictive simulations. 
layout(1)
plot_SS_srrt.dmc(data=samples2$data,sim=pp2)
#Median signal-respond RTs 
#SSD nrt     p n.sim
#1 0.20   1 0.607   745
#2 0.25   9 0.517  1000
#3 0.30  29 0.721  1000
#4 0.35  44 0.664  1000
#5 0.40  29 0.483  1000
#6 0.45  11 0.369  1000
#7 0.50   1 0.286   872

### Extension exercise 1 ----

# Compare the observed vs. predicted error rates on Go trials.
pp.go <- subset(pp2,SS=="GO")
data.go <- subset(samples2$data,SS=="GO")

# Predicted:
pp.error <- with(pp.go,tapply((as.numeric(S)!=(as.numeric(R)-1) & !is.na(RT)),reps,function(x) mean(x)))
# Observed:
error <- mean(as.numeric(data.go$S)!=(as.numeric(data.go$R)-1) & !is.na(data.go$RT))
# Plot
layout(1)
hist(pp.error)
points(error,1,col="red",pch=16,cex=2)
# Posterior predictive p-value:
mean(pp.error>=error)

### Extension exercise 2 ----

# Fit the data from sample2 with a misspecified model that does not account for errors (see fisrt stop-signal lesson).
# How does this effect the parameter estimates?

# Specify the model without errors:
modelMis <- model.dmc(
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),
  responses=c("NR","r1","r2"),
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")),
  p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1"),
  # No errors:
  constants=c(mu.false=1e6,sigma.false=.001,tau.false=.001),
  type="exgss")

p.vectorMis  <- c(mu.true=.5,muS=.2,sigma.true=.05,sigmaS=.03,tau.true=.08,tauS=.05,tf=qnorm(.1),gf=qnorm(.1))

# Remove error RTs:
remove <- as.numeric(samples2$data$S)!=(as.numeric(samples2$data$R)-1) & !is.na(samples2$data$RT)
dataMis <- data.model.dmc(samples2$data[!remove,],modelMis)

# Sanity checks;
# Plot go RT distribution with accuracy=1
correctMis <- as.numeric(dataMis$S)==(as.numeric(dataMis$R)-1)
layout(1)
plot.cell.density(dataMis[dataMis$SS=="GO",],
                  C=correctMis[dataMis$SS=="GO"],
                  ymax=6,main="Go RTs")

# Overall accuracy on go task (i.e., proportion go omissions in this case; remember, no errors!):
tapply(as.numeric(dataMis$S)==(as.numeric(dataMis$R)-1),dataMis$SS,mean,na.rm=TRUE)["GO"]

# Show response rate;
# Go omission rate changes relative to samples2$data simply because number of Go trials changes;
# Stop response rate decreases because we remove only the signal-respond RTs:
tapply(!is.na(dataMis$RT),dataMis[,c("SS")],mean)
tapply(!is.na(samples2$data$RT),samples2$data[,c("SS")],mean)

# As before, truncated normal priors, but now we don't need prior for mu.false, sigma.false, and tau.false:
p1Mis <- p.vectorMis; p1Mis[1:length(p1Mis)] <- 1; p1Mis[7:8] <- 3
p.priorMis <- prior.p.dmc(
  dists = rep("tnorm",length(p1Mis)),p1=p.vectorMis,p2=p1Mis, 
  lower=c(rep(0,length(p1Mis)-2),-12,-12),upper=c(rep(2,2),rep(.5,4),rep(12,2)) 
)
par(mfcol=c(2,6)); for (i in names(p.priorMis)) plot.prior(i,p.priorMis)

# Start sampling:
samples3 <- samples.dmc(nmc=100,p.priorMis,dataMis)
samples3 <- run.unstuck.dmc(samples3,report=10,cores=4,p.migrate=.05,verbose=TRUE)
layout(1)
plot.dmc(samples3,pll.chain=TRUE,density=FALSE,smooth=FALSE)

samples4 <- run.converge.dmc(samples.dmc(samples=samples3,nmc=50,thin=5),report=10,cores=4,cut=1.1,verbose=TRUE,nmc=50,minN=500,max.try=20)

# Posterior loglikelihood looks nice:
plot.dmc(samples4,pll.chain=TRUE,density=FALSE,smooth=FALSE)
# Parameter chains look like fat hairy catapillars:
plot.dmc(samples4,layout=c(2,4),smooth=FALSE,density=FALSE)

# Rhat looks good:
gelman.diag.dmc(samples4)

# Biased estimates:
summary.dmc(samples4)
p.vector

# Priors are updated but posteriors are off:
plot.dmc(samples4,layout=c(2,4),p.prior=p.priorMis)
p.vector

# Compare to true model:
# The asymptotic underestimation of tauS is less obvious here beacuse of the relatively
# wide posteriors
# Note that gf changes relative to samples2 simply because number of Go trials changes
par(mfcol=c(2,4))
for(i in names(p.vectorMis)){
  hist(samples2$theta[,i,],main=i,freq=F,breaks="fd")
  lines(density(samples4$theta[,i,]),col="red",lwd=2)
}

ppMis <- post.predict.dmc(samples4,n.post=1000)
ppMis2 <- post.predict.dmc(samples4,n.post=1000,save.simulation=TRUE)

# Good fit to observed correct go RTs:
layout(1)
plot.pp.dmc(ppMis,style="cdf",layout=c(2,2),ylim=c(0,1),fits.pcol="grey",
            model.legend=FALSE,dname="",aname="")

# Misspecified model often produces more variable predictions:
# and more extreme p values (but not always...)
par(mfcol=c(1,2))
plot_SS_if.dmc(data=samples4$data,sim=ppMis2)
plot_SS_if.dmc(data=samples2$data,sim=pp2)

par(mfcol=c(1,2))
plot_SS_srrt.dmc(data=samples4$data,sim=ppMis2,ylim=c(.3,.8))
plot_SS_srrt.dmc(data=samples2$data,sim=pp2,ylim=c(.3,.8))

# save_data(data,samples,samples2,pp,pp2,dataMis,samples3,samples4,ppMis,ppMis2,file="dmc_6_5.RData")

#### Extension Exercise 3: Parameter-recovery ----

# This extension has its own data file so lets start again.
rm(list=ls())
source ("dmc/dmc.R")
load_model ("EXG-SS","exgSS.R")
# load_data("dmc_6_EXGRecovery.RData")

# Specify model:
model <- model.dmc(
    factors=list(S=c("s1","s2"),SS=c("GO","SS")),
    responses=c("NR","r1","r2"),
    match.map=list(M=list(s1="r1",s2="r2",s1="NR")),
    p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1"),
    type="exgss")
  
# This parameter setting will produce about 20% (relatively slow) errors
p.vector  <- c(mu.true=.45,mu.false=.572,muS=.2,sigma.true=.05,sigma.false=.019,sigmaS=.03,
               tau.true=.08,tau.false=.03,tauS=.05,tf=.1,gf=.1)
check.p.vector(p.vector,model)

# Specify the number of go and stop-signal trials
# Here we use 600 go trials (300 for each stimuli) and 200 stop-signal trials
n.trials <- c(300,300,100,100)

# Specify the number of replicated data sets
# This is just a toy example, in a real recovery study n.rep should me *much* higher
# We recommend at least 100 replications
n.rep <- 10

# First simulate n.rep replications, each with identical parameters but with
# different random realizations of data (see tutorial dmc_4_1 for more information) 
# Here we use staircase-tracking with step size .05s starting at .25s.
raw.data <- h.simulate.dmc(model,ps=p.vector,ns=n.rep,n=matrix(rep(n.trials,each=n.rep),ncol=4),
                       SSD = matrix(rep(c(Inf,Inf,.25,.25),each=n.rep),ncol=4),staircase=0.05)
                                                           
###------ Some checks
# Check generating parameter values
attributes(raw.data)$parameters
# Check number of trials
table(raw.data[,c("S","SS","s")])
# Check number of stop-signal trials per SSD
table(raw.data[,c("SS","SSD","s")])
# Check go omission rates on go trials and inhibition rates on stop-signal trials
tapply(is.na(raw.data$RT),raw.data[,c("s","SS")],mean)
# Plot Go RT distribution collapsed across n.rep
# Also show accuracy and go omission rate (P(NA))
correct <- as.numeric(raw.data$S)==(as.numeric(raw.data$R)-1)
layout(1)
plot.cell.density(raw.data[raw.data$SS=="GO",],
                   C=correct[raw.data$SS=="GO"],
                   xlim=c(0,3000),main="Go RTs")

# The data-model object is now a list, with one replication in each slot
data <- data.model.dmc(raw.data,model)

# Set up priors; here we will use uniform priors
# (see dmc_2_1 for more information)
p.prior <- prior.p.dmc(
  dists = rep("beta",length(p.vector)),
  p1= c(mu.true=1,mu.false=1,muS=1,sigma.true=1,sigma.false=1,sigmaS=1,tau.true=1,tau.false=1,tauS=1,tf=1,gf=1),
  p2= rep(1,length(p.vector)),
  lower=rep(0,length(p.vector)),upper=c(rep(2,length(p.vector)-2),rep(1,2)))

# Plot priors
par(mfcol=c(2,6)); for (i in names(p.prior)) plot.prior(i,p.prior,ylim = c(0,4))

# Set up object ready for sampling and generate start values for each replication 
# For more information, see dmc_4_4
samples <- h.samples.dmc(nmc=120,p.prior,data,thin=10)

# Run the sampling for each replication 
# (see DMCpaper1.R for more information) 
# By default, we use one core per replication, which is more efficient than spreading
# chains across cores
samples2 <- h.RUN.dmc(samples,cores=10,meanN=2000,use.effectiveSize=FALSE)

### Check recovery performance for each replication; for instance, for n.rep=2

# Do the MCMC chains look like fat hairy caterpillars?
plot.dmc(samples2[[2]],layout=c(2,6))
# Does Rhat look good?
gelman.diag.dmc(samples2[[2]])
# Do the posterior distributions look OK (e.g., not truncated by arbitrary prior bounds)?
plot.dmc(samples2,layout=c(2,6),p.prior=p.prior)
# Do the parameter correlations look reasonable (e.g., no perfect correlations)?
# Note that in the ex-Gaussian distribution,
# some of the parameter correlations are quite strong,
# typically between mu and tau parameters of the same type
pairs.dmc(samples2[[2]],thin=25)

# Check parameter recovery
# The first table shows for each parameter:
# the true data-generating value, 
# the 2.5th, 50th, and 97.5th percentile of the posterior distribution,
# and the difference between the posterior median and the true value
# The second table shows for each parameter 
# the 25th and the 50th percentile of the posterior distribution
check.recovery.dmc(samples2[[2]],p.vector)

### Check average recovery performance

# Does Rhats (multivariate psrf) look good?
gelman.diag.dmc(samples2)
# Plot parameter correlations
pairs.dmc(samples2,thin=25)

# Check parameter recovery
# h.check.recovery applies check.recovery to each fit then takes the average. It
# also calculates the "coverage" of the 50% and 95% credible interval, i.e., the percentage 
# of repetitions for which the true value falls in the 50% and 95% credible interval, respectively. 
# If uncertainty is being properly quantified these should be around 50% and 95%.
# Note this can be a little slow to run with high n.rep! 
h.check.recovery.dmc(samples2,p.vector,do.coverage=TRUE)

# save_data(samples2,file="dmc_6_EXGRecovery.RData")
