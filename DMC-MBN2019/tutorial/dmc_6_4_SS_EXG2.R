##################  DMC Lesson 6: More Models

### Lesson 6.4: ExGaussian stop-signal model with 1 stop & 1 go accumulator (BEESTS)

# An EXG stop-signal example with TRIGGER and GO FAILURE and CONTEXT 
# INDEPENDENT parametrization (i.e., same go parameters on go and stop trials).
# 2-accumulator race suitable for high accuracy data like BEESTS. 

rm(list=ls())
# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("EXG-SS","exgSS.R")

# load_data("dmc_6_4.Rdata")

# We will look at the most realistic scenario:
is.tf <- TRUE
is.gf <- TRUE
use.staircase <- TRUE

if (!is.tf & !is.gf) {
  
  model <- model.dmc(
    # SS stands for trial type (GO or Stop-signal [SS]):
    factors=list(S=c("s1","s2"),SS=c("GO","SS")), 
    # NR stands for "No response", i.e., go omission & successful inhibitions:
    responses=c("NR","r1","r2"), 
    # Match scores correct responses for each GO stimulus as usual, and scores 
    # the "correct" stimulus corresponding to an NR response, but the latter has
    # no effect (it just avoids a standard check making sure that each response is scored): 
    match.map=list(M=list(s1="r1",s2="r2",s1="NR")), 
    p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1"),
    # No errors, no trigger failures, and no go failures: 
    constants=c(mu.false=1e6,sigma.false=.001,tau.false=.001,tf=0,gf=0),
    type="exgss")
  
  # This gives mean GoRT of .5 + .08 = .58 and SD Go RT of sqrt(.05^2+.08^2) = 0.09 
  # and mean SSRT of .2+.05 = 0.25 and SD SSRT of sqrt(.03^2+.05^2) = 0.06:
  p.vector  <- c(mu.true=.5,muS=.2,sigma.true=.05,sigmaS=.03,tau.true=.08,tauS=.05) 
}

if (is.tf & !is.gf) {
  
  model <- model.dmc(
    # SS stands for trial type (GO or Stop-signal [SS]):
    factors=list(S=c("s1","s2"),SS=c("GO","SS")),
    # NR stands for "No response", i.e., go omission & successful inhibitions:
    responses=c("NR","r1","r2"),
    # Match scores correct responses for each GO stimulus as usual, and scores 
    # the "correct" stimulus corresponding to an NR response, but the latter has
    # no effect (it just avoids a standard check making sure that each response is scored): 
    match.map=list(M=list(s1="r1",s2="r2",s1="NR")),
    p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1"),
    # No errors and no go failures:
    constants=c(mu.false=1e6,sigma.false=.001,tau.false=.001,gf=0),
    type="exgss")
  
  # This gives mean GoRT of .5 + .08 = .58 and SD Go RT of sqrt(.05^2+.08^2) = 0.09 
  # and mean SSRT of .2+.05 = 0.25 and SD SSRT of sqrt(.03^2+.05^2) = 0.06:
  p.vector  <- c(mu.true=.5,muS=.2,sigma.true=.05,sigmaS=.03,tau.true=.08,tauS=.05,tf=.1) 
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
    # No errors and no trigger failures:
    constants=c(mu.false=1e6,sigma.false=.001,tau.false=.001,tf=0),
    type="exgss")
  
  # This gives mean GoRT of .5 + .08 = .58 and SD Go RT of sqrt(.05^2+.08^2) = 0.09 
  # and mean SSRT of .2+.05 = 0.25 and SD SSRT of sqrt(.03^2+.05^2) = 0.06:
  p.vector  <- c(mu.true=.5,muS=.2,sigma.true=.05,sigmaS=.03,tau.true=.08,tauS=.05,gf=.1) 
}

if (is.tf & is.gf) {
  
  model <- model.dmc(
    # SS stands for trial type (GO or Stop-signal [SS]):
    factors=list(S=c("s1","s2"),SS=c("GO","SS")),
    # NR stands for "No response", i.e., go omission & successful inhibitions:
    responses=c("NR","r1","r2"),
    # Match scores correct responses for each GO stimulus as usual, and scores 
    # the "correct" stimulus corresponding to an NR response, but the latter has
    # no effect (it just avoids a standard check making sure that each response is scored): 
    match.map=list(M=list(s1="r1",s2="r2",s1="NR")),
    p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1"),
    # No errors:
    constants=c(mu.false=1e6,sigma.false=.001,tau.false=.001),
    type="exgss")
  
  # This gives mean GoRT of .5 + .08 = .58 and SD Go RT of sqrt(.05^2+.08^2) = 0.09 
  # and mean SSRT of .2+.05 = 0.25 and SD SSRT of sqrt(.03^2+.05^2) = 0.06:
  p.vector  <- c(mu.true=.5,muS=.2,sigma.true=.05,sigmaS=.03,tau.true=.08,tauS=.05,tf=.1,gf=.1) 
}

check.p.vector(p.vector,model)

# This's how (toy) stop-signal data look like with fixed SSDs:
n <- 6
# SSD must be a scalar, or a vector the same length as the number of cells, 
# or the same length as the data and have Inf in all go cells.
# SSD = .32 will result in response rate of about 50%, but only for THIS particular parameter setting:
data_fixed <- data.model.dmc(simulate.dmc(p.vector,model,n=n,SSD=c(Inf,Inf,.32,.32)),model)
data_fixed

# This's how (toy) stop-signal data look like with staircase tracking:
# SSD gives start of the tracking algorithm, then moves stop-signal delay 
# up und down by .05ms contingent on performance:
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
                  ymax=5,main="Go RTs")

# Overall accuracy on go task (i.e., proportion go omissions in this case; remember, no errors!):
tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]
# GO 
# 0.884

# Show the different SSDs:
sort(tapply(as.character(data$SSD),data[,c("SS")],unique)$SS)
# "0.1"  "0.15" "0.2"  "0.25" "0.3"  "0.35" "0.4"  "0.45" "0.5"  

# Show the number of trials for each SSD:
Ns = tapply(data$RT,data$SSD,length)
Ns
# 0.1 0.15  0.2 0.25  0.3 0.35  0.4 0.45  0.5  Inf 
#   1    6   18   46   76   67   28    7    1  750

# Show response rate:
tapply(!is.na(data$RT),data[,c("SS")],mean)
#    GO    SS 
# 0.884 0.492  

# Response rate broken down by SSD & corresponding inhibition function:
tapply(!is.na(data$RT),data[,c("SS","SSD")],mean)
plot_SS_if.dmc(data)  #P(Respond) increases as a function of SSD, as it should

# Plot median signal-respond RT per SSD:
tapply(data$RT,data$SSD,median,na.rm=TRUE) 
plot_SS_srrt.dmc(data) # Median SRRT increases as a function of SSD, as it should

# Show number of signal-respond RTs per SSD:
Nr = tapply(!is.na(data$RT),data[,c("SS","SSD")],sum)[2,]
Nr
# 0.1 0.15  0.2 0.25  0.3 0.35  0.4 0.45  0.5  Inf 
#  0    1    5   13   31   44   22    6    1   NA 

# Signal-respond RTs should be faster than go RTs:
hist(data$RT[data$SS=="GO"],breaks="fd",main="SRRT vs. Go RT",freq=F,ylim=c(0,8))
lines(density(data$RT[data$SS=="SS"],na.rm=T),col="red",lwd=2)

###----------------------------------------- Let's start fitting

# Uniform (scaled beta) priors:
p1 <- p.vector; p1[1:length(p1)] <- 1; p1
p.prior <- prior.p.dmc(
  dists = rep("beta",length(p1)),p1=p1,p2=rep(1,length(p1)), # Uniform(0,1)
  lower=rep(0,length(p1)),upper=c(rep(2,2),rep(.5,4),rep(1,2)) # Scale to Uniform(lower,upper)
)
par(mfcol=c(2,4)); for (i in names(p.prior)) plot.prior(i,p.prior,ylim = c(0,4))

# Start sampling
samples <- samples.dmc(nmc=100,p.prior,data)
samples <- run.unstuck.dmc(samples,report=10,cores=4,p.migrate=.05,verbose=TRUE)
layout(1)
plot.dmc(samples,pll.chain=TRUE)

samples2 <- run.converge.dmc(samples.dmc(samples=samples,nmc=50,thin=5),
                             report=10,cores=4,cut=1.1,verbose=TRUE,nmc=50,minN=500,max.try=20)

# Posterior loglikelihood looks nice:
layout(1)
plot.dmc(samples2,pll.chain=TRUE)
# Parameter chains look like fat hairy caterpillars:
plot.dmc(samples2,layout=c(2,4))

# Rhat looks good:
gelman.diag.dmc(samples2)
#Potential scale reduction factors:
#  
#  Point est. Upper C.I.
#mu.true          1.03       1.04
#sigma.true       1.03       1.05
#tau.true         1.02       1.03
#tf               1.01       1.02
#muS              1.01       1.02
#sigmaS           1.03       1.04
#tauS             1.01       1.02
#gf               1.01       1.02
#
#Multivariate psrf
#
#1.05

# Good parameter recovery:
summary.dmc(samples2)
p.vector
#1. Empirical mean and standard deviation for each variable,
#plus standard error of the mean:
#  
#  Mean       SD  Naive SE Time-series SE
#mu.true    0.49891 0.004133 6.888e-05      0.0001674
#sigma.true 0.04897 0.002820 4.700e-05      0.0001125
#tau.true   0.08324 0.004853 8.088e-05      0.0002165
#tf         0.10586 0.055423 9.237e-04      0.0023759
#muS        0.20162 0.025448 4.241e-04      0.0011498
#sigmaS     0.03379 0.022154 3.692e-04      0.0012352
#tauS       0.04818 0.032775 5.463e-04      0.0015743
#gf         0.08682 0.009922 1.654e-04      0.0004765
#
#2. Quantiles for each variable:
#  
#  2.5%     25%     50%     75%   97.5%
#mu.true    0.4908928 0.49609 0.49880 0.50175 0.50693
#sigma.true 0.0436874 0.04701 0.04888 0.05086 0.05461
#tau.true   0.0739193 0.07980 0.08322 0.08655 0.09260
#tf         0.0095323 0.06593 0.10293 0.14164 0.22274
#muS        0.1475448 0.18492 0.20296 0.22085 0.24441
#sigmaS     0.0008918 0.01762 0.03191 0.04636 0.08317
#tauS       0.0026869 0.02271 0.04343 0.06915 0.12211
#gf         0.0676455 0.08001 0.08651 0.09335 0.10714

# Parameter recovery can be summarized as follows:
check.recovery.dmc(samples2,p.vector)
#mu.true sigma.true tau.true   tf  muS sigmaS  tauS    gf
#True              0.50       0.05     0.08 0.10 0.20   0.03  0.05  0.10
#2.5% Estimate     0.49       0.04     0.07 0.01 0.15   0.00  0.00  0.07
#50% Estimate      0.50       0.05     0.08 0.10 0.20   0.03  0.04  0.09
#97.5% Estimate    0.51       0.05     0.09 0.22 0.24   0.08  0.12  0.11
#Median-True       0.00       0.00     0.00 0.00 0.00   0.00 -0.01 -0.01
#attr(,"ci50")
#
#25%        75%
#  mu.true    0.49609155 0.50175421
#sigma.true 0.04701213 0.05086053
#tau.true   0.07980446 0.08654633
#tf         0.06593021 0.14164063
#muS        0.18491995 0.22084618
#sigmaS     0.01761506 0.04635668
#tauS       0.02270701 0.06914597
#gf         0.08001261 0.09335299

# Priors are nicely updated
# Go parameters are better constrained than stop parameters (more data = less uncertainty)
plot.dmc(samples2,layout=c(2,4),p.prior=p.prior)

# Good fits, as expected for simulated data.
pp <- post.predict.dmc(samples2,n.post=1000)
# NB: For Stop-Signal fits need to save off the raw simulations!
pp2 <- post.predict.dmc(samples2,n.post=1000,save.simulation=TRUE)

# Compare cdf of observed and predicted go RTs:
layout(1)
plot.pp.dmc(pp,style="cdf",layout=c(2,2),ylim=c(0,1),fits.pcol="grey",
            model.legend=FALSE,dname="",aname="")

# Then compare observed (black dots) and predicted (grey violins) response rates on the 
# different SSDs using posterior-predictive simulations.
# Also print the number of trials/SSD (n) and corresponding posterior-predictive 
# p value (p = mean(predicted response rate >= observed response rate))

layout(1)
plot_SS_if.dmc(data=samples2$data,sim=pp2)
#Inhibition function 
#SSD  n     p
#1 0.10  2 1.000
#2 0.15  4 0.079
#3 0.20 11 0.606
#4 0.25 37 0.693
#5 0.30 70 0.803
#6 0.35 77 0.700
#7 0.40 39 0.053
#8 0.45  8 0.796
#9 0.50  2 0.719

# Note that this function also provides plots of inhibition functions averaged
# over subjects, in which case the data argument should be a samples object for
# multiple subjects. By default the averaging is of percentiles (i.e., for each 
# subject get the percentile cut points of their SSD distribution, then average 
# particiapnts' probabilty of inhibition within each percentile range). This 
# default choice was guided by the observation that individual differences cause 
# averages on absolute SSD values to be flat, unlike the inhibition functions 
# for any individual.

# Finally compare observed (black dots) and predicted (grey violins) median signal-respond 
# RTs using posterior-predictive simulations. 

# For each SSD for which the median signal-respond RT can be computed,
# also print the number of trials/SSD (n), the number of observed signal-respond RTs/SSD (nrt),
# the corresponding posterior-predictive p value (p), and the number of simulated 
# datasets for which median signal-respond RT can be computed:
layout(1)
plot_SS_srrt.dmc(data=samples2$data,sim=pp2)
#Median signal-respond RTs 
#SSD nrt     p n.sim
#1 0.15   2 0.436   436
#2 0.20   2 0.581   890
#3 0.25   9 0.801  1000
#4 0.30  26 0.555  1000
#5 0.35  43 0.181  1000
#6 0.40  33 0.789  1000
#7 0.45   6 0.514  1000
#8 0.50   2 0.085   976

# In this case the default for averaging is to pool all of the SSDs across 
# subjects then get the cut points between which to calculate each participant’s
# SRRTs from the percentiles of that pooled SSD distribution. The percentile.av
# method used for the inhibition function is also available here but it is not the default
# as it tends to produce SRRT functions that tend to be flatter than observed for
# any indvidual.

###----------------------------------------------------------------------In-class exercise 1:

# Compare the observed vs. predicted omission rates on Go trials.

pp.go <- subset(pp2,SS=="GO")
data.go <- subset(samples2$data,SS=="GO")
# Predicted:
pp.omission <- with(pp.go,tapply(RT,reps,function(x) mean(is.na(x))))
# Observed
omission <- mean(is.na(data.go$RT))
# Plot
layout(1)
hist(pp.omission)
points(omission,1,col="red",pch=16,cex=2)
# Posterior predictive p-value:
mean(pp.omission>=omission)

###----------------------------------------------------------------------In class excercise 2:

# Fit the data from sample2 with a misspecified model that does not feature trigger failures.
# How does this effect the parameter estimates?

# Specify the model without trigger failures:
modelMis <- model.dmc(
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),
  responses=c("NR","r1","r2"),
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")),
  p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1"),
  constants=c(mu.false=1e6,sigma.false=.001,tau.false=.001,tf=0), #No trigger failures
  type="exgss")

p.vectorMis <- c(mu.true=NA,muS=NA,sigma.true=NA,sigmaS=NA,tau.true=NA,tauS=NA,gf=NA) 

# Specify the data:
dataMis <- data.model.dmc(samples2$data,modelMis)

# As before, uniform (scaled beta) priors, but now we don't need prior for tf:
p1Mis <- p.vectorMis; p1Mis[1:length(p1Mis)] <- 1; p1Mis
p.priorMis <- prior.p.dmc(
  dists = rep("beta",length(p1Mis)),p1=p1Mis,p2=rep(1,length(p1Mis)), # Uniform(0,1)
  lower=rep(0,length(p1Mis)),upper=c(rep(2,2),rep(.5,4),1) # Scale to Uniform(lower,upper)
)
par(mfcol=c(2,4)); for (i in names(p.priorMis)) plot.prior(i,p.priorMis,ylim = c(0,4))

# Start sampling:
samples3 <- samples.dmc(nmc=100,p.priorMis,dataMis)
samples3 <- run.unstuck.dmc(samples3,report=10,cores=4,p.migrate=.05,verbose=TRUE)
layout(1)
plot.dmc(samples3,pll.chain=TRUE,density=FALSE,smooth=FALSE)

samples4 <- run.converge.dmc(samples.dmc(samples=samples3,nmc=50,thin=5),report=10,cores=4,cut=1.1,verbose=TRUE,nmc=50,minN=500,max.try=20)

# Posterior loglikelihood looks nice:
layout(1)
plot.dmc(samples4,pll.chain=TRUE,density=FALSE,smooth=FALSE)
# Parameter chains look like fat hairy catapillars:
plot.dmc(samples4,layout=c(2,4),smooth=FALSE,density=FALSE)

# Rhat looks good:
gelman.diag.dmc(samples4)

# Biased stop estimates:
summary.dmc(samples4)
p.vector

# Priors are updated but posteriors don't recover true values very well:
plot.dmc(samples4,layout=c(2,4),p.prior=p.priorMis)
p.vector

# Compare to true model:
par(mfcol=c(2,4))
for(i in names(p.vectorMis)){
  hist(samples2$theta[,i,],main=i,freq=F,breaks="fd")
  lines(density(samples4$theta[,i,]),col="red",lwd=2)
}

ppMis <- post.predict.dmc(samples4,n.post=1000)
ppMis2 <- post.predict.dmc(samples4,n.post=1000,save.simulation=TRUE)

# Similar predictions for Go cdfs (no effect on go parameters):
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

# save_data(data,samples,samples2,pp,pp2,dataMis,samples3,samples4,ppMis,ppMis2,file="dmc_6_4b.RData")


