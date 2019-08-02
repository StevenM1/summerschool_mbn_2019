##################  DMC Lesson 6: More Models

### Lesson 6.5: ExGaussian stop-signal model with 1 stop & 2 go accumulators 


### An EX-Gaussian stop-signal example with TRIGGER and GO FAILURE and 
# CONTEXT INDEPENDENT paramterization (i.e., seperate stop accumulator parameters).
# This example uses a 3 accumualtor race suitable for high error rates on the go task.
# It also uses a probit scale for the tf and gf parametrs, and truncated normal
# priors (this could also be done with the 2 accumualtor modle), which work best 
# with heirarchical models (although these are not fit here).
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
# probit scale we set it to -6 (pnorm(-6)=9.865876e-10)
if (!is.tf & !is.gf) {
  
  model <- model.dmc(
    # SS stands for trial type (GO or Stop-signal [SS]):
    factors=list(S=c("s1","s2"),SS=c("GO","SS")),
    # NR stands for "No response", i.e., go omission & succesfull inhibitions:
    responses=c("NR","r1","r2"),
    # Match scores correct responses for each GO stimulus as usual, and also scores 
    # the "correct" stimulus corresponding to an NR response, but the latter has
    # no effect (it just avoids a standard check making sure that each response is scored): 
    match.map=list(M=list(s1="r1",s2="r2",s1="NR")),
    p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1"),
    # No trigger failures and no go failures:            
    constants=c(tf=-6,gf=-6),
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
    # NR stands for "No response", i.e., go omission & succesfull inhibitions:
    responses=c("NR","r1","r2"),
    # Match scores correct responses for each GO stimulus as usual, and also scores 
    # the "correct" stimulus corresponding to an NR response, but the latter has
    # no effect (it just avoids a standard check making sure that each response is scored): 
    match.map=list(M=list(s1="r1",s2="r2",s1="NR")),
    p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1"),
    # No go failures: 
    constants=c(gf=-6),
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
    # NR stands for "No response", i.e., go omission & succesfull inhibitions:
    responses=c("NR","r1","r2"),
    # Match scores correct responses for each GO stimulus as usual, and also scores 
    # the "correct" stimulus corresponding to an NR response, but the latter has
    # no effect (it just avoids a standard check making sure that each response is scored): 
    match.map=list(M=list(s1="r1",s2="r2",s1="NR")),
    p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1"),
    # No trigger failures: 
    constants=c(tf=-6),
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
    # NR stands for "No response", i.e., go omission & succesfull inhibitions:
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
# 0.15  0.2      0.25       0.3      0.35       0.4      0.45        Inf 
# NA  0.4169316 0.5367189 0.5176544 0.5497371 0.5468397 0.5735215 0.5588533

# Plot median signal-respond RT per SSD:
plot_SS_srrt.dmc(data) # Median SRRT increases as a function of SSD, as it should

# Show number of signal-respond RT per SSD:
Nr = tapply(!is.na(data$RT),data[,c("SS","SSD")],sum)[2,]
Nr
# 0.15  0.2 0.25  0.3 0.35  0.4 0.45  Inf 
#   0    1   10   32   47   30    5   NA 

# Signal-respond RTs should be faster than go RTs:
hist(data$RT[data$SS=="GO"],breaks="fd",main="SRRT vs. Go RT",freq=F,ylim=c(0,8))
lines(density(data$RT[data$SS=="SS"],na.rm=T),col="red",lwd=2)

###----------------------------------------- Let's start fitting

# Truncated normal priors
p1 <- p.vector; p1[1:length(p1)] <- 1; p1[10:11] <- 3
p.prior <- prior.p.dmc(
  dists = rep("tnorm",length(p1)),p1=p.vector,p2=p1, 
  lower=c(rep(0,length(p1)-2),-6,-6),upper=c(rep(2,3),rep(.5,6),rep(6,2)) 
)
par(mfcol=c(2,6)); for (i in names(p.prior)) plot.prior(i,p.prior)

# Start sampling
samples <- samples.dmc(nmc=100,p.prior,data)
samples <- run.unstuck.dmc(samples,report = 10,cores=33,p.migrate=.05,verbose=TRUE)
layout(1)
plot.dmc(samples,pll.chain=TRUE)

samples2 <- run.converge.dmc(samples.dmc(samples=samples,nmc=50,thin=5),
  report=10, cores=33,cut=1.1,verbose=TRUE,nmc=50,minN=500,max.try=20)

# Posterior loglikelihood looks nice:
layout(1)
plot.dmc(samples2,pll.chain=TRUE,density=FALSE,smooth=FALSE)
# Parameter chains look like fat hairy catapillars:
plot.dmc(samples2,layout=c(2,6))

# Rhat looks good:
gelman.diag.dmc(samples2)
#             Point est. Upper C.I.
# mu.true           1.05       1.07
# mu.false          1.03       1.04
# sigma.true        1.06       1.08
# sigma.false       1.04       1.05
# tau.true          1.05       1.06
# tau.false         1.04       1.05
# tf                1.05       1.07
# muS               1.04       1.06
# sigmaS            1.02       1.03
# tauS              1.03       1.05
# gf                1.04       1.05
# 
# Multivariate psrf
# 
# 1.08

# Good parameter recovery except tf
check.recovery.dmc(samples2,p.vector)
#                mu.true mu.false sigma.true sigma.false tau.true
# True              0.50     0.60       0.05        0.03     0.08
# 2.5% Estimate     0.50     0.59       0.05        0.02     0.05
# 50% Estimate      0.51     0.60       0.05        0.03     0.07
# 97.5% Estimate    0.52     0.61       0.06        0.03     0.09
# Median-True       0.01     0.00       0.00        0.00    -0.01
#                tau.false    tf  muS sigmaS tauS    gf
# True                0.04 -1.28 0.20   0.03 0.05 -1.28
# 2.5% Estimate       0.02 -5.72 0.17   0.00 0.01 -1.47
# 50% Estimate        0.04 -2.89 0.21   0.02 0.05 -1.35
# 97.5% Estimate      0.05 -1.04 0.24   0.07 0.10 -1.23
# Median-True         0.00 -1.61 0.01  -0.01 0.00 -0.07

# Note as tf is on the probit scale the miss is 
round(pnorm(c(-5.72,-2.89,-1.04)),4)
# [1] 0.0000 0.0019 0.1492


# Priors are nicely updated, except tf
plot.dmc(samples2,layout=c(2,6), p.prior=p.prior, show.obs=FALSE)

# Fits 
pp <- post.predict.dmc(samples2,n.post=200,save.simulation=TRUE)

layout(1)
plot_SS_if.dmc(data=samples2$data,sim=pp)
# Inhibition function 
#    SSD  n     p
# 1 0.20 14 1.000
# 2 0.25 43 0.140
# 3 0.30 78 0.940
# 4 0.35 78 0.530
# 5 0.40 32 0.465
# 6 0.45  5 0.550

layout(1)
plot_SS_srrt.dmc(data=samples2$data,sim=pp)
# Median signal-respond RTs 
#    SSD nrt     p n.sim
# 1 0.25  15 0.535   200
# 2 0.30  27 0.480   200
# 3 0.35  51 0.705   200
# 4 0.40  27 0.000   200
# 5 0.45   5 0.250   200


# save_data(data,samples,samples2,pp,file="dmc_6_5_EXG3probit.RData")
