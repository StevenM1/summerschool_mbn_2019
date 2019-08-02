##################  DMC Lesson 5: Miscellaneous models


### Lesson 4.5 constant: Hierarchical model where some hyper-parameters 
### are fixed: LNR example.  

rm(list=ls())

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model("LNR","lnr.R")

# load_data("dmc_4_5constant.RData")

# Setup a simple 4 parameter model
model <- model.dmc(p.map=list(meanlog="M",t0="1",sdlog="1",st0="1"),
                   match.map=list(M=list(s1=1,s2=2)),constants=c(st0=0),
                   factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="lnr")

# Setup a population distribution
p.mu  <- c(meanlog.true=-1,meanlog.false=0,t0=.3,sdlog=1) 
p.sigma <- c(.2,.2,.05,.3); names(p.sigma) <- names(p.mu)
p.prior <- prior.p.dmc(p1=p.mu,p2=p.sigma,lower=c(NA,NA,0.1,0))
# par(mfcol=c(2,2)); for (i in names(p.prior)) plot.prior(i,p.prior)

# Simulate a sample of 40 subjects with 500 data points each.
raw.data <- h.simulate.dmc(model,p.prior=p.prior,n=250,ns=40)

# parameters
ps <- attr(raw.data,"parameters")
round(ps,2)
round(apply(ps,2,mean),3)
#  meanlog.true meanlog.false            t0         sdlog 
#        -1.038        -0.005         0.288         0.962 
round(apply(ps,2,sd),3)
#  meanlog.true meanlog.false            t0         sdlog 
#         0.211         0.194         0.055         0.285 

# Combine with true model.
data.model <- data.model.dmc(raw.data,model)


# When a hyper-parameter is designated as a constant its value is taken from the
# hyper prior. In this example we set the usual broad hyper-priors for all 
# except the scale of sdlog, which we set to exactly match the simulated 
# population value, as we will fix this value then examine parameter recovery. 

# Location prior, somewhat different from true parameters
p.mu.mu <- c(meanlog.true=0,meanlog.false=0,t0=0.2,sdlog=0.5)
p.mu.sigma <- p.mu.mu; p.mu.sigma[1:4] <- c(3,3,1,3)
mu.prior <- prior.p.dmc(p1=p.mu.mu,p2=p.mu.sigma,lower=c(NA,NA,0,0))
par(mfcol=c(2,2)); for (i in names(mu.prior)) plot.prior(i,mu.prior)

# We indicate the sdlog has prior distribution of type "constant". For this 
# distribution only the first (p1) parameter has an effect, determining the 
# fixed value (a value must be specified for p2 but it is ignored, here we use 
# NA). We use exponential priors for remaining parameters.
p.sigma.shape <- p.sigma.scale <- p.mu.mu; 
p.sigma.shape[1:4] <- c(rep(1,3),.3)
p.sigma.scale[1:4] <- c(1,1,.3,NA)
sigma.prior <- prior.p.dmc(p1=p.sigma.shape,p2=p.sigma.scale,
                           dists=c(rep("gamma",3),"constant"))
# Note that plot.prior will ignore a constant prior
par(mfcol=c(2,2)); for (i in names(sigma.prior)) plot.prior(i,sigma.prior)

# Make a hyper-prior list
pp.prior=list(mu.prior,sigma.prior)


hsamples <- h.samples.dmc(nmc=100,data=data.model,
                          thin=5,p.prior=p.prior,pp.prior=pp.prior)

hsamples <- h.run.unstuck.dmc(hsamples,cores=4,report=1,
                              p.migrate=0.05,h.p.migrate=0.05)
# Run without migration
hsamples1 <- h.run.converge.dmc(h.samples.dmc(nmc=100,samples=hsamples),
                                nmc=50,cores=4,verbose=TRUE)

h.gelman.diag.dmc(hsamples1)
# hyper     5    16    18    28    25    32    26    19     4    36    17 
#  1.03  1.02  1.02  1.02  1.02  1.03  1.03  1.03  1.03  1.03  1.03  1.03 
#    37    13    35    24    14    20    11     8    12     6     3    34 
#  1.03  1.03  1.04  1.04  1.04  1.04  1.04  1.04  1.04  1.04  1.04  1.04 
#    40    10    23    39    33     9    38    21    22     7    30     2 
#  1.04  1.04  1.04  1.04  1.04  1.04  1.04  1.04  1.04  1.04  1.05  1.05 
#    27    29    31    15     1 
#  1.05  1.05  1.06  1.06  1.06

# Stable chains.
plot.dmc(hsamples1,hyper=TRUE)

# Good number of samples
effectiveSize.dmc(hsamples1,hyper=TRUE)
#  meanlog.true.h1 meanlog.false.h1            t0.h1         sdlog.h1 
#              777              839              873             1122 
#  meanlog.true.h2 meanlog.false.h2            t0.h2 
#              860              912              861 

# Good fit
pp <- h.post.predict.dmc(hsamples1)
plot.pp.dmc(pp)

# Good parameter recovery
hest <- summary.dmc(hsamples1,hyper=TRUE)
# Mean estimates 
round(hest$statistics[,"Mean"][1:4],3)
#  meanlog.true.h1 meanlog.false.h1            t0.h1         sdlog.h1 
#           -1.042            0.002            0.288            0.957 
# accurate
round(hest$statistics[,"Mean"][1:4]-apply(ps,2,mean),3)
#  meanlog.true.h1 meanlog.false.h1            t0.h1         sdlog.h1 
#           -0.004            0.007            0.000           -0.004 

# And scale 
round(hest$statistics[,"Mean"][5:7],3)
#  meanlog.true.h2 meanlog.false.h2            t0.h2 
#            0.230            0.213            0.057 
# quite accurate
round(hest$statistics[,"Mean"][5:7]-apply(ps,2,sd)[-4],3)
#  meanlog.true.h2 meanlog.false.h2            t0.h2 
#            0.019            0.019            0.002 

### Can also make a location parameter a constant, here we do it with t0, setting
# it at the true value in the prior
p.mu.mu <- c(meanlog.true=0,meanlog.false=0,t0=0.3,sdlog=0.5)
# Broad scales
p.mu.sigma <- p.mu.mu; p.mu.sigma[1:4] <- c(3,3,NA,3)
mu.prior <- prior.p.dmc(p1=p.mu.mu,p2=p.mu.sigma,lower=c(NA,NA,0,0),
                        dists=c("tnorm","tnorm","constant","tnorm"))
# par(mfcol=c(2,2)); for (i in names(mu.prior)) plot.prior(i,mu.prior)

# We keep the fixed sdlog scale as before, showing that arbitary mixtures of
# constants work.
p.sigma.shape <- p.sigma.scale <- p.mu.mu; 
p.sigma.shape[1:4] <- c(rep(1,3),.3)
p.sigma.scale[1:4] <- c(1,1,.3,NA)
sigma.prior <- prior.p.dmc(p1=p.sigma.shape,p2=p.sigma.scale,
                           dists=c(rep("gamma",3),"constant"))
# par(mfcol=c(2,2)); for (i in names(sigma.prior)) plot.prior(i,sigma.prior)

# Make a hyper-prior list
pp.priorA=list(mu.prior,sigma.prior)


hsamplesA <- h.samples.dmc(nmc=50,data=data.model,
                           thin=5,p.prior=p.prior,pp.prior=pp.priorA)
hsamplesA <- h.run.unstuck.dmc(hsamplesA,cores=4,report=1,
                               p.migrate=0.05,h.p.migrate=0.05)
# Run without migration, verobse=TRUE provides the summary of convergence 
# copied below.
hsamplesA1 <- h.run.converge.dmc(h.samples.dmc(nmc=100,samples=hsamplesA),
                                 nmc=50,cores=4,verbose=TRUE)
# Iterations = 100, Effective N = NA, Hyper mpsrf = 1.04
# Subject mpsrf achieved (sorted):
#   15   32   29    7   33   13    9   34   10   39   19   21   36   38 
# 1.08 1.08 1.07 1.07 1.07 1.07 1.07 1.07 1.07 1.07 1.06 1.06 1.06 1.06 
#    1   23    4   28   40   22    2    8   18   24   11    3   12    5 
# 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.05 1.05 1.05 1.05 1.05 
#   25    6   26   20   35   37   14   30   17   16   31   27 
# 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.04 1.04 1.04 1.04 

# Chains are nice
plot.dmc(hsamplesA1,hyper=TRUE,layout=c(2,3))

# Fits are good
ppA <- h.post.predict.dmc(hsamplesA1)
plot.pp.dmc(ppA)

hest <- summary.dmc(hsamplesA1,hyper=TRUE)
# Mean estimates accurate
round(hest$statistics[,"Mean"][1:3],3)
#  meanlog.true.h1 meanlog.false.h1         sdlog.h1 
#           -1.040            0.004            0.958 
round(hest$statistics[,"Mean"][1:3]-apply(ps,2,mean)[-3],3)
#  meanlog.true.h1 meanlog.false.h1         sdlog.h1 
#           -0.002            0.009           -0.004 

# And scale 
round(hest$statistics[,"Mean"][5:7],3)
#  meanlog.true.h2 meanlog.false.h2            t0.h2 
#            0.230            0.213            0.058 
round(hest$statistics[,"Mean"][5:7]-apply(ps,2,sd)[-4],3)
#  meanlog.true.h2 meanlog.false.h2            t0.h2 
#            0.019            0.019            0.004 


### Can even have both a scale and location parameter a constant, here t0, again
# setting the prior to the true value to check recovery.
# RUNNING THIS EXAMPLE IS LEFT AS AN EXERCISE
p.mu.mu <- c(meanlog.true=0,meanlog.false=0,t0=0.3,sdlog=0.5)
# Broad scales
p.mu.sigma <- p.mu.mu; p.mu.sigma[1:4] <- c(3,3,NA,3)
mu.prior <- prior.p.dmc(p1=p.mu.mu,p2=p.mu.sigma,lower=c(NA,NA,0,0),
                        dists=c("tnorm","tnorm","constant","tnorm"))
# par(mfcol=c(2,2)); for (i in names(mu.prior)) plot.prior(i,mu.prior)

# We keep the fixed sdlog scale as before.
p.sigma.shape <- p.sigma.scale <- p.mu.mu; 
p.sigma.shape[1:4] <- c(rep(1,2),.05,.3)
p.sigma.scale[1:4] <- c(1,1,NA,NA)
sigma.prior <- prior.p.dmc(p1=p.sigma.shape,p2=p.sigma.scale,
                           dists=c(rep("gamma",2),"constant","constant"))
# par(mfcol=c(2,2)); for (i in names(sigma.prior)) plot.prior(i,sigma.prior)

# Make a hyper-prior list
pp.priorB=list(mu.prior,sigma.prior)

hsamplesB <- h.samples.dmc(nmc=50,data=data.model,
                           thin=5,p.prior=p.prior,pp.prior=pp.priorB)
hsamplesB <- h.run.unstuck.dmc(hsamplesA,cores=12,report=1,
                               p.migrate=0.05,h.p.migrate=0.05)
hsamplesB1 <- h.run.converge.dmc(h.samples.dmc(nmc=100,samples=hsamplesB),
                                 nmc=50,cores=12,verbose=TRUE)

plot.dmc(hsamplesB1,hyper=TRUE,layout=c(2,3))

ppB <- h.post.predict.dmc(hsamplesB1)
plot.pp.dmc(ppB)


hest <- summary.dmc(hsamplesB1,hyper=TRUE)
# Mean estimates 
round(hest$statistics[,"Mean"][1:3],3)
round(hest$statistics[,"Mean"][1:3]-apply(ps,2,mean)[-3],3)
# And scale 
round(hest$statistics[,"Mean"][4:6],3)
round(hest$statistics[,"Mean"][4:6]-apply(ps,2,sd)[-c(3:4)],3)


save_data(raw.data,pp.prior,hsamples, hsamples1, 
          pp.priorA, hsamplesA, hsamplesA1, pp, ppA, file="dmc_4_5constant.RData")
