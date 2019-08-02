##################  DMC Lesson 6: More Models

### Lesson 6.2  Sampling and assessing a single Wald subject  (with go failure)

# The shifted Wald model is discussed in Heathcote, A. (2004). Fitting the Wald 
# and Ex-Wald Distributions to Response Time Data. Behaviour Research Methods, 
# Instruments & Computers, 36, 678-694. 
# http://www.weebly.com/uploads/4/9/3/3/49339445/26_.pdf
# We also added "go failure" to this model (see dmc_6_1_LBA_GF.R).

rm(list=ls()) 

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("Wald","wald.R")

# load_data ("dmc_6_2.RData")


# Start with no go failure model 
model <- model.dmc(p.map=list(A="1",B="1",v="M",t0="1",gf="1"),
                   constants=c(gf=-7),match.map=list(M=list(s1=1,s2=2)),
                   factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="wald")
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"       "B"       "v.true"  "v.false" "t0"     
# 
# Constants are (see attr(,"constants") ):
# gf 
# -7 
# 
# Model type = wald 


# Simulate some data, with around 75% accuracy
p.vector  <- c(A=.5,B=1,v.true=2,v.false=.75,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)

# Look reasonable
par(mfrow=c(1,2))
plot.cell.density(data.cell=data.model[data.model$S=="s1",],C="r1",xlim=c(0,2))
plot.cell.density(data.cell=data.model[data.model$S=="s2",],C="r2",xlim=c(0,2))

crct <- (data.model$S=="s1" & data.model$R=="r1") |
  (data.model$S=="s2" & data.model$R=="r2")
round(tapply(crct,list(data.model$S),mean),2)
#   s1   s2 
# 0.76 0.77 
# Fast errors I CANT FIND ANY SETTING WHERE THEY AR SLOW EVEN A=0, bigger A
# makes them faster and faster.
round(tapply(data.model$RT,list(data.model$S,C=crct),mean),2)
#     C
#      FALSE TRUE
#   s1  0.72 0.73
#   s2  0.73 0.73

head(likelihood.dmc(p.vector,data.model))
# [1] 0.4683722 1.3918407 1.3745176 0.2904727 0.2519452 0.7266601

# Check profiles, they all look OK
par(mfrow=c(2,3))
profile.dmc("A",      .1,  2,p.vector,data.model,ylim=NA)
profile.dmc("B",      .1,  2,p.vector,data.model,ylim=NA)
profile.dmc("v.true",  1,  3,p.vector,data.model,ylim=NA)
profile.dmc("v.false",.5, 1.5,p.vector,data.model,ylim=NA)
profile.dmc("t0",     .01,.4,p.vector,data.model,ylim=NA)


# Give t0 a uniform prior from 0.1-1s, other priors normal, truncated
# below for A and B as they must be positive, unbounded for the v
p.prior <- prior.p.dmc(
  dists = rep("tnorm",5),
  p1=c(A=1,B=1,v.true=1,v.false=.5,t0=.5),                           
  p2=c(1,1,3,3,1),lower=c(0,0,0,0,.1),upper=c(NA,NA,NA,NA,1)
)
# par(mfcol=c(2,3)); for (i in names(p.prior)) plot.prior(i,p.prior)

# Try estimation
samples <- samples.dmc(nmc=100,p.prior,data.model)
samples1 <- run.unstuck.dmc(samples,p.migrate=.05,verbose=TRUE,cores=15)

samples2 <- run.converge.dmc(samples.dmc(samples=samples1,nmc=100),
                             nmc=50,verbose=TRUE,max.try=10,cores=15)
# [1] "Final multivariate psrf = 1.07547366981695"
# Effective sample size
#       A       B  v.true v.false      t0 
#     235     241     196     236     220 
samples3 <- run.dmc(samples.dmc(samples=samples2,nmc=100,thin=5),report=1,cores=15)    

plot.dmc(samples3,pll.chain=TRUE)
plot.dmc(samples3,layout=c(2,3))
plot.dmc(samples3,layout=c(2,3),p.prior=p.prior)

# Poor recovery of A
check.recovery.dmc(samples3,p.vector)
#                    A    B v.true v.false    t0
# True            0.50 1.00   2.00    0.75  0.20
# 2.5% Estimate   0.05 0.92   1.95    0.70  0.18
# 50% Estimate    0.38 1.06   1.99    0.75  0.19
# 97.5% Estimate  0.62 1.25   2.03    0.79  0.21
# Median-True    -0.12 0.06  -0.01    0.00 -0.01\


# Very high negative correlation between A and B and positive between A and t0
pairs.dmc(samples3)

#######  GO FAILURE

modelGF <- model.dmc(p.map=list(A="1",B="1",v="M",t0="1",gf="1"),
                     match.map=list(M=list(s1=1,s2=2)),
                     factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="wald")
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"       "B"       "v.true"  "v.false" "t0"      "gf"     
# 
# Constants are (see attr(,"constants") ):
# numeric(0)
# 
# Model type = wald 


# Simulate some data, with around 75% accuracy, 15.866% go failure
p.vector  <- c(A=.5,B=1,v.true=2,v.false=.75,t0=.2,gf=-1)
data.model <- data.model.dmc(simulate.dmc(p.vector,modelGF,n=1e4),modelGF)

# Look reasonable
par(mfrow=c(1,2))
plot.cell.density(data.cell=data.model[data.model$S=="s1",],C="r1",xlim=c(0,2))
plot.cell.density(data.cell=data.model[data.model$S=="s2",],C="r2",xlim=c(0,2))


# Check profiles, they all look OK
par(mfrow=c(2,3))
profile.dmc("A",      .1,  2,p.vector,data.model,ylim=NA)
profile.dmc("B",      .1,  2,p.vector,data.model,ylim=NA)
profile.dmc("v.true",  1,  3,p.vector,data.model,ylim=NA)
profile.dmc("v.false",.5, 1.5,p.vector,data.model,ylim=NA)
profile.dmc("t0",     .01,.4,p.vector,data.model,ylim=NA)
profile.dmc("gf",     -1.5,-.5,p.vector,data.model,ylim=NA)


# Give t0 a uniform prior from 0.1-1s, other priors normal, truncated
# below for A and B as they must be positive, unbounded for the v
p.prior <- prior.p.dmc(
  dists = rep("tnorm",6),
  p1=c(A=1,B=1,v.true=1,v.false=.5,t0=.5,gf=-1.5),                           
  p2=c(1,1,3,3,1,2),lower=c(0,0,0,0,.1,NA),upper=c(NA,NA,NA,NA,1,NA)
)
# par(mfcol=c(2,3)); for (i in names(p.prior)) plot.prior(i,p.prior)

# Try estimation
samplesGF <- samples.dmc(nmc=100,p.prior,data.model)
samplesGF1 <- run.unstuck.dmc(samplesGF,p.migrate=.05,verbose=TRUE,cores=18)

samplesGF2 <- run.converge.dmc(samples.dmc(samples=samplesGF1,nmc=100),
                               nmc=50,verbose=TRUE,max.try=10,cores=18)
# [1] "Final multivariate psrf = 1.09879647913292"
# Effective sample size
#       A       B  v.true v.false      t0      gf 
#     465     480     505     525     497     527

samplesGF3 <- run.dmc(samples.dmc(samples=samplesGF2,nmc=100,thin=5),report=1,cores=18)    

plot.dmc(samplesGF3,pll.chain=TRUE)
plot.dmc(samplesGF3,layout=c(2,3))
plot.dmc(samplesGF3,layout=c(2,3),p.prior=p.prior)

# Better A recovery on average but CI still wide. 
check.recovery.dmc(samplesGF3,p.vector)
#                   A     B v.true v.false   t0    gf
# True           0.50  1.00   2.00    0.75 0.20 -1.00
# 2.5% Estimate  0.13  0.86   1.93    0.69 0.19 -1.02
# 50% Estimate   0.54  0.96   1.98    0.74 0.21 -0.99
# 97.5% Estimate 0.73  1.18   2.03    0.79 0.22 -0.98
# Median-True    0.04 -0.04  -0.02   -0.01 0.01  0.01

# Very high negative correlation between A and B and positive between A and t0
pairs.dmc(samplesGF3)


######### SIMPLE RT = 100% accuracy case (with two stimuli)

modelS <- model.dmc(p.map=list(A="1",B="1",v="M",t0="1",gf="1"),
                    match.map=list(M=list(s1=1,s2=2)),
                    constants=c(v.false=-1,gf=-7), # ANY negative value => Inf RT/never wins
                    factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="wald")
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"      "B"      "v.true" "t0"    
# 
# Constants are (see attr(,"constants") ):
# v.false      gf 
#      -1      -7 
# 
# Model type = wald 


# Simulate some data
p.vector  <- c(A=.5,B=1,v.true=2,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),modelS)

# Look reasonable
par(mfrow=c(1,2))
plot.cell.density(data.cell=data.model[data.model$S=="s1",],C="r1",xlim=c(0,2))
plot.cell.density(data.cell=data.model[data.model$S=="s2",],C="r2",xlim=c(0,2))


head(likelihood.dmc(p.vector,data.model))
# [1] 0.4683722 1.3918407 1.3745176 0.2904727 0.2519452 0.7266601

# Check profiles, they all look OK
par(mfrow=c(2,2))
profile.dmc("A",      .1,  2,p.vector,data.model,ylim=NA)
profile.dmc("B",      .1,  2,p.vector,data.model,ylim=NA)
profile.dmc("v.true",  1,  3,p.vector,data.model,ylim=NA)
profile.dmc("t0",     .01,.4,p.vector,data.model,ylim=NA)


# Give t0 a uniform prior from 0.1-1s, other priors normal, truncated
# below for A and B as they must be positive, unbounded for the v
p.prior <- prior.p.dmc(
  dists = rep("tnorm",4),
  p1=c(A=1,B=1,v.true=1,t0=.5),                           
  p2=c(1,1,3,1),lower=c(0,0,0,.1),upper=c(NA,NA,NA,1)
)
# par(mfcol=c(2,2)); for (i in names(p.prior)) plot.prior(i,p.prior)

# Try estimation
samplesS <- samples.dmc(nmc=100,p.prior,data.model)
samplesS1 <- run.unstuck.dmc(samplesS,p.migrate=.05,verbose=TRUE,cores=12)

samplesS2 <- run.converge.dmc(samples.dmc(samples=samplesS1,nmc=100),
                              nmc=50,verbose=TRUE,max.try=10,cores=12)
# [1] "Final multivariate psrf = 1.09758215308722"
# Effective sample size
#      A      B v.true     t0 
#    185    188    221    190

samplesS3 <- run.dmc(samples.dmc(samples=samplesS2,nmc=100,thin=5),report=1,cores=12)    

plot.dmc(samplesS3,pll.chain=TRUE)
plot.dmc(samplesS3,layout=c(2,2))
plot.dmc(samplesS3,layout=c(2,2),p.prior=p.prior)

# Poor recovery of A
check.recovery.dmc(samplesS3,p.vector)
#                    A    B v.true    t0
# True            0.50 1.00   2.00  0.20
# 2.5% Estimate   0.03 0.88   1.97  0.18
# 50% Estimate    0.40 1.06   2.00  0.19
# 97.5% Estimate  0.68 1.27   2.04  0.21
# Median-True    -0.10 0.06   0.00 -0.01

# Very high negative correlation between A and B and positive between A and t0
pairs.dmc(samplesS3)


# save_data (samples3,samplesGF3,samplesGF3,file="dmc_6_2.RData")
