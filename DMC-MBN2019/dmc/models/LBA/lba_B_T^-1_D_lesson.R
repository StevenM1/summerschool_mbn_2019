### Hierarchical LBA Model 6 distance/time cell levels passed as constants.


rm(list=ls())

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("LBA","lba_B_T^-1_D.R")

# load_data("T^-1_D.RData")

# There are 6 cells that will be mapped to mean_v in lba_B_T^-1_D.R
DT <- c("1112","1121","1122","1221","1222","2121")
# Each digit in these levels is used to lookup a numeric value for covariates 
# used to calcuate mean_v for each accumulator in the cell. In this case the
# the first two digits (1000s and 100s) refer to T(ime) and D(istance)
# covariates for accumulator A, and the last two (10s and units) the same for
# accumulator B.
consts <- as.numeric(DT); names(consts) <- paste("mean_v",DT,sep=".")
consts <- c(consts,c(sd_v=1,st0=0))

# Note that no match.map is specified as mapping to the appropriate 
# accumulator (A or B) handled in lba_B_T^-1_D.R. 
model <- model.dmc(p.map = list(A="1",B="R",t0="1",sd_v="1",st0="1",
                                Wt="1",Wd="1",mean_v="DT"), 
          factors=list(DT=DT),
          constants = consts, 
          responses = c("r1","r2"),
          type="norm")
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"    "B.r1" "B.r2" "t0"   "Wt"   "Wd"  
# 
# Constants are (see attr(,"constants") ):
# mean_v.1112 mean_v.1121 mean_v.1122 mean_v.1221 mean_v.1222 mean_v.2121 
#        1112        1121        1122        1221        1222        2121 
#        sd_v         st0 
#           1           0 
# 
# Model type = norm (posdrift= TRUE ) 
# 
p.vector <- c(A=1,B.r1=1,B.r2=2,t0=.2,Wt=1,Wd=.1)

# Simulate large data set
data <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
          
# Likelihoods look sensible
head(likelihood.dmc(p.vector,data))

# Profiles recover parameters
ylim=c(-50000,-25000)
par(mfrow=c(2,3))
profile.dmc("A",.1,2,p.vector,data,ylim=ylim)
profile.dmc("B.A",.1, 2,p.vector,data,ylim=ylim)
profile.dmc("B.B",1,3,p.vector,data,ylim=ylim)
profile.dmc("t0",.1,.3,p.vector,data,ylim=ylim)
profile.dmc("Wt",.1,2,p.vector,data,ylim=ylim)
profile.dmc("Wd",.01,.2,p.vector,data,ylim=ylim)

### Do sampling

# Vague prior
p.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","beta","tnorm","tnorm"),
  p1=c(A=.5,B.r1=1.5,B.r2=1.5,t0=1,Wt=.5,Wd=.5),                           
  p2=c(2,2,2,1,1,1),lower=c(0,0,0,.1,0,0),upper=c(NA,NA,NA,1,NA,NA)
)
par(mfcol=c(2,3)); for (i in names(p.prior)) plot.prior(i,p.prior)

# Setup samples
samples <- samples.dmc(nmc=400,p.prior,data)

# Run with migration
samples <- run.dmc(samples, report = 5, cores=12,p.migrate=.05)
plot.dmc(samples,pll.chain=TRUE)
# Looks OK from 360 or so
plot.dmc(samples,pll.chain=TRUE,start=350)
plot.dmc(samples,start=360)

# Turn off migration and thin a bit
samples1 <- run.dmc(samples.dmc(nmc=100,samples=samples,thin=5), 
                   cores=12,report=5)

gelman.diag.dmc(samples1)
#      Point est. Upper C.I.
# A          1.06       1.09
# B.r1       1.07       1.10
# B.r2       1.06       1.09
# t0         1.06       1.10
# Wt         1.06       1.10
# Wd         1.06       1.09
# 
# Multivariate psrf
# 
# 1.09

effectiveSize.dmc(samples1)
#    A B.r1 B.r2   t0   Wt   Wd 
#  375  437  450  414  425  486 


plot.dmc(samples1,pll.chain=TRUE,density=TRUE)
plot.dmc(samples1)

pp <- post.predict.dmc(samples1)
plot.pp.dmc(pp,"cdf",layout=c(2,3))


# Recovery looks good
check.recovery.dmc(samples1,p.vector)
#                   A  B.r1  B.r2   t0   Wt  Wd
# True           1.00  1.00  2.00 0.20 1.00 0.1
# 2.5% Estimate  0.98  0.94  1.93 0.19 0.98 0.1
# 50% Estimate   1.03  0.98  1.97 0.20 1.01 0.1
# 97.5% Estimate 1.09  1.02  2.02 0.21 1.06 0.1
# Median-True    0.03 -0.02 -0.03 0.00 0.01 0.0

### Approach - Avoidance Factor

modelA <- model.dmc(p.map = list(A="1",B="1",t0="1",sd_v="1",st0="1",
                                Wt="ApAv",Wd="ApAv",mean_v="DT"), 
          factors=list(DT=DT,ApAv=c("Approach","Avoid")),
          constants = consts, 
          responses = c("r1","r2"),
          type="norm")
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"           "B"           "t0"          "Wt.Approach" "Wt.Avoid"    "Wd.Approach" "Wd.Avoid"   
# 
# Constants are (see attr(,"constants") ):
# mean_v.1112 mean_v.1121 mean_v.1122 mean_v.1221 mean_v.1222 mean_v.2121        sd_v         st0 
#        1112        1121        1122        1221        1222        2121           1           0 
# 
# Model type = norm (posdrift= TRUE ) 


p.vectorA <- c(A=1,B=1,t0=.2,Wt.Approach=1,Wt.Avoid=.5,Wd.Approach=.1,Wd.Avoid=.2)

# Simulate large data set
dataA <- data.model.dmc(simulate.dmc(p.vectorA,model1,n=1e4),modelA)
          
# Likelihoods look sensible
head(likelihood.dmc(p.vector1,data1))

# Profiles recover parameters
ylim=c(-50000,20000)
par(mfrow=c(2,4))
profile.dmc("A",.1,2,p.vector1,data1,ylim=ylim)
profile.dmc("B",.1, 2,p.vector1,data1,ylim=ylim)
profile.dmc("t0",.1,.3,p.vector1,data1,ylim=ylim)
profile.dmc("Wt.Approach",.1,2,p.vector1,data1,ylim=ylim)
profile.dmc("Wd.Approach",.01,.3,p.vector1,data1,ylim=ylim)
profile.dmc("Wt.Avoid",.1,2,p.vector1,data1,ylim=ylim)
profile.dmc("Wd.Avoid",.01,.3,p.vector1,data1,ylim=ylim)

### Do sampling

# Vague prior
p.priorA <- prior.p.dmc(
  dists = c("tnorm","tnorm","beta","tnorm","tnorm","tnorm","tnorm"),
  p1=c(A=.5,B=1.5,t0=1,Wt.Approach=.5,Wd.Approach=.5,Wt.Avoid=.5,Wd.Avoid=.5),                           
  p2=c(2,2,1,1,1,1,1),lower=c(0,0,.1,0,0,0,0),upper=c(NA,NA,1,NA,NA,NA,NA)
)
par(mfcol=c(2,4)); for (i in names(p.prior1)) plot.prior(i,p.prior1)

# Setup samples
samplesA <- samples.dmc(nmc=400,p.priorA,dataA)

# Run with migration
samplesA <- run.dmc(samplesA, report = 5, cores=12,p.migrate=.05)
samplesA <- run.dmc(samples.dmc(samples=samplesA,add=TRUE,nmc=100), report = 5, cores=12,p.migrate=.05)

plot.dmc(samplesA,pll.chain=TRUE)
# Good by 440
plot.dmc(samplesA,pll.chain=TRUE,start=450)

# Turn off migration and thin a bit
samplesA1 <- run.dmc(samples.dmc(nmc=100,samples=samplesA,thin=5), 
                   cores=12,report=5)


# Stuck chain pulls in by 50
plot.dmc(samplesA1,pll.chain=TRUE)
plot.dmc(samplesA1,layout=c(2,4))

gelman.diag.dmc(samplesA1)
#             Point est. Upper C.I.
# A                 1.03       1.05
# B                 1.01       1.02
# t0                1.02       1.03
# Wt.Approach       1.02       1.04
# Wt.Avoid          1.03       1.04
# Wd.Approach       1.04       1.06
# Wd.Avoid          1.04       1.07
# 
# Multivariate psrf
# 
# 1.08

effectiveSize.dmc(samplesA1)
#           A           B          t0 Wt.Approach    Wt.Avoid Wd.Approach    Wd.Avoid 
#         321         314         305         316         401         420         332 



ppA <- post.predict.dmc(samplesA1)
plot.pp.dmc(ppA,"cdf",layout=c(2,4))


# Recovery looks good
check.recovery.dmc(samplesA1,p.vector1)
#                   A    B  t0 Wt.Approach Wt.Avoid Wd.Approach Wd.Avoid
# True           1.00 1.00 0.2        1.00     0.50         0.1      0.2
# 2.5% Estimate  0.99 0.99 0.2        0.95     0.48         0.1      0.2
# 50% Estimate   1.01 1.00 0.2        0.98     0.49         0.1      0.2
# 97.5% Estimate 1.04 1.01 0.2        1.01     0.50         0.1      0.2
# Median-True    0.01 0.00 0.0       -0.02    -0.01         0.0      0.0


save_data (samples,samples1,pp,samplesA,samplesA1,ppA,file="T^-1_D.RData")


########################################################################
################ Parameter recovery study ##############################
########################################################################
# 40 subjects, 120 in each of approach and avoid.

pop.mean <- c(A=1, B=1,t0=.3, 
  Wt.Approach=1,Wt.Avoid=.5,Wd.Approach=.1,Wd.Avoid=.2)
pop.scale <-c(A=.1, B=.1, t0=.05,
  Wt.Approach=.1,Wt.Avoid=.05,Wd.Approach=.01,Wd.Avoid=.02)  
pop.prior <- prior.p.dmc(
  dists = rep("tnorm",7),
  p1=pop.mean,p2=pop.scale,
  lower=c(0,0,.1,0,0,0,0),upper=c(NA,NA,1,NA,NA,NA,NA)
)
 
##  Check population distributions
par(mfcol=c(2,4)); for (i in names(pop.prior)) plot.prior(i,pop.prior)


# Simulate some data, 12 cells so 240 in total per subject
raw.data <- h.simulate.dmc(modelA, p.prior = pop.prior, n = 20, ns = 40)
data.model <- data.model.dmc(raw.data, modelA)

# Take a look at parameters
ps <- attr(raw.data, "parameters")
round(apply(ps,2,mean),2)
#           A           B          t0 Wt.Approach    Wt.Avoid Wd.Approach    Wd.Avoid 
#        1.00        0.99        0.30        1.00        0.50        0.10        0.20 
round(apply(ps,2,sd),2)
#           A           B          t0 Wt.Approach    Wt.Avoid Wd.Approach    Wd.Avoid 
#        0.08        0.10        0.05        0.10        0.05        0.01        0.02 
       
### FIT FIXED EFFECTS

# specify a broader prior than the true population distribution
p.priorG <- prior.p.dmc(
  dists = rep("tnorm",7),
  p1=pop.mean,                           
  p2=pop.scale*5,
  lower=c(0,0,.1,0,0,0,0),upper=c(NA,NA,1,NA,NA,NA,NA)
)
##  Check population distributions
par(mfcol=c(2,4)); for (i in names(p.priorG)) plot.prior(i,p.priorG)

# Start with thin of 5
samplesF  <- h.samples.dmc(nmc = 50, p.priorG, data.model, thin = 5)

# The sampling was run from the command line from the file "recovery.R"
# save(raw.data,samplesF,file="recovery.RData")
samplesF  <- h.run.unstuck.dmc(samplesF, p.migrate = .05, cores = 40)
samplesF1 <- h.run.converge.dmc(h.samples.dmc(nmc=100, samples=samplesF), 
  nmc=50,cores=40)

load("recovery.RData")

# All chains are well converged. 
gelman.diag.dmc(samplesF1)
#   27   13    9   22   32    4   36   16   25    8   15   30   37    1    7   33   23 
# 1.07 1.07 1.08 1.08 1.08 1.08 1.08 1.08 1.09 1.09 1.09 1.09 1.09 1.09 1.09 1.09 1.09 
#   20    6   17   26   38   21   19   28   39   29   10    3   34   24   18   31   11 
# 1.09 1.09 1.09 1.09 1.09 1.09 1.09 1.09 1.09 1.09 1.10 1.10 1.10 1.10 1.10 1.10 1.10 
#   40    5    2   14   12   35 
# 1.10 1.10 1.10 1.10 1.10 1.10 
# Mean
# [1] 1.09

# Paramter chains look generally well converged and mixed, this just shows the
# first 
for (i in 1:1) {
  plot.dmc(samplesF1,density=FALSE,smooth=FALSE,subject=i,layout=c(2,4))
  plot.dmc(samplesF1[[1]],pll.chain=TRUE)
#  Sys.sleep(0.1) 
}

es <- effectiveSize.dmc(samplesF1)
round(apply(data.frame(es),1,mean))        
#           A           B          t0 Wt.Approach    Wt.Avoid Wd.Approach    Wd.Avoid 
#         400         409         414         377         374         411         409 
round(apply(data.frame(es),1,min))
#           A           B          t0 Wt.Approach    Wt.Avoid Wd.Approach    Wd.Avoid 
#         303         306         295         274         247         302         310 


# Good fits on average
ppF <- h.post.predict.dmc(samplesF1) 
plot.pp.dmc(ppF,layout=c(3,4))

# Here we use an lapply to get all fits, the assign to tmp is just to stop
# an empty list output. Most fits are pretty good, even subject 1. 
tmp <- lapply(ppF, function(x){plot.pp.dmc(x, style="cdf") })

# Check parameter recovery.
estF <- summary.dmc(samplesF1)  
# Calculates check.recovery.dmc for each subject and reports the average
h.check.recovery(samplesF1,ps)
#                     A     B    t0 Wt.Approach Wt.Avoid Wd.Approach Wd.Avoid
# True            0.979 1.000 0.297       1.007    0.501       0.097    0.199
# 2.5% Estimate   0.483 0.779 0.248       0.700    0.359       0.071    0.152
# 50% Estimate    0.976 1.024 0.298       1.099    0.495       0.100    0.199
# 97.5% Estimate  1.418 1.322 0.338       1.795    0.719       0.131    0.257
# Median-True    -0.003 0.023 0.001       0.092   -0.007       0.003   -0.001

save_data (samples,samples1,pp,samplesA,samplesA1,ppA,
  raw.data,samplesF,samplesF1,ppF,file="T^-1_D.RData")


### FIT RANDOM EFFECTS

# Use all truncated normal priors for locations
mu.prior <- prior.p.dmc(
  dists = rep("tnorm",7),
  p1=pop.mean,                           
  p2=c(2,2,1,2,2,.5,.5),
  lower=c(0,0,.1,0,0,0,0),upper=c(NA,NA,1,NA,NA,NA,NA)
)
par(mfcol=c(2,4)); for (i in names(mu.prior)) plot.prior(i,mu.prior)

# Uniform for SD
sigma.prior <- prior.p.dmc(
  dists = rep("beta", 7),
  c(A=1, B=1, t0=1,Wt.Approach=1,Wt.Avoid=1,Wd.Approach=1,Wd.Avoid=1),p2=rep(1,7)
)
par(mfcol=c(2,4)); for (i in names(sigma.prior)) plot.prior(i,sigma.prior)

# Make a hyper-prior list
pp.prior <- list(mu.prior, sigma.prior)

# The sampling was run from the command line from the file "recovery.R"
hsamples <- h.samples.dmc(nmc = 100, p.priorG, data.model, thin = 5, pp.prior = pp.prior)
# save(raw.data,samplesF,samplesF1,hsamples,file="recovery.RData")

hsamples <- h.run.unstuck.dmc(hsamples, cores=21, report=10, 
                      p.migrate = .05, h.p.migrate = .05)
hsamples1 <- h.run.converge.dmc(h.samples.dmc(nmc=100, samples=hsamples), 
  nmc=50,cores=21,report=10,minN=500)

# Get it back from external run
# load("recovery.RData")

# Chains look very nice, the initial decrease in pll is due to migration
plot.dmc(hsamples1,hyper=TRUE,pll.chain=TRUE)
# Also evident in some moves in parameters
plot.dmc(hsamples1,hyper=TRUE,layout=c(2,7))

hsamples2 <- h.samples.dmc(samples=hsamples1,nmc=0,remove=1:200,add=TRUE)

h.gelman.diag.dmc(hsamples2)
# hyper     7    24    27    29     6     3    14    10     1    23    32     8 
#  1.01  1.01  1.01  1.01  1.01  1.01  1.01  1.01  1.01  1.01  1.01  1.01  1.01 
#    34     2    17    39    37    38    40    15    33    36    19    20    22 
#  1.01  1.01  1.01  1.01  1.01  1.01  1.01  1.01  1.01  1.01  1.01  1.01  1.01 
#    26    31    16    11    35     4    13    25    28    12     9    30     5 
#  1.01  1.01  1.02  1.02  1.02  1.02  1.02  1.02  1.02  1.02  1.02  1.02  1.02 
#    21    18 
#  1.02  1.02 


# OK yield at hyper given nominal 21*100=2100
effectiveSize.dmc(hsamples2,hyper=TRUE)
#           A.h1           B.h1          t0.h1 Wt.Approach.h1    Wt.Avoid.h1 
#           2484           3298           3727           2576           2578 
# Wd.Approach.h1    Wd.Avoid.h1           A.h2           B.h2          t0.h2 
#           3219           3120           1806           2955           3653 
# Wt.Approach.h2    Wt.Avoid.h2 Wd.Approach.h2    Wd.Avoid.h2 
#           2076           2668           3195           2950 
              

# Fits are failry nice.
hpp <- h.post.predict.dmc(hsamples2)
plot.pp.dmc(hpp,layout=c(3,4))

# Parameter recovery

# First check population means (ptype=1 picks out means), failry good, coverage
# acheived for all parameters
h.check.recovery(hsamples2,ps=pop.mean,hyper=TRUE,ptype=1)
#                    A      B     t0 Wt.Approach Wt.Avoid Wd.Approach Wd.Avoid
# True           1.000  1.000  0.300       1.000    0.500       0.100    0.200
# 2.5% Estimate  0.967  0.946  0.276       0.992    0.462       0.097    0.191
# 50% Estimate   1.004  0.991  0.294       1.061    0.489       0.102    0.197
# 97.5% Estimate 1.040  1.035  0.311       1.135    0.516       0.107    0.203
# Median-True    0.004 -0.009 -0.006       0.061   -0.011       0.002   -0.003

# ptype=2 picks out scales, not so good for the Wt parameters
h.check.recovery(hsamples2,ps=pop.scale,hyper=TRUE,ptype=2)
#                     A     B    t0 Wt.Approach Wt.Avoid Wd.Approach Wd.Avoid
# True            0.100 0.100 0.050       0.100    0.050       0.010    0.020
# 2.5% Estimate   0.064 0.083 0.044       0.125    0.050       0.010    0.012
# 50% Estimate    0.088 0.110 0.055       0.172    0.067       0.014    0.016
# 97.5% Estimate  0.121 0.149 0.071       0.235    0.090       0.018    0.021
# Median-True    -0.012 0.010 0.005       0.072    0.017       0.004   -0.004

# However, failry good at individual participant level
h.check.recovery(hsamples2,ps)
#                    A     B     t0 Wt.Approach Wt.Avoid Wd.Approach Wd.Avoid
# True           1.001 0.989  0.295       1.004    0.504       0.097    0.200
# 2.5% Estimate  0.851 0.837  0.264       0.798    0.398       0.084    0.173
# 50% Estimate   1.004 0.990  0.294       1.054    0.486       0.102    0.197
# 97.5% Estimate 1.155 1.153  0.321       1.366    0.592       0.119    0.220
# Median-True    0.003 0.001 -0.001       0.049   -0.017       0.004   -0.004

save_data (raw.data,samples,samples1,pp,samplesA,samplesA1,ppA,
  data,samplesF,samplesF1,ppF,hsamples,hsamples1,hsamples2,hpp,file="T^-1_D.RData")





