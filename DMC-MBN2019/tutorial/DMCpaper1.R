rm(list=ls()) 

# Current working directory must be set to the top-level folder  
# containing the dmc subfolder 
source ("dmc/dmc.R")

# Let's work with the LBA model
load_model ("LBA","lba_B.R")

# # Load precomputed results of time-consuming computations. 
# load_data("DMCpaper1.RData")

# Note that in the examples below we often use large numbers of cores. If you
# are on a system with a small number of cores this can cause it to freeze.
# If you run these commands adjust the number of cores to suit your system.
 
### Example 1.0: Non-identified model ----

# Set up an LBA model for a binary choice in the simplest possible choice design 
# with two stimuli (s1 and s2) mapping on two responses (r1 and r2). We assume
# there is no non-decision time noise (the standard LBA assumption) and estimate 
# a single value for every LBA parameter except mean_v, where we allow one value
# for responses that match the stimulus (i.e., r1 to s1 and r2 to s2) and a 
# second smaller value for mismatches in order to allow for above chance 
# accuracy. 

modelLBA0 <- model.dmc(
  p.map = list(A="1",B="1",mean_v="M",sd_v="1",t0="1",st0="1"),
  constants=c(st0=0),                 
  match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2")),
  responses=c("r1","r2"),
  type="norm")

# Select some parameter that result in RTs of about half a second and accuracy
# around 75%. This is the same for both stimuli as they have identical parameters.
pLBA0  <- c(A=.5,B=1,sd_v=1,mean_v.true=4,mean_v.false=3,t0=.2)

# Create a large amount of data (10,000 responses for each stimuli) so we can
# look at the large-sample behaviour of the model (i.e., behaviour that will 
# change little if you run the example again yourself). Note that response time
# in DMC is assumed to be on the seconds scale.
data <- simulate.dmc(pLBA0,modelLBA0,n=1e4)

# Take a look at the simulated data. Note that the commented-out code saves a 
# pdf of the plot. The appearance options change the defaults for a journal 
# figure without colour (the defaults use colour).

# pdf("DMCpaper/LBA0density.pdf",height=8,width=10)
par(mfrow=c(1,2))
plot.cell.density(data.cell=data[data$S=="s1",],
                  C="r1",xlim=c(0,1),main="Choice s1",
                  CorrectError.col=c("black","black"),CorrectError.lty=c(1,3))
plot.cell.density(data.cell=data[data$S=="s2",],
                  C="r2",xlim=c(0,1),main="Choice s2",
                  CorrectError.col=c("black","black"),CorrectError.lty=c(1,3))
# dev.off()

# Bind the data and model together ready for sampling.
dmLBA0 <- data.model.dmc(data,modelLBA0)

# Specify fairly uninformative normal priors, truncated below at zero for A and 
# B as they must be positive, unbounded for the v, bounded below by 0.1 for t0
# (on the assumption that a choice RT less than 0.1s is implausible).
p.priorLBA0 <- prior.p.dmc(
  p1=c(A=1,B=1,sd_v=1,mean_v.true=2,mean_v.false=1,t0=.3),                           
  p2=c(1,1,1,3,3,.25),lower=c(0,0,0,NA,NA,.1),upper=rep(NA,6)
)

# Plot the priors to check they look reasonable.
par(mfcol=c(2,3)); for (i in names(p.priorLBA0)) plot.prior(i,p.priorLBA0)

# Combine the prior with the data and model to create an object ready to contain
# the results of sampling. It also gets samples for the first iteration by 
# sampling from the prior (see tutorial dmc_3_1_LNRsampling.R for other options)
# until it gets values that have well defined likelihoods.

# In the example below, nmc is the number of iterations the object is ready
# to accept. This can change during sampling when RUN.dmc is used to do sampling
# see below, and when using RUN.dmc nmc is best set to 120 in the first instance. We
# also specify thinning of 10 so that we store less redundant information (this
# value was chosen from previous experience with this model, see tutorial 
# dmc_3_2_LNRassesing.R for details on thinning). 
sLBA0 <- samples.dmc(nmc=120,p.priorLBA0,dmLBA0,thin=10)


# Now run the sampling using RUN.dmc, which tries to deal with any convergence 
# problems automatically. This is an extension of the automatic convergence
# algorithms described in tutorial dmc_3_3_LBA.R, run.unstuck.dmc, which 
# deals with chains that become stuck in bad areas, and run.converge.dmc, which
# can add new samples and remove old samples in an attempt to obtain good chain 
# mixing. RUN.dmc also checks to see if the samples are stationary (not
# changing in location and spread), by comparing the first and last third of
# the chains, throwing away the first third if the checks fail and adding the
# same number of new samples. It does all checks repeatedly, starting
# with stuck chains then moving on to mixing and stationarity, but going back to
# to stuck chains if necessary until it obtains both convergence and, if 
# requested, a minimum number of samples. 

# Where verbose=TRUE, as here, the test results are reported. Here we also 
# request a minimum sample size of 2700 (when use.effectiveSize = TRUE the 
# size estimate also takes account of autocorrelation between samples, see
# tutorial dmc_3_2_LNRassesing.R). We use 18 cores, one for each chain (by 
# default there are three times as many chains as parameters), for speed.
sLBA0 <- RUN.dmc(sLBA0,cores=18,verbose=TRUE,
                 minN=2700,use.effectiveSize = FALSE)

# Check chain mixing, using Brooks & Gelman (1998) "looks good R hat" measure
# (see the coda packages ?gelman.diag help and tutorial dmc_3_2_LNRassesing.R)
gelman.diag.dmc(sLBA0)
# Potential scale reduction factors:
# 
#              Point est. Upper C.I.
# A                  1.02       1.03
# B                  1.01       1.02
# mean_v.true        1.02       1.03
# mean_v.false       1.02       1.03
# sd_v               1.02       1.03
# t0                 1.01       1.01
# 
# Multivariate psrf
# 
# 1.03

# Posterior log-likelihoods are flat indicating convergence.
plot.dmc(sLBA0,pll.chain=TRUE)


# The parameter chains are also flat fat hairy caterpillars. 
plot.dmc(sLBA0,layout=c(2,3))

# However, (apart from non-decision time) the median parameter estimates are off 
# by a long way, and they are very uncertain (wide credible intervals).
recover0 <- check.recovery.dmc(sLBA0,pLBA0)
#                   A    B mean_v.true mean_v.false sd_v   t0
# True           0.50 1.00        4.00         3.00 1.00 0.20
# 2.5% Estimate  0.31 0.70        2.77         2.05 0.69 0.18
# 50% Estimate   0.58 1.23        4.85         3.61 1.22 0.20
# 97.5% Estimate 0.91 1.88        7.39         5.48 1.86 0.21
# Median-True    0.08 0.23        0.85         0.61 0.22 0.00

# This is because we have failed to fix a scaling parameter! It is evident if
# we look at the correlations among the different parameters in the posterior.
# The thin=10 argument plots every 10th recorded sample to speed things up; in the paper 
# it was set to 1 to give a more detailed picture but the correlations don’t 
# change much.
pairs.dmc(sample=sLBA0,thin=10)   

### Example 1.1: Conventional model ----

# Now let’s fit the conventional identified model where we fix sd_v, here at the
# true value using the constants argument. This takes sd_v out of sampling but
# uses the value specified in constants when calculating the likelihood.
modelLBA1 <- model.dmc(
  p.map = list(A="1",B="1",mean_v="M",sd_v="1",t0="1",st0="1"),
  constants=c(st0=0,sd_v=1),
  match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2")),
  responses=c("r1","r2"),
  type="norm")

# We use exactly the same data as before.
dmLBA1 <- data.model.dmc(data,modelLBA1)

# We also use the same prior, except for sd_v which is no longer sampled.
sLBA1 <- samples.dmc(nmc=120,p.priorLBA0[-3],dmLBA1,thin=10)

# Run the sampling
sLBA1 <- RUN.dmc(sLBA1,cores=15,minN=2700,use.effectiveSize = FALSE)

# Again converged well.
gelman.diag.dmc(sLBA1)
plot.dmc(sLBA1,pll.chain=TRUE)
plot.dmc(sLBA1,layout=c(2,3))

# Now recovery is extremely good and uncertainty is very small.
pLBA1 <- pLBA0[-3] # remove sd_v
recover1 <- check.recovery.dmc(sLBA1,pLBA1)
#                    A    B mean_v.true mean_v.false   t0
# True            0.50 1.00        4.00         3.00 0.20
# 2.5% Estimate   0.37 0.91        3.87         2.84 0.18
# 50% Estimate    0.48 1.01        3.99         2.97 0.20
# 97.5% Estimate  0.55 1.17        4.11         3.10 0.21
# Median-True    -0.02 0.01       -0.01        -0.03 0.00


# There is still quite some correlation, although not as extreme as before, 
# except for B and t0, which are now much more highly correlated. 
pairs.dmc(sLBA1,thin=10)  


### Example 1.2: Conventional model with different fixed value ----

# Double the value of sd_v
modelLBA2 <- model.dmc(
  p.map = list(A="1",B="1",mean_v="M",sd_v="1",t0="1",st0="1"),
  constants=c(st0=0,sd_v=2),
  match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2")),
  responses=c("r1","r2"),
  type="norm")

# We use exactly the same data as before.
dmLBA2 <- data.model.dmc(data,modelLBA2)

# We also use the same prior, except for sd_v which is no longer sampled.
sLBA2 <- samples.dmc(nmc=120,p.priorLBA0[-3],dmLBA2,thin=10)
sLBA2 <- RUN.dmc(sLBA2,cores=15,minN=2700,use.effectiveSize = FALSE)


# Again converged well.
gelman.diag.dmc(sLBA2)
plot.dmc(sLBA2,pll.chain=TRUE)
plot.dmc(sLBA2,layout=c(2,3))

# The recovered accumulator parameters (all but t0) just double to account for
# the doubling of sd_v.
pLBA2 <- c(pLBA1[-5]*2,pLBA1[5]) # double accumulator parameters
check.recovery.dmc(sLBA2,pLBA2)
#                    A    B mean_v.true mean_v.false   t0
# True            1.00 2.00        8.00         6.00 0.20
# 2.5% Estimate   0.76 1.79        7.71         5.65 0.18
# 50% Estimate    0.96 2.00        7.94         5.91 0.20
# 97.5% Estimate  1.10 2.28        8.18         6.17 0.21
# Median-True    -0.04 0.00       -0.06        -0.09 0.00


### Example 1.3: Fixing other parameters ----

# Fix B as reference
modelLBA3a <- model.dmc(
  p.map = list(A="1",B="1",mean_v="M",sd_v="1",t0="1",st0="1"),
  constants=c(st0=0,B=1),
  match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2")),
  responses=c("r1","r2"),
  type="norm")

# We use exactly the same data as before.
dmLBA3a <- data.model.dmc(data,modelLBA3a)
sLBA3a <- samples.dmc(nmc=120,p.priorLBA0[-2],dmLBA3a,thin=10)
sLBA3a <- RUN.dmc(sLBA3a,cores=15,minN=2700,use.effectiveSize = FALSE)

# Again converged well.
gelman.diag.dmc(sLBA3a)
plot.dmc(sLBA3a,pll.chain=TRUE)
plot.dmc(sLBA3a,layout=c(2,3))

# Generally greater bias, except A and t0, and uncertainty except t0
recover3a <- check.recovery.dmc(sLBA3a,pLBA0[-2])
#                    A mean_v.true mean_v.false sd_v   t0
# True            0.50        4.00         3.00 1.00 0.20
# 2.5% Estimate   0.36        3.58         2.70 0.89 0.18
# 50% Estimate    0.48        3.98         2.96 1.00 0.20
# 97.5% Estimate  0.60        4.36         3.19 1.11 0.21
# Median-True    -0.02       -0.02        -0.04 0.00 0.00

# % Absolute Bias: fix B vs. fix sd_v 
round(100*(abs(apply(recover3a,2,function(x){abs(x[5])})[-4])-
             abs(apply(recover1,2,function(x){abs(x[5])})[-2]))/
        recover1[1,-2],1)
#            A  mean_v.true mean_v.false           t0 
#         -1.5          0.2          0.3         -1.1 

# Uncertainty Ratio: fix B vs. fix sd_v         
round(apply(recover3a,2,function(x){x[4]-x[2]})[-4]/
        apply(recover1,2,function(x){x[4]-x[2]})[-2],1)
#            A  mean_v.true mean_v.false           t0 
#          1.3          3.2          1.9          0.9 

# A pairs plot shows this is due to much higher parameter correlations.
pairs.dmc(sLBA3a,thin=10) #Error


# Fix mean_v.true as reference
modelLBA3b <- model.dmc(
  p.map = list(A="1",B="1",mean_v="M",sd_v="1",t0="1",st0="1"),
  constants=c(st0=0,mean_v.true=4),
  match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2")),
  responses=c("r1","r2"),
  type="norm")

# We use exactly the same data as before.
dmLBA3b <- data.model.dmc(data,modelLBA3b)
sLBA3b <- samples.dmc(nmc=120,p.priorLBA0[-4],dmLBA3b,thin=10)
sLBA3b <- RUN.dmc(sLBA3b,cores=15,minN=2700,use.effectiveSize = FALSE)

# Again converged well.
gelman.diag.dmc(sLBA3b)
plot.dmc(sLBA3b,pll.chain=TRUE)
plot.dmc(sLBA3b,layout=c(2,3))

# A little better than fixing sd_v
recover3b <- check.recovery.dmc(sLBA3b,pLBA0[-4])
#                    A    B mean_v.false sd_v   t0
# True            0.50 1.00         3.00 1.00 0.20
# 2.5% Estimate   0.36 0.92         2.92 0.97 0.18
# 50% Estimate    0.48 1.01         2.98 1.00 0.20
# 97.5% Estimate  0.55 1.16         3.04 1.04 0.21
# Median-True    -0.02 0.01        -0.02 0.00 0.00

# %Bias: fix mean_v.true vs. fix sd_v 
round(100*(abs(apply(recover3b,2,function(x){abs(x[5])})[-4])-
             abs(apply(recover1,2,function(x){abs(x[5])})[-3]))/
        recover1[1,-3],1)
#            A            B mean_v.false           t0 
#         -0.7          0.1         -0.2         -0.2 

# Uncertainty ratio: fix mean_v.true vs. fix sd_v         
round(apply(recover3b,2,function(x){x[4]-x[2]})[-4]/
        apply(recover1,2,function(x){x[4]-x[2]})[-3],1)
#            A            B mean_v.false           t0 
#          1.0          0.9          0.5          1.0

# A high correlation between B and t0 returns
pairs.dmc(sLBA3b,thin=10)

# Fix mean_v.false as reference
modelLBA3c <- model.dmc(
  p.map = list(A="1",B="1",mean_v="M",sd_v="1",t0="1",st0="1"),
  constants=c(st0=0,mean_v.false=3),
  match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2")),
  responses=c("r1","r2"),
  type="norm")

# We use exactly the same data as before.
dmLBA3c <- data.model.dmc(data,modelLBA3c)
sLBA3c <- samples.dmc(nmc=120,p.priorLBA0[-5],dmLBA3c,thin=10)
sLBA3c <- RUN.dmc(sLBA3c,cores=15,minN=2700,use.effectiveSize = FALSE)

# Again converged well.
gelman.diag.dmc(sLBA3c)
plot.dmc(sLBA3c,pll.chain=TRUE)
plot.dmc(sLBA3c,layout=c(2,3))


recover3c <- check.recovery.dmc(sLBA3c,pLBA0[-5])
#                    A    B mean_v.true sd_v   t0
# True            0.50 1.00        4.00 1.00 0.20
# 2.5% Estimate   0.37 0.94        3.95 0.97 0.18
# 50% Estimate    0.48 1.02        4.03 1.01 0.20
# 97.5% Estimate  0.56 1.14        4.10 1.06 0.21
# Median-True    -0.02 0.02        0.03 0.01 0.00

# % Bias: fix mean_v.false vs. fix sd_v 
round(100*(abs(apply(recover3c,2,function(x){abs(x[5])})[-4])-
             abs(apply(recover1,2,function(x){abs(x[5])})[-4]))/
        recover1[1,-4],1)
#           A           B mean_v.true          t0 
#        -0.9         1.1         0.4        -0.2 

# Uncertainty ratio: fix mean_v.false vs. fix sd_v         
round(apply(recover3c,2,function(x){x[4]-x[2]})[-4]/
        apply(recover1,2,function(x){x[4]-x[2]})[-4],1)
#           A           B mean_v.true          t0 
#         1.1         0.7         0.6         1.0

# A pairs plot shows this is due to much higher parameter correlations.
pairs.dmc(sLBA3c,thin=10)

# Fix A as reference
modelLBA3d <- model.dmc(
  p.map = list(A="1",B="1",mean_v="M",sd_v="1",t0="1",st0="1"),
  constants=c(st0=0,A=0.5),
  match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2")),
  responses=c("r1","r2"),
  type="norm")

# We use exactly the same data as before.
dmLBA3d <- data.model.dmc(data,modelLBA3d)
sLBA3d <- samples.dmc(nmc=120,p.priorLBA0[-1],dmLBA3d,thin=10)
sLBA3d <- RUN.dmc(sLBA3d,cores=15,minN=2700,use.effectiveSize = FALSE)


# Again converged well.
gelman.diag.dmc(sLBA3d)
plot.dmc(sLBA3d,pll.chain=TRUE)
plot.dmc(sLBA3d,layout=c(2,3))

# Posterior medians are accurate but uncertainty is greater than in the 
# conventional case for accumulator parameters.

recover3d <- check.recovery.dmc(sLBA3d,pLBA0[-1])
#                   B mean_v.true mean_v.false sd_v   t0
# True           1.00        4.00         3.00 1.00 0.20
# 2.5% Estimate  0.87        3.69         2.73 0.92 0.18
# 50% Estimate   1.08        4.22         3.14 1.06 0.20
# 97.5% Estimate 1.39        5.03         3.79 1.25 0.21
# Median-True    0.08        0.22         0.14 0.06 0.00

# Much worse than fixing sd_v
round(100*(abs(apply(recover3d,2,function(x){abs(x[5])})[-4])-
             abs(apply(recover1,2,function(x){abs(x[5])})[-1]))/
        recover1[1,-1],1)
#            B  mean_v.true mean_v.false           t0 
#          6.7          5.1          3.8          0.3
round(apply(recover3d,2,function(x){x[4]-x[2]})[-4]/
        apply(recover1,2,function(x){x[4]-x[2]})[-1],1)
#            B  mean_v.true mean_v.false           t0 
#          2.0          5.6          4.1          0.8

# And much worse than fixing B.
round(100*(abs(apply(recover3d,2,function(x){abs(x[5])})[-1])-
             abs(apply(recover3a,2,function(x){abs(x[5])})[-1]))/
        recover3a[1,-1],1)
#   mean_v.true mean_v.false         sd_v           t0 
#          4.9          3.4          5.6          1.4
round(apply(recover3d,2,function(x){x[4]-x[2]})[-1]/
        apply(recover3a,2,function(x){x[4]-x[2]})[-1],1)
#  mean_v.true mean_v.false         sd_v           t0 
#          1.7          2.2          1.5          1.0 

# A pairs plot shows this is due to much higher parameter correlations.
pairs.dmc(sLBA3d,thin=10)

### Example 1.4: Fixing only one sd_v ----

# Allow sd_v to vary with match, fixing only the value for the mismatch
modelLBA4 <- model.dmc(
  p.map = list(A="1",B="1",mean_v="M",sd_v="M",t0="1",st0="1"),
  constants=c(st0=0,sd_v.false=1),
  match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2")),
  responses=c("r1","r2"),
  type="norm")


# Set sd_v.true to less than the fixed value sd_v.false=1 as is often seen
# in data (although usually a bit larger than the value of 0.5 used here)
pLBA4  <- c(A=.5,B=1,sd_v.true=0.5,mean_v.true=4,mean_v.false=3,t0=.2)
data4 <- simulate.dmc(pLBA4,modelLBA4,n=1e4)

# Note that there is a slight increase in accuracy as the lower 
par(mfrow=c(1,2))
plot.cell.density(data.cell=data4[data4$S=="s1",],
                  C="r1",xlim=c(0,1),main="Choice s1",lwd=2)
plot.cell.density(data.cell=data4[data4$S=="s2",],
                  C="r2",xlim=c(0,1),main="Choice s2",lwd=2)

# Note that the greater variability of the mismatch distribution tends to speed
# up errors relative to corrects.

# score the two data sets with sd_v the same for match and mismatch (data)
# and mismatch > match (data4)
crct <- as.numeric(data$S)==as.numeric(data$R)
crct4 <- as.numeric(data4$S)==as.numeric(data4$R)

# Errors slower than correct on average.
round(tapply(data$RT,crct,mean),3)
# FALSE  TRUE 
# 0.527 0.502

# With sd_v larger for mismatch than match errors are now faster on average.
round(tapply(data4$RT,crct4,mean),3)
# FALSE  TRUE 
# 0.496 0.507


# Make a prior (the same as LBA0 but with a different name for the sampled sd_v 
# parameter, and sample the new model.
p.priorLBA4 <- prior.p.dmc(
  p1=c(A=1,B=1,sd_v.true=1,mean_v.true=2,mean_v.false=1,t0=.3),                           
  p2=c(1,1,1,3,3,.25),lower=c(0,0,0,NA,NA,.1),upper=rep(NA,6)
)
# par(mfcol=c(2,3)); for (i in names(p.priorLBA4)) plot.prior(i,p.priorLBA4)

dmLBA4 <- data.model.dmc(data4,modelLBA4)
sLBA4 <- samples.dmc(nmc=120,p.priorLBA4,dmLBA4,thin=10)
sLBA4 <- RUN.dmc(sLBA4,cores=18,minN=2700,use.effectiveSize = FALSE)


# Again converged well.
gelman.diag.dmc(sLBA4)
plot.dmc(sLBA4,pll.chain=TRUE)
plot.dmc(sLBA4,layout=c(2,3))

# Parameter recovery is comparable to the simpler standard model
recover4 <- check.recovery.dmc(sLBA4,pLBA4)
#                   A     B mean_v.true mean_v.false sd_v.true   t0
# True           0.50  1.00        4.00         3.00      0.50 0.20
# 2.5% Estimate  0.46  0.84        3.66         2.67      0.47 0.19
# 50% Estimate   0.50  0.94        3.89         2.90      0.49 0.21
# 97.5% Estimate 0.53  1.06        4.11         3.11      0.51 0.22
# Median-True    0.00 -0.06       -0.11        -0.10     -0.01 0.01


### Example 1.5: Small sample parameter-recovery study ----

# Let's use DMC's ability to simulate and fit multiple subjects to
# do a parameter-recovery study 

# First simulate 200 subjects (i.e., replications), each with identical parameters but with
# different random realizations of data (see dmc_4_1_Simulating.R for 
# how to use this function) 
data5 <- h.simulate.dmc(modelLBA4,ps=pLBA4,ns=200,n=1e2)

# The data-model object is now a list, with one subject in each slot.
dLBA5 <- data.model.dmc(data5,modelLBA4)

# Set up object ready for sampling (see dmc_4_4_LNRfixed for more details).
sLBA5 <- h.samples.dmc(nmc=120,p.priorLBA4,dLBA5,thin=10)


# Run the sampling, h.RUN.dmc just applies RUN.dmc to each subject. By
# default it uses one core per subject, which is more efficient than spreading
# chains across cores.
sLBA5 <- h.RUN.dmc(sLBA5,cores=60,meanN=2700,use.effectiveSize = FALSE)

# h.check.recovery applies check.recovery to each fit then takes the average. It
# also calculates "coverage", the percentage of fits for which the true value 
# falls in the 95% credible interval. If uncertainty is being properly 
# quantified this should be around 95%.

# Note this can be a little slow to run with 200 fits! The average recovery
# is quite good, indicating no tendency to the development of small-sample 
# biases, but as expected, average uncertainty is quite high. 
recover5 <- h.check.recovery.dmc(sLBA5,pLBA4)
#                     A      B mean_v.true mean_v.false sd_v.true     t0
# True            0.500  1.000       4.000        3.000     0.500  0.200
# 2.5% Estimate   0.197  0.506       2.905        1.838     0.362  0.133
# 50% Estimate    0.500  0.989       4.078        3.030     0.525  0.212
# 97.5% Estimate  0.814  1.616       5.487        4.378     0.723  0.287
# Median-True     0.000 -0.011       0.078        0.030     0.025  0.012
# Coverage       93.000 95.000      94.500       94.500    90.000 93.000

# Small increase in % bias
round(100*(abs(apply(recover5,2,function(x){abs(x[3])}))-
             abs(apply(recover4,2,function(x){abs(x[3])})))/
        recover4[1,],1)
#            A            B  mean_v.true mean_v.false    sd_v.true           t0 
#          1.0          4.8          4.7          4.5          7.1          2.0 
#          
# Bias around 1-5%
round(100*recover5[5,]/recover5[1,],1)
#            A            B  mean_v.true mean_v.false    sd_v.true           t0 
#          0.1         -1.1          1.9          1.0          5.0          5.8 

# The gain in precision of the 95% CI for the increase in sample size by a 
# factor of 100 (i.e., 20,000 vs. 200) varies from close to a square-root law 
# (sqrt(100)=10) to about half that value. 
round(apply(recover5,2,function(x){x[4]-x[2]})/
        apply(recover4,2,function(x){x[4]-x[2]}),1)
#            A            B  mean_v.true mean_v.false    sd_v.true           t0 
#          9.5          5.0          5.8          5.8          9.8          5.1 

# Assessing posterior correlations is a little more difficult here as we need to
# somehow aggregate over subjects to produce a single plot (by default DMC
# produces one plot per subjects. By default, DMC does this by first 
# standardizing each subjects' parameters (i.e., subtracting the mean and 
# dividing by the standard deviation) so correlations are not deflated by 
# individual differences, but this also means that the natural scale of the 
# parameters is lost in the plots. Here where all subjects have identical
# data generating distributions this is not necessary so we turn off the 
# default to retain the scale. NOTE this can be very slow so we turn thinning
# up high, although this still means there are 5400 points in each panel of the
# plot!

pairs.dmc(sLBA5,thin=100, scale.subjects = FALSE)

# Note that the spread of the t0 estimate now become sufficient that its
# distribution is clearly truncated by the prior lower bound of 0.1s. In real
# small-sample data such low estimates of t0 are common. 

# Compare updating between small and large-sample cases
plot.dmc(sLBA4,layout=c(2,3),p.prior=p.priorLBA4)
plot.dmc(sLBA5[[2]],layout=c(2,3),p.prior=p.priorLBA4)


### Example 1.6: Simple RT/no errors ----

# Setup a model with constants so that the incorrect response can never win.
modelLBA6 <- model.dmc(
  p.map = list(A="1",B="1",mean_v="M",sd_v="M",t0="1",st0="1"),
  constants=c(st0=0,sd_v.true=1,mean_v.false=-100,sd_v.false=1e-3),
  match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2")),
  responses=c("r1","r2"),
  type="norm")

pLBA6  <- c(A=.5,B=1,mean_v.true=4,t0=.2)
data6 <- simulate.dmc(pLBA6,modelLBA6,n=1e4)
dmLBA6 <- data.model.dmc(data6,modelLBA6)

# Only one response to each stimulus 
par(mfrow=c(1,2))
plot.cell.density(data.cell=data6[data6$S=="s1",],xlim=c(0,1),
                  main="Simple s1")
plot.cell.density(data.cell=data6[data6$S=="s2",],
                  xlim=c(0,1),main="Choice s2")


# Make prior and sample
p.priorLBA6 <- prior.p.dmc(
  p1=c(A=1,B=1,mean_v.true=2,t0=.3),                           
  p2=c(1,1,3,.25),lower=c(0,0,NA,.1),upper=rep(NA,4)
)
sLBA6 <- samples.dmc(nmc=120,p.priorLBA6,dmLBA6,thin=10)
sLBA6 <- RUN.dmc(sLBA6,cores=12,minN=2700,use.effectiveSize = FALSE)


# Again converged well.
gelman.diag.dmc(sLBA6)
plot.dmc(sLBA6,pll.chain=TRUE)
plot.dmc(sLBA6,layout=c(2,2))

# Recovery is similar to the choice case but with a bigger miss in B and t0.
# Levels of uncertainty are similar.
check.recovery.dmc(sLBA6,pLBA6)
#                   A     B mean_v.true   t0
# True           0.50  1.00        4.00 0.20
# 2.5% Estimate  0.48  0.78        3.86 0.19
# 50% Estimate   0.56  0.89        3.97 0.22
# 97.5% Estimate 0.62  1.05        4.09 0.23
# Median-True    0.06 -0.11       -0.03 0.02

# As in the choice case there is a strong correlation between B and t0 
pairs.dmc(sLBA6,thin=10)

### Example 1.7: Simple RT with LATER ----

# Setup a model with constants so that the incorrect response can never win.
modelLBA7 <- model.dmc(
  p.map = list(A="1",B="1",mean_v="M",sd_v="M",t0="1",st0="1"),
  constants=c(st0=0,sd_v.true=1,mean_v.false=-100,sd_v.false=1e-3,A=0),
  match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2")),
  responses=c("r1","r2"),
  type="norm")

pLBA7  <- c(B=1,mean_v.true=4,t0=.2)
data7 <- simulate.dmc(pLBA7,modelLBA7,n=1e4)
dmLBA7 <- data.model.dmc(data7,modelLBA7)


# Make prior and sample
p.priorLBA7 <- prior.p.dmc(
  p1=c(B=1,mean_v.true=2,t0=.3),                           
  p2=c(1,3,.25),lower=c(0,NA,.1),upper=rep(NA,3)
)
sLBA7 <- samples.dmc(nmc=120,p.priorLBA7,dmLBA7,thin=10)
sLBA7 <- RUN.dmc(sLBA7,cores=9,minN=2700,use.effectiveSize = FALSE)


# Again converged well.
gelman.diag.dmc(sLBA7)
plot.dmc(sLBA7,pll.chain=TRUE)
plot.dmc(sLBA7,layout=c(2,2))

# Similar recovery to LBA 
check.recovery.dmc(sLBA7,pLBA7,digits=3)
#                    B mean_v.true     t0
# True           1.000       4.000  0.200
# 2.5% Estimate  1.012       4.021  0.187
# 50% Estimate   1.060       4.114  0.193
# 97.5% Estimate 1.110       4.213  0.198
# Median-True    0.060       0.114 -0.007

# Even stronger correlations than the LBA
pairs.dmc(sLBA7,thin=10)

### Example 1.8 Standard model with low error rate; asymptotic ----

# Use same model as Example 4 but reduce mean_v.false so error rate is ~ 2.5%

pLBA8  <- c(A=.5,B=1,sd_v.true=0.5,mean_v.true=4,mean_v.false=1.5,t0=.2)
data8 <- simulate.dmc(pLBA8,modelLBA4,n=1e4)

# As expected the error rate is an order of magnitude smaller than Example 4
par(mfrow=c(1,2))
plot.cell.density(data.cell=data8[data8$S=="s1",],
                  C="r1",xlim=c(0,1),main="Choice s1",lwd=2)
plot.cell.density(data.cell=data8[data8$S=="s2",],
                  C="r2",xlim=c(0,1),main="Choice s2",lwd=2)
mean(as.numeric(data8$S)==as.numeric(data8$R))
# [1] 0.97635

# Note that with the lower error rate errors become slow.
crct8 <- as.numeric(data8$S)==as.numeric(data8$R)
round(tapply(data8$RT,crct8,mean),3)
# FALSE  TRUE 
# 0.534 0.516

dmLBA8 <- data.model.dmc(data8,modelLBA4)
sLBA8 <- samples.dmc(nmc=120,p.priorLBA4,dmLBA8,thin=10)
sLBA8 <- RUN.dmc(sLBA8,cores=18,minN=2700,use.effectiveSize = FALSE)


# Again converged well.
gelman.diag.dmc(sLBA8)
plot.dmc(sLBA8,pll.chain=TRUE)
plot.dmc(sLBA8,layout=c(2,3))

# Parameter recovery. 
recover8 <- check.recovery.dmc(sLBA8,pLBA8)
#                    A     B mean_v.true mean_v.false sd_v.true   t0
# True            0.50  1.00        4.00         1.50      0.50 0.20
# 2.5% Estimate   0.44  0.84        3.53         1.06      0.45 0.18
# 50% Estimate    0.48  0.96        3.85         1.38      0.49 0.20
# 97.5% Estimate  0.52  1.09        4.17         1.67      0.53 0.22
# Median-True    -0.02 -0.04       -0.15        -0.12     -0.01 0.00

# Bias is not much affected relative to the 25% error case, sometimes more, 
# sometimes less.
round(100*(abs(apply(recover8,2,function(x){abs(x[5])}))-
             abs(apply(recover4,2,function(x){abs(x[5])})))/
        recover4[1,],1)
#            A            B  mean_v.true mean_v.false    sd_v.true           t0 
#          3.8         -1.5          1.0          0.6          0.3         -2.6

# However credible intervals are generally increased anywhere from 10% to 
# over 100%. 
round(apply(recover8,2,function(x){x[4]-x[2]})/
        apply(recover4,2,function(x){x[4]-x[2]}),2)
#            A            B  mean_v.true mean_v.false    sd_v.true           t0 
#         1.28         1.15         1.46         1.40         2.05         1.12 

### Example 1.9: Standard model with low-error; small-sample parameter recovery study ----

# Simulate 200 subjects, with ~ 2.5% errors on average
data9 <- h.simulate.dmc(modelLBA4,ps=pLBA8,ns=200,n=1e2)

mean(as.numeric(data9$S)==as.numeric(data9$R))
# [1] 0.976825

# The data-model object is now a list, with one subject in each slot.
dLBA9 <- data.model.dmc(data9,modelLBA4)

# Set up object ready for sampling (see dmc_4_4_LNRfixed for more details).
sLBA9 <- h.samples.dmc(nmc=120,p.priorLBA4,dLBA9,thin=10)


# Run the sampling, h.RUN.dmc just applies RUN.dmc to each subject. By
# default it uses one core per subject, which is more efficient than spreading
# chains across cores.
sLBA9 <- h.RUN.dmc(sLBA9,cores=60,meanN=2700,use.effectiveSize = FALSE)

# For all but t0 uncertainty is increased a little and bias is increased 
# substantially, but coverage remains good.
recover9 <- h.check.recovery.dmc(sLBA9,pLBA8)
#                     A      B mean_v.true mean_v.false sd_v.true     t0
# True            0.500  1.000       4.000        1.500     0.500  0.200
# 2.5% Estimate   0.245  0.545       2.866       -0.043     0.357  0.129
# 50% Estimate    0.589  1.086       4.721        1.803     0.624  0.214
# 97.5% Estimate  0.992  1.839       6.996        3.634     1.012  0.296
# Median-True     0.089  0.086       0.721        0.303     0.124  0.014
# Coverage       95.000 99.000      95.000       96.500    91.000 94.000

# Large increase in bias with the exception of t0
round(100*(abs(apply(recover9,2,function(x){abs(x[5])}))-
             abs(apply(recover5,2,function(x){abs(x[5])})))/
        recover5[1,],1)
#            A            B  mean_v.true mean_v.false    sd_v.true           t0 
#         17.7          7.5         16.1          9.1         19.9          1.3

# Credible intervals are generally increased by 20-80%.
round(apply(recover9,2,function(x){x[4]-x[2]})/
        apply(recover5,2,function(x){x[4]-x[2]}),1)
#            A            B  mean_v.true mean_v.false    sd_v.true           t0 
#          1.2          1.2          1.6          1.4          1.8          1.1

# Correlations are of comparable magnitude to 25% error case, although their 
# pattern shifts in some cases, and again the t0 distribution is clearly 
# truncated below by the prior.

pairs.dmc(sLBA9,thin=100,scale.subjects = FALSE)

### Example 1.10: Combined low/high error small-sample parameter recovery study ---- 

modelLBA10 <- model.dmc(
  p.map = list(A="1",B="1",mean_v=c("F","M"),sd_v="M",t0="1",st0="1"),
  constants=c(st0=0,sd_v.false=1),
  match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2"),F=c("f1","f2")),
  responses=c("r1","r2"),
  type="norm")


# Set sd_v.true to less than the fixed value sd_v.false=1 as is often seen
# in data (although usually a bit larger than the value of 0.5 used here)
pLBA10  <- c(A=.5,B=1,sd_v.true=0.5,mean_v.f1.true=4,mean_v.f2.true=4,
             mean_v.f1.false=3,mean_v.f2.false=1.5,t0=.2)

# Simulate 200 subjects, with ~ 2.5% errors on average in the low-error condition and
# ~ 25% errors on average in the high-error condition, 200 trials in total, 100
# per condition, and 50 in each cell of the 2*2 design.
data10 <- h.simulate.dmc(modelLBA10,ps=pLBA10,ns=200,n=1e2/2)

dLBA10 <- data.model.dmc(data10,modelLBA10)


p.priorLBA10 <- prior.p.dmc(
  p1=c(A=1,B=1,sd_v.true=1,mean_v.f1.true=2,mean_v.f2.true=2,
       mean_v.f1.false=1,mean_v.f2.false=1,t0=.3),                           
  p2=c(1,1,1,3,3,3,3,.25),lower=c(0,0,0,NA,NA,NA,NA,.1),upper=rep(NA,8)
)

# Set up objects ready for sampling (see dmc_4_4_LNRfixed for more details).
sLBA10 <- h.samples.dmc(nmc=120,p.priorLBA10,dLBA10,thin=10)


# Run the sampling, h.RUN.dmc just applies RUN.dmc to each subject. By
# default it uses one core per subject, which is more efficient than spreading
# chains across cores.
sLBA10 <- h.RUN.dmc(sLBA10,cores=60,meanN=2700,use.effectiveSize = FALSE)

# Coverage nominal.
recover10 <- h.check.recovery.dmc(sLBA10,pLBA10)
#                     A      B mean_v.f1.true mean_v.f2.true mean_v.f1.false mean_v.f2.false sd_v.true     t0
# True            0.500  1.000          4.000          4.000           3.000           1.500     0.500  0.200
# 2.5% Estimate   0.237  0.458          2.814          2.814           1.731          -0.715     0.361  0.143
# 50% Estimate    0.515  0.897          4.021          4.018           2.943           1.057     0.537  0.227
# 97.5% Estimate  0.810  1.493          5.447          5.441           4.295           2.545     0.766  0.300
# Median-True     0.015 -0.103          0.021          0.018          -0.057          -0.443     0.037  0.027
# Coverage       94.500 92.500         96.500         96.500          95.000          91.500    94.000 91.000

# The following compares results to those from 1.5 and 1.9 on shared parameters.

# Some increase in bias relative to high-error case in Example 1.5.
round(100*(abs(apply(recover10[,-6],2,function(x){abs(x[5])}))-
             abs(apply(recover5[,c(1:3,3,4:6)],2,function(x){abs(x[5])})))/
        recover5[1,c(1:3,3,4:6)],1)
#               A               B  mean_v.f1.true  mean_v.f2.true mean_v.f1.false       sd_v.true              t0 
#             2.9             9.2            -1.4            -1.5             0.9             2.4             7.6

# Larger decreases in bias relative to low-error case in Example 1.9 
round(100*(abs(apply(recover10[,-5],2,function(x){abs(x[5])}))-
             abs(apply(recover9[,c(1:3,3,4:6)],2,function(x){abs(x[5])})))/
        recover9[1,c(1:3,3,4:6)],1)
#               A               B  mean_v.f1.true  mean_v.f2.true mean_v.f2.false       sd_v.true              t0 
#           -14.7             1.7           -17.5           -17.6             9.4           -17.4             6.3 

# Uncertainty is generally similar to that in high-error case in Example 1.5 
round(apply(recover5[,c(1:3,3,4:6)],2,function(x){x[4]-x[2]})/
        apply(recover10[,-6],2,function(x){x[4]-x[2]}),1)
#            A            B  mean_v.true  mean_v.true mean_v.false    sd_v.true           t0 
#          1.1          1.1          1.0          1.0          1.0          0.9          1.0  

# However it is substantially better compared to low-error case in Example 1.9 
# (10-60% more uncertainty in Example 1.9)           
round(apply(recover9[,c(1:3,3,4:6)],2,function(x){x[4]-x[2]})/
        apply(recover10[,-5],2,function(x){x[4]-x[2]}),2)
#            A            B  mean_v.true  mean_v.true mean_v.false    sd_v.true           t0 
#         1.30         1.25         1.57         1.57         1.13         1.62         1.07

### Parameter testing

# By default compare.p tests the average over the 200 subjects of a contrast
# defined by two parameter names. It can also compute a contrast on a user 
# defined function. Here we compare mean_v as a function of factor F for the
# match and mismatch case. As expected, the there is no evidence for a difference
# between the matching parameters but strong evidence for the mismatching,
# although, as shown in the recovery table above (recover10), the f2 mismatching value is
# underestimated.

par(mfrow=c(1,2))
compare.p(sLBA10,pnames=c("mean_v.f1.true","mean_v.f2.true"),
          pretty=c("v.f1","v.f2"),show.plot=TRUE)
#        v.f1  v.f2 contrast
# 2.5%  3.955 3.951   -0.012
# 50%   4.050 4.046    0.003
# 97.5% 4.138 4.133    0.020
# p.gt.0 
#  0.674
compare.p(sLBA10,pnames=c("mean_v.f1.false","mean_v.f2.false"),
          pretty=c("v.f1","v.f2"),show.plot=TRUE)
#        v.f1 v.f2 contrast
# 2.5%  2.870 0.90    1.851
# 50%   2.961 1.02    1.940
# 97.5% 3.046 1.14    2.037
# p.gt.0 
#      1 

# Here we get the p values for each individual subject. This requires giving
# compare.p a list of length one containing each of the subjects values.
p.v.true.f1.gt.f2 <- unlist(lapply(sLBA10,function(x){
  attr(compare.p(list(x),pnames=c("mean_v.f1.true","mean_v.f2.true"),verbose=FALSE),"p.gt.0")}))

# Consistent with a well calibrated test, the distribution of p values is 
# uniform for this null comparison
hist(p.v.true.f1.gt.f2,breaks=seq(0,1,by=.05),main="mean_v.true",
     xlab="p(f1>f2)")


# For the clear difference between false rates, even the least favourable p
# value strongly supports f1 > f2.
p.v.false.f1.gt.f2 <- unlist(lapply(sLBA10,function(x){
  attr(compare.p(list(x),pnames=c("mean_v.f1.false","mean_v.f2.false"),verbose=FALSE),"p.gt.0")}))
min(p.v.false.f1.gt.f2)
# [1] 0.9982639

### Example 1.11: Fits to low and high error data with simpler (true) model ----

# In this model, levels of mean_v are set as a constant=Inf to be replaced with
# a separately estimated parameter mean_v1. In this case, we use it to replace 
# mean_v.f1.true and mean_v.f2.true with the same value, corresponding to the
# data generating model where these two values are identical.
load_model ("LBA","lba_B2v.R")

modelLBA11 <- model.dmc(
  p.map = list(A="1",B="1",mean_v=c("F","M"),mean_v1="1",sd_v="M",t0="1",st0="1"),
  constants=c(st0=0,sd_v.false=1,mean_v.f1.true=Inf,mean_v.f2.true=Inf),
  match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2"),F=c("f1","f2")),
  responses=c("r1","r2"),
  type="norm")
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"               "B"               "mean_v.f1.false" "mean_v.f2.false"
# [5] "mean_v1"         "sd_v.true"       "t0"             
# 
# Constants are (see attr(,"constants") ):
#            st0     sd_v.false mean_v.f1.true mean_v.f2.true 
#              0              1            Inf            Inf 
# 
# Model type = norm (posdrift= TRUE ) 


dLBA11 <- data.model.dmc(data10,modelLBA11)

p.priorLBA11 <- prior.p.dmc(
  p1=c(A=1,B=1,sd_v.true=1,mean_v1=2,mean_v.f1.false=1,mean_v.f2.false=1,t0=.3),                           
  p2=c(1,1,1,3,3,3,.25),lower=c(0,0,0,NA,NA,NA,.1),upper=rep(NA,7)
)
par(mfcol=c(2,4)); for (i in names(p.priorLBA11)) plot.prior(i,p.priorLBA11)

# Run the sampling.
sLBA11 <- h.samples.dmc(nmc=120,p.priorLBA11,dLBA11,thin=10)
sLBA11 <- h.RUN.dmc(sLBA10,cores=60,meanN=2700,use.effectiveSize = FALSE)


# Coverage is again nominal
pLBA11  <- c(A=.5,B=1,mean_v.f1.false=3,mean_v.f2.false=1.5,mean_v1=4,
             sd_v.true=0.5,t0=.2)
recover11 <- h.check.recovery.dmc(sLBA11,pLBA11)
#                     A      B mean_v.f1.false mean_v.f2.false mean_v1 sd_v.true     t0
# True            0.500  1.000           3.000           1.500   4.000     0.500  0.200
# 2.5% Estimate   0.227  0.491           1.766          -0.640   2.847     0.361  0.134
# 50% Estimate    0.519  0.962           3.046           1.165   4.106     0.538  0.218
# 97.5% Estimate  0.834  1.597           4.495           2.719   5.625     0.767  0.294
# Median-True     0.019 -0.038           0.046          -0.335   0.106     0.038  0.018
# Coverage       94.000 95.500          95.500          90.500  95.500    93.500 93.500


# Generally less bias in simpler model
round(100*(abs(apply(recover11,2,function(x){abs(x[5])}))-
             abs(apply(recover10[,c(1:2,5:6,3,7:8)],2,function(x){abs(x[5])})))/
        recover11[1,],1)
#               A               B mean_v.f1.false mean_v.f2.false         mean_v1       sd_v.true              t0 
#             0.7            -6.6            -0.4            -7.2             2.1             0.2            -4.6 

# Uncertainty is similar 
round(apply(recover11,2,function(x){x[4]-x[2]})/
        apply(recover10[,c(1:2,5:6,3,7:8)],2,function(x){x[4]-x[2]}),2)
#               A               B mean_v.f1.false mean_v.f2.false         mean_v1       sd_v.true              t0 
#            1.06            1.07            1.06            1.03            1.05            1.01            1.02 

### MODEL RECOVERY STUDY

# Both the models in 1.10 and 1.11 are "true", but 1.11 is simpler. Model 
# selection criteria such as DIC and WAIC should prefer the simpler model.
# Lets check how well they do. Note this is a difficult case as 1.10 must fit
# at least as well, and probably a bit better as it can soak up some noise.

# These functions report the summed DIC for all 200 subject and save off each
# individual DIC. In terms of the sum DIC clearly picks the simpler model by 100.
DIC11 <- h.IC.dmc(sLBA11,DIC=TRUE)
# Summed Minimum Deviance and DIC
# [1] -98710.39 -96943.90
DIC10 <- h.IC.dmc(sLBA10,DIC=TRUE)
# Summed Minimum Deviance and DIC
# [1] -98938.54 -96842.38

# DIC based model selection was successful 75% of the time.
mean((DIC11[,2]-DIC10[row.names(DIC11),2])<0)
# [1] 0.77


BPIC11 <- h.IC.dmc(sLBA11)
# Summed Minimum Deviance and BPIC
# [1] -98710.39 -95935.04
BPIC10 <- h.IC.dmc(sLBA10)
# Summed Minimum Deviance and BPIC
# [1] -98938.54 -95697.53

# BPIC based model selection was successful even more.
mean((BPIC11[,2]-BPIC10[row.names(DIC11),2])<0)
# [1] 0.83


# WAIC is harder to get, we must first compute "pointwise likelihoods", which 
# produces a large object, so we thin to do this only for some posterior 
# samples (here we choose the thinning to make the total size equal and use
# many cores to speed it up).

# Divide 21 chains x 160/chain by 70 = 48 samples per person x 200 data points = 
# 9600 point-wise likelihoods per participant x 200 participants = 1.92m values.
pll.11 <- group_trial_log_likes(sLBA11,thin_pointwise = 70,cores=50)
# Divide 24 chains x 120/chain by 60 = 48 samples per person x 200 data points = 
# 9600 point-wise likelihoods per participant x 200 participants = 1.92m values
pll.10 <- group_trial_log_likes(sLBA10,thin_pointwise = 60,cores=50)

# This function gets the sum over subjects and picks the simpler model by 100 
waic.dmc(pll.11,digits=2,save=TRUE)
#       p    se_p    waic se_waic 
#     972      15  -96890     369
waic.dmc(pll.10,digits=2,save=TRUE)
#       p    se_p    waic se_waic 
#    1115      15  -96795     369

# We have to get the indivdiual WAICs by repeating the whole operation using
# the following function that first gets the pointwise likelihoods for each
# subejct.
get.waic <- function(samples,thin_pointwise=1)
  waic.dmc(trial_log_likes(samples,thin_pointwise = thin_pointwise),save=TRUE)$waic   

# This is slow as it is not parallelized. 
WAICs10 <- unlist(lapply(sLBA10,get.waic,thin_pointwise=24))
WAICs11 <- unlist(lapply(sLBA11,get.waic,thin_pointwise=28))

# WAIC is poor at picking the simpler model (on some other runs it was as
# as low as .51).
mean((WAICs11-WAICs10)<0)
# [1] 0.67

save_data(sLBA0,sLBA1,sLBA2,sLBA3a,sLBA3b,sLBA3c,sLBA3d,
          sLBA4,sLBA5,recover5,sLBA6,sLBA7,sLBA8,recover8,sLBA9,recover9,
          sLBA10,recover10,sLBA11,recover11,WAICs10,WAICs11,
          file="DMCpaper1.RData")

