##################  DMC Lesson 3: Sampling

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# THESE LESSONS REQUIRE PACKAGES # 
#      snowfall AND rlecuyer     #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

### Lesson 3.1:  Sampling a single LNR subject 

# In a Bayesian context we fit a model and obtain estimates of its parameters
# that best describe the data by sampling from the model. This lesson explains
# how to sample. It assumes you know about things like "burnin" (the need to get
# some early samples that are later thrown away) "chains" (series of sampled 
# parameters) and the DE-MCMC method of obtaining samples, see Turner, B. M., 
# Sederberg, P. B., Brown, S. D., & Steyvers, M. (2013). A method for 
# efficiently sampling from distributions with correlated dimensions. 
# Psychological Methods, 18(3), 368–384. http://doi.org/10.1037/a0032222 

rm(list=ls())
# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("LNR","lnr.R")

# load_data("dmc_3_1.RData")

# A simple 5 parameter LNR model with normal/gamma/uniform priors for 
# meanlog/sdlog/t0 respectively
model <- model.dmc(p.map=list(meanlog="M",sdlog="M",t0="1",st0="1"),
                   match.map=list(M=list(s1=1,s2=2)),constants=c(st0=0),
                   factors=list(S=c("s1","s2")),responses=c("r1","r2"),
                   type="lnr")
p.vector  <- c(meanlog.true=-1,meanlog.false=0,
               sdlog.true=1,sdlog.false=1,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
p.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","gamma","gamma","beta"),
  p1=c(meanlog.true=-1,meanlog.false=0, sdlog.true=2,sdlog.false=2,t0=1),                           
  p2=c(1,1,1,1,1),lower=c(NA,NA,NA,NA,.1),upper=c(NA,NA,NA,NA,1)
)
par(mfcol=c(2,3)); for (i in names(p.prior)) plot.prior(i,p.prior)


### SETUP SAMPLING

# Sets up an object to contain the results of taking n.chains of length nmc 
# samples based on p.prior and data.model (the data+model object). By default 
# uses 3*n.pars chains (set with n.chains argument). A rule of thumb is that the
# number of chains should be at least 2*n.pars + 1 (apart from computation,
# more is usually better).

samples <- samples.dmc(nmc=100,p.prior,data.model)

# Samples is a list with components: 
# 1) theta: array (n.chains x n.pars x nmc) of samples
# 2) summed_log_prior: nmc x n.chains matrix of log-prior summed over data points
# 3) log_likelihoods: nmc x n.chains matrix of log-likelihoods summed over data points
# 4) data and p.prior objects
# 5) rp = random perturbation used by DEMCMC, default = .001
# 6) p.names: names of parameters
# 7) constants: n.chains, n.pars, nmc, start (index to start sampling after)

# Sometimes you might want to get start points more tightly clustered 
# around reasonable values than will be obtained from the prior. To do 
# this you can use the argument "start.prior" This has the same form as
# p.prior but you can specify tighter distributions.

start.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","gamma","gamma","beta"),
  p1=c(meanlog.true=-1,meanlog.false=0, sdlog.true=2,sdlog.false=2,t0=1),                           
  p2=c(1,1,1,1,1)/5,lower=c(NA,NA,NA,NA,.1),upper=c(NA,NA,NA,NA,1)
)
samples1 <- samples.dmc(nmc=100,p.prior,data.model,start.prior=start.prior)

# NB1: By default samples.dmc samples from the prior until it gets parameters 
#      that return finite posterior likelihood values (and returns an error if
#      it fails after 10,000 attempts). If the argument theta1 (matrix nchains x
#      n.pars) is non-null these parameters are used to try and initialise the 
#      chains instead. Parameters names must be set for thhe column of theta1.
theta1 <- matrix(rep(p.vector,each=15),nrow=15,
                 dimnames=dimnames(samples$theta)[1:2])
samples_theta1 <- samples.dmc(nmc=100,p.prior,data.model,n.chains=15,theta1=theta1)

# RUN SAMPLING

# run.dmc runs nmc sampling steps, reporting each report (default = 10) update.

# The DEMCMC algorithm it uses usually needs no tuning. The only tuning 
# argument is "gamma.mult" with default value 2.38. The DEMCMC gamma 
# parameter is set to gamma.mult/sqrt(2*np) where np is the number of 
# parameters being sampled. Smaller values of gamma.mult cause smaller changes
# in parameters and may be useful if chains are unchanging due to trying large
# jumps that usually get rejected. gamma.mult=NA causes gamma to be sampled
# uniformly between 0.5 - 1, which can sometimes be useful.
#
# If argument farjump is not na (default) on every farjump iteration 
# gamma.mult=.98. Ter Baark (2006) recommended this for multi-modal posteriors
# with farjump=10 (i.e., every 10th iteration).


# Turn on migration step with probability p.migrate = 0.05 to assist burnin, 
# helps to pull in "stuck" chains that wander far from the other chains.
#
# NB: Setting p.migrate to high can cause all chains to quickly converge to a 
#     fairly good but sub-optimal solution that then takes quite a bit more 
#     sampling before the chains move to a better solution (and during that
#     time you can be fooled that convergence has been achieved.

# You can record how long this takes using "system.time"

system.time({
  samples <- run.dmc(samples, report = 10,p.migrate=.05)
})

# "run.dmc" can also do chain updates in parallel using the "snowfall" 
# package. Here is run using 4 cores.

system.time({
  samples1 <- run.dmc(samples, report = 10, cores=4,p.migrate=.05)
})

# The reduction in time is never by a factor = cores. The inefficiency
# is because of the cost of starting up the parallel cluster (which is less
# important for longer runs, but here it dominates) and of scattering and 
# gathering the chain updates.

# Use the "coda" package to view chains.
plot.dmc(samples)

# By default CODA usually draws smooth lines through chains and also plots 
# corresponding densities. These are turned off by dmc but can be turned
# on as follows. Note that with the extra density plots you will get
# two pages of output.
plot.dmc(samples,smooth=TRUE,density=TRUE)

# Any arguments to the coda plot function can be passed in this way. For example
# the ask argument can be used to pause after each page is plotted.
plot.dmc(samples,smooth=TRUE,density=TRUE,ask=TRUE)

# Alternately you can use the layout command to get it all on one page 
plot.dmc(samples,smooth=TRUE,density=TRUE,layout=c(2,5))

# NB: You can turn the parameter samples (i.e., samples$theta) into a coda 
#     mcmc.list object using the theta.as.mcmc.list(samples) function, and
#     directly apply coda's extensive range of diagnostic functions.  
# Here you see that plot now recognizes the coda object and plots with defaults.
plot(theta.as.mcmc.list(samples))

# Returning to our example, the effect of the narrower priors can be seen in 
# less initial variability and quicker convergence in samples1, although the
# effect is not very marked.
plot.dmc(samples1)

# If the samples argument is provided with add=TRUE then run.dmc adds nmc 
# samples, using the last values obtained in the previous sampling as the 
# starting point
#
samples2 <- run.dmc(samples.dmc(nmc=100,samples=samples1,add=TRUE), 
                    p.migrate=.05)

# Use the start argument to look at only the extra samples (you can also specify 
# an upper bound to plot using "end"). 
plot.dmc(samples2,start=101)

# It is often useful to plot the posterior log-likelihood for each chain. We 
# expect this should increase steadily become stable.
plot.dmc(samples2,pll.chain=TRUE)

# If we focus on the final 100 we can see migration pulling in a wayward chain.
plot.dmc(samples2,pll.chain=TRUE,start=101)

# If we focus on the last few iterations everything looks good.
plot.dmc(samples2,start=175)

# Given no chains appear stuck we can now turn of migration. Migration speeds up
# convergence to the posterior mode and helps pull in stuck chains but can 
# result in samples not being a proper representation of the posterior (usually
# by reducing their variability). Hence we now add 500 more samples, but discard 
# the first 200 (add = FALSE is the default, and replaces the existing samples). 
samples3 <- run.dmc(samples.dmc(nmc=500,samples=samples2), 
                    cores=4,report=50)

# The visual appearance of the chains ("fat hairy caterpillars") indicates
# convergence.
plot.dmc(samples3)

# The next lesson indicates how to test convergence more formally and how to 
# use the information in the chains.

# Before we finish note that samples.dmc allows you selectively remove samples
# using the "remove" argument, which specifies a range of samples that will be
# omitted. 

# To preserve the remaining samples and not include room for any new samples 
# you use add=TRUE with nmc=0 (of course if you want to add room for new samples
# then set nmc appropriately). 

# For example, to remove the first two hundred samples:
samples3.short <- samples.dmc(nmc=0,samples=samples3,add=TRUE,remove=1:200)
plot.dmc(samples3.short)

# You can remove from anywhere, for example:
samples3.short <- samples.dmc(nmc=0,samples=samples3,add=TRUE,remove=c(201:300,450:490))
plot.dmc(samples3.short)

# save_data (p.prior,data.model,samples,samples1,samples2,samples3,file="dmc_3_1.RData")

### ADVANCED CONTENT, NOT NEEDED THE FIRST TIME AROUND.

# On occasion you may have knowledge of the likely range of values of a parameter
# (e.g., from prior literature) and want to assume that rather than sampling the
# parameter. In this case you need to incorporate a "constant.prior" in the model.

# Suppose you wanto to fix meanlog.true ~ N(-1,.01) and meanlog.false ~ N(0,.01)
constant.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm"),
  p1=c(meanlog.true=-1,meanlog.false=0),                           
  p2=c(.1,.1),lower=c(NA,NA),upper=c(NA,NA)
)

par(mfcol=c(1,2)); for (i in names(constant.prior)) plot.prior(i,constant.prior)

# On each iteration of sampling the values of contant prior paramaters are
# replaced with samples form the constant prior, e.g., 
rprior.dmc(constant.prior)
#      meanlog.true meanlog.false
# [1,]    -1.009186  -0.005842251  

# Now incorporate it into the model
modelCP <- model.dmc(p.map=list(meanlog="M",sdlog="M",t0="1",st0="1"),
                     match.map=list(M=list(s1=1,s2=2)),constants=c(st0=0),
                     constant.prior=constant.prior,
                     factors=list(S=c("s1","s2")),responses=c("r1","r2"),
                     type="lnr")
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "sdlog.true"  "sdlog.false" "t0"         
# 
# Constants are (see attr(,"constants") ):
# st0 
#   0 
# 
# Constant prior for parameters: meanlog.true meanlog.false 
# 
# Model type = lnr 

# The required prior for sampling omits the constant prior parameters
p.priorCP <- prior.p.dmc(
  dists = c("gamma","gamma","beta"),
  p1=c(sdlog.true=2,sdlog.false=2,t0=1),                           
  p2=c(1,1,1),lower=c(NA,NA,.1),upper=c(NA,NA,1)
)


# Make up a data model object *using data from full model*
data.modelCP <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),modelCP)

samplesCP <- samples.dmc(nmc=200,p.priorCP,data.modelCP)
samplesCP <- run.dmc(samplesCP,p.migrate=.05,cores=4)
plot.dmc(samplesCP)
# Seems to have increased autocorrelation quite a lot, so thin
samplesCP <- run.dmc(samples.dmc(samples=samplesCP,nmc=200,thin=10),cores=4)
plot.dmc(samplesCP)

# Good recovery
check.recovery.dmc(samplesCP,p.vector[-c(1:2)])
#                sdlog.true sdlog.false   t0
# True                 1.00        1.00 0.20
# 2.5% Estimate        0.81        0.89 0.18
# 50% Estimate         0.98        1.02 0.20
# 97.5% Estimate       1.20        1.23 0.21
# Median-True         -0.02        0.02 0.00

# Still very high autocorrelation
acf.dmc(samplesCP)

# Similar but slightly worse performance (as expected) if data is from true model,
# as follows. Checking is left as an exercise.
data.modelCP <- data.model.dmc(simulate.dmc(p.vector[-c(1:2)],modelCP,n=1e4),modelCP)


