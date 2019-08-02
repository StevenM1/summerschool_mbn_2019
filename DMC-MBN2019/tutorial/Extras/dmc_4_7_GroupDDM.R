##################  DMC Lesson 4: Multiple Subjects

### Lesson 4.7groupDDM: Hierarchical DDM with sampling failure & group parameters

# Shows failure of hierarchical DDM for a simple one factor 2 level with the 
# same number of samples overall as 4_7. Group parameter for both sv and sz
# fix these problems. Note that this example uses a small sv population mean,
# similar problems don’t occur with larger sv. Chain update mixing in DMC 
# (defualt random.phi=TRUE and random.theta=TRUE) solves the convergence 
# problem, but parameter recovery still difficult.


rm(list=ls())

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("DDM","ddm.R") 

# load_data ("dmc_4_7group_fixed.RData")
# load_data ("dmc_4_7group_random.RData")

## Set up a DDM Model, exactly as Lesson 3.4
model <- model.dmc(
  p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
  match.map = list(M=list(s1="r1",s2="r2")),
  factors   = list(S=c("s1","s2")),
  constants = c(st0=0,d=0),
  responses = c("r1","r2"),
  type = "rd")

# Population distribution

pop.mean <- c(a=1,  v=1, z=0.5,sv=0.25, sz=0.3, t0=0.3)
pop.scale <-c(a=0.2,v=.2,z=0.1,sv=.05,sz=0.05,t0=0.05)
pop.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm","tnorm","tnorm"),
  p1=pop.mean,                           
  p2=pop.scale,
  lower=c(0,-5, 0, 0, 0, 0),
  upper=c(2, 5, 1, 2, 1, 1)
)
##  Check population distributions
# par(mfcol=c(2,3)); for (i in names(pop.prior)) plot.prior(i,pop.prior)


# Simulate some data
raw.data <- h.simulate.dmc(model, p.prior = pop.prior, n = 500, ns = 40)
data.model <- data.model.dmc(raw.data, model)

# # Take a look at the first 10 subjects data
# par(mfrow=c(2,5)) # Row 1 = subjects 1..5, row2 = 6..10
# for (i in 1:10)
#   plot.cell.density(data.cell=data.model[[i]][data.model[[i]]$S=="s1",],C="r1")


# Take a look at parameters
ps <- round( attr(raw.data, "parameters"), 2); 
round(apply(ps,2,mean),2)
#    a    v    z   sv   sz   t0 
# 1.03 1.00 0.50 0.25 0.30 0.31 
round(apply(ps,2,sd),2)
#    a    v    z   sv   sz   t0 
# 0.19 0.19 0.09 0.04 0.05 0.05 


### FIT FIXED EFFECTS

# specify a broader prior than the true population distribution
p.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm","tnorm","tnorm"),
  p1=pop.mean,                           
  p2=pop.scale*5,
  lower=c(0,-5, 0, 0, 0, 0),
  upper=c(2, 5, 1, 2, 1, 1)
)

# Start with thin of 5 assuming at least as bad as LNR  
samples  <- h.samples.dmc(nmc = 50, p.prior, data.model, thin = 5)

save(raw.data,samples,file="4_7group.RData")

# Remember there is no output as this runs, except after each subject is done
# the number of iteraitons for that subject is printed; it will take a while!
samples  <- h.run.dmc(samples, p.migrate = .05, cores = 4)

# No stuck chains by 50
unlist(lapply(samples,pick.stuck.dmc,start=50))

# Running the following for each chain suggests there might be a little 
# autocorrelation left at lag=2, increase thin to 10
for (i in 1:40) acf.dmc(samples[[i]],start=50,par=NA,chain=1)

# Now run without migration, and turning up thinning to 10
samples1 <- h.run.dmc(h.samples.dmc(nmc=500, samples=samples, thin=10), cores=4)

# Nothing has got stuck
unlist(lapply(samples1,pick.stuck.dmc,start=1))

# Paramter chains look well converged and mixed
plot.dmc(samples1,density=FALSE,smooth=FALSE,subject=1,layout=c(2,3))

# All chains are well converged. 
gelman.diag.dmc(samples1)
#   38   12    1   35   15   23    8   26   37   32   31   20   18   19    9 
# 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.01 1.01 1.01 1.01 1.01 1.01 
#    2   16   40    3   36    5   34   13   24   28   33   29   17   25   22 
# 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 
#   21   10    6    4    7   11   39   27   14   30 
# 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 
# Mean
# [1] 1.01


# Effective sample size: good yield given 18 x 500 = 8000 samples
es <- effectiveSize.dmc(samples1)
apply(do.call(rbind,es),2,min)
#    a    v    z   sz   sv   t0 
# 2905 2979 2846 1975 2844 2289 

# Good fits?
pp <- h.post.predict.dmc(samples1) 
tmp <- lapply(pp, function(x){plot.pp.dmc(x, style="cdf") })

# Check against true values: fairly good parameter recovery.
est <- t(data.frame(lapply(summary.dmc(samples1),function(x){
  x$statistics[,"Mean"]})))[,dimnames(ps)[[2]]]

# Pretty good, excep sv and sz
round(apply(est,2,mean),3)
#     a     v     z    sv    sz    t0 
# 1.027 0.985 0.499 0.311 0.252 0.312 
round(apply(est,2,mean)-apply(ps,2,mean),3)
#      a      v      z     sv     sz     t0 
#  0.001 -0.019  0.002  0.065 -0.044 -0.001 

# Scale for sz and sv off by a bit
round(apply(est,2,sd),3)
#     a     v     z    sv    sz    t0 
# 0.189 0.193 0.089 0.063 0.062 0.054 
round(apply(est,2,sd)-apply(ps,2,sd),3)
#      a      v      z     sv     sz     t0 
#  0.001 -0.001  0.000  0.020  0.010  0.000 

save_data(raw.data,samples,samples1,pp,est,file="dmc_4_7group_fixed.RData")

### FIT RANDOM EFFECTS

# Use all truncated normal priors for locations
mu.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm","tnorm","tnorm"),
  p1=pop.mean,                           
  p2=pop.scale*5,
  lower=c(0,-5, 0, 0, 0, 0),
  upper=c(2, 5, 1, 2, 1, 1)
)
# par(mfcol=c(2,3)); for (i in names(mu.prior)) plot.prior(i,mu.prior)

# As in the LBA example suppose we have little idea about individual variation 
# apart from being sure that standard deviations are not huge. Again, we use
# uniform priors for the hyper-scale parameters on the range 0-1. Note that
# similar results to those reported below were obtained with a range of
# different vague priors.
sigma.prior <- prior.p.dmc(
  dists = rep("beta", length(p.prior)),
  p1=c(a=1, v=1,z=1,sv=1, sz=1, t0=1),p2=c(1,1,1,1,1,1)
)
# par(mfcol=c(2,3)); for (i in names(sigma.prior)) plot.prior(i,sigma.prior)

# Make a hyper-prior list
pp.prior <- list(mu.prior, sigma.prior)


# Make a new samples object from the 40 subject data, longer and with more
# immediate thinning that LNR given likely more autocorrelation
hsamples <- h.samples.dmc(nmc = 50, p.prior = p.prior, pp.prior=pp.prior, 
                          data = data.model, thin = 10)

save_data(hsamples,file="4_7group_test.RData")
load_data("4_7group_test.RData")

### Fit with default setting of random.phi=TRUE, and random.theta=TRUE, which 
### mixes up links between hyper and lower level chains in update of hyper and 
### data data parameters respectively. This is the default setting for DMC as it 
### fixes convergence problems identified in the next section.

hsamples.good <- hsamples

# Fit with migration at both levels
hsamples.good <- h.run.dmc(hsamples, cores=8, report=1, 
                           p.migrate = .05, h.p.migrate = .05)
hsamples.good <- h.run.dmc(h.samples.dmc(nmc=150,samples=hsample.good,add=TRUE),
                           p.migrate = .05, h.p.migrate = .05,cores=4,report=10)

# Hyper likelihoods improve gradually
plot.dmc(hsamples.good,hyper=TRUE,smooth=FALSE,density=FALSE,pll.chain=TRUE,start=11)

# Try another, 500 but as no stuck chains, so turn off migration.
hsamples1.good <- h.run.dmc(h.samples.dmc(nmc = 500, samples = hsamples.good, thin = 10), 
                            cores = 8, report = 1)
# Big drop in likeihood, so add another 500
hsamples1.good <- h.run.dmc(h.samples.dmc(nmc = 500, samples = hsamples1.good, add = TRUE), 
                            cores = 8, report = 1)

# After drop evens out form 300 
plot.dmc(hsamples1.good,hyper=TRUE,smooth=FALSE,density=FALSE,pll.chain=TRUE)
# Mixing nicely but sv very spread out
plot.dmc(hsamples1.good,hyper=TRUE,layout=c(2,6),smooth=FALSE,density=FALSE,start=300)

# Everything looks good 
gelman.diag.dmc(hsamples1.good,hyper=TRUE)
#       Point est. Upper C.I.
# a.h1        1.00       1.00
# v.h1        1.00       1.00
# z.h1        1.00       1.00
# sv.h1       1.00       1.00
# sz.h1       1.01       1.01
# t0.h1       1.00       1.00
# a.h2        1.00       1.00
# v.h2        1.00       1.00
# z.h2        1.00       1.00
# sv.h2       1.00       1.00
# sz.h2       1.00       1.00
# t0.h2       1.00       1.00
# 
# Multivariate psrf
# 
# 1

hest.good <- summary.dmc(hsamples1.good,hyper=TRUE)
# Mean recovery not to bad except for sv greatly under-estimated.
round(hest.good$statistics[,"Mean"][1:6],3)
#  a.h1  v.h1  z.h1 sv.h1 sz.h1 t0.h1 
# 1.023 0.973 0.499 0.121 0.252 0.312 
round(hest.good$statistics[,"Mean"][1:6]-apply(ps,2,mean),3)
#   a.h1   v.h1   z.h1  sv.h1  sz.h1  t0.h1 
# -0.003 -0.032  0.002 -0.125 -0.044 -0.002 

# OK except sv greatly over-estimated, sz a bit over.
round(hest.good$statistics[,"Mean"][7:12],3)
#  a.h2  v.h2  z.h2 sv.h2 sz.h2 t0.h2 
# 0.192 0.185 0.092 0.234 0.076 0.056 
round(hest.good$statistics[,"Mean"][7:12]-apply(ps,2,sd),3)
#   a.h2   v.h2   z.h2  sv.h2  sz.h2  t0.h2 
#  0.004 -0.009  0.003  0.191  0.023  0.002 


### Run with random.phi=FALSE and random.theta=FALSE, keeps chains at both 
### levels in a 1:1 relationship

# Fit with migration at both levels
hsamples <- h.run.dmc(hsamples, cores=8, report=1,
                      random.phi=FALSE,random.theta=FALSE,p.migrate = .05, h.p.migrate = .05)
hsamples <- h.run.dmc(h.samples.dmc(nmc=150,samples=hsample,add=TRUE),
                      p.migrate = .05, h.p.migrate = .05,cores=4,report=10,
                      random.phi=FALSE,random.theta=FALSE)

# Hyper likelihoods improve gradually
plot.dmc(hsamples,hyper=TRUE,smooth=FALSE,density=FALSE,pll.chain=TRUE,start=11)

# Try another, 500 but as no stuck chains, so turn off migration.
hsamples1 <- h.run.dmc(h.samples.dmc(nmc = 500, samples = hsamples, thin = 10), 
                       cores = 8, report = 1, random.phi=FALSE,random.theta=FALSE)

# Mixing fails dramatically. 
plot.dmc(hsamples1,hyper=TRUE,smooth=FALSE,density=FALSE,pll.chain=TRUE)
# For sv and sz some stay in zero variance trap, some gradually move out.
plot.dmc(hsamples1,hyper=TRUE,layout=c(2,6),smooth=FALSE,density=FALSE)

# Everything looks good except sv and sz
gelman.diag.dmc(hsamples1,hyper=TRUE)
#       Point est. Upper C.I.
# a.h1        1.00       1.00
# v.h1        1.01       1.01
# z.h1        1.00       1.00
# sv.h1       1.28       1.48
# sz.h1       1.20       1.45
# t0.h1       1.00       1.00
# a.h2        1.00       1.00
# v.h2        1.01       1.01
# z.h2        1.00       1.00
# sv.h2       1.70       2.21
# sz.h2       1.89       2.61
# t0.h2       1.00       1.00
# 
# Multivariate psrf
# 
# 2

hest <- summary.dmc(hsamples1,hyper=TRUE)
# Mean recovery not to bad except for sv, much worse than random.phi/theta=TRUE
round(hest$statistics[,"Mean"][1:6],3)
#  a.h1  v.h1  z.h1 sv.h1 sz.h1 t0.h1 
# 1.019 0.961 0.499 0.046 0.265 0.312 
round(hest$statistics[,"Mean"][1:6]-apply(ps,2,mean),3)
#   a.h1   v.h1   z.h1  sv.h1  sz.h1  t0.h1 
# -0.007 -0.043  0.002 -0.200 -0.030 -0.001

# OK for all but sz, recovery of sv actually better than for random.phi/theta=TRUE
round(hest$statistics[,"Mean"][7:12],3)
#  a.h2  v.h2  z.h2 sv.h2 sz.h2 t0.h2 
# 0.191 0.186 0.092 0.041 0.012 0.056 
round(hest$statistics[,"Mean"][7:12]-apply(ps,2,sd),3)
#   a.h2   v.h2   z.h2  sv.h2  sz.h2  t0.h2 
#  0.003 -0.008  0.003 -0.002 -0.041  0.002 


### FIT RANDOM EFFECTS WITH sv and sz AS GROUP PARAMETERS
### NB: Run with random.phi/theta=FALSE to observe if group parameters fix mixing.

# Set sv and sz to a constant in model (value doesn’t matter)
model.group <- model.dmc(
  p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
  match.map = list(M=list(s1="r1",s2="r2")),
  factors   = list(S=c("s1","s2")),
  constants = c(st0=0,d=0,sz=.1,sv=.1),
  responses = c("r1","r2"),
  type = "rd")


# Data level prior drops sz and sv
p.prior.group <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm"),
  p1=pop.mean[-c(4:5)],                           
  p2=pop.scale[-c(4:5)]*5,
  lower=c(0,-5, 0, 0),
  upper=c(5, 7, 1, 1)
)

# Drop sv and sz in sigma prior
sigma.prior.group <- prior.p.dmc(
  dists = rep("beta", length(p.prior.group)),
  p1=c(a=1, v=1,z=1,t0=1),p2=c(1,1,1,1),
  upper=c(2,2,1,1)
)
# par(mfcol=c(2,2))
# for (i in names(sigma.prior.group)) plot.prior(i,sigma.prior.group)

data.model.group <- data.model.dmc(raw.data, model.group)
pp.prior.group <- list(mu.prior, sigma.prior.group)


hsamples.group <- h.samples.dmc(nmc = 5, p.prior=p.prior.group, 
                                data=data.model.group, pp.prior=pp.prior.group,thin=10,
                                random.phi=TRUE,random.theta=FALSE)
hsamples.group <- h.run.dmc(hsamples.group, cores=8, report=1, 
                            p.migrate = .05, h.p.migrate = .05,random.phi=TRUE,random.theta=FALSE)


hsamples.group <- h.run.dmc(h.samples.dmc(nmc=195,samples=hsamples.group,add=TRUE),
                            cores=8,report=1,p.migrate = .05, h.p.migrate = .05,
                            random.phi=TRUE,random.theta=FALSE)

# Stable by 100
plot.dmc(hsamples.group,hyper=TRUE,smooth=FALSE,density=FALSE,pll.chain=TRUE,start=10)
# Chains look nice
plot.dmc(hsamples.group,hyper=TRUE,layout=c(2,5),smooth=FALSE,density=FALSE,start=100)

# Try another, 500 but as no stuck chains, so turn off migration.
hsamples1.group <- h.run.dmc(h.samples.dmc(nmc = 500, samples = hsamples.group, thin = 10), 
                             cores = 8, report = 1,random.phi=TRUE,random.theta=FALSE)

# Mixing good. 
plot.dmc(hsamples1.group,hyper=TRUE,smooth=FALSE,density=FALSE,pll.chain=TRUE)
# Chains look nice.
plot.dmc(hsamples1.group,hyper=TRUE,layout=c(2,5),smooth=FALSE,density=FALSE)

# Everything looks good
gelman.diag.dmc(hsamples1.group,hyper=TRUE)
#       Point est. Upper C.I.
# a.h1        1.00       1.00
# v.h1        1.00       1.00
# z.h1        1.00       1.00
# sv.h1       1.01       1.01
# sz.h1       1.02       1.04
# t0.h1       1.00       1.00
# a.h2        1.00       1.00
# v.h2        1.00       1.00
# z.h2        1.00       1.00
# t0.h2       1.00       1.00
# 
# Multivariate psrf
# 
# 1.02

gest <- summary.dmc(hsamples1.group,hyper=TRUE)
# Mean recovery not too bad except for sv (a bit worse than individual) and sz
# (a bit better than indivdiual)
round(gest$statistics[,"Mean"][1:6],3)
#  a.h1  v.h1  z.h1 sv.h1 sz.h1 t0.h1 
# 1.025 0.983 0.500 0.349 0.299 0.313 
round(hest$statistics[,"Mean"][1:6]-apply(ps,2,mean),3)
#   a.h1   v.h1   z.h1  sv.h1  sz.h1  t0.h1 
# -0.001 -0.021  0.002  0.103  0.003  0.000 

# OK for all, similar to individual
round(gest$statistics[,"Mean"][7:10],3)
#  a.h2  v.h2  z.h2 t0.h2 
# 0.194 0.191 0.091 0.056 
round(gest$statistics[,"Mean"][7:10]-apply(ps,2,sd)[-c(4,5)],3)
#   a.h2   v.h2   z.h2  t0.h2 
#  0.005 -0.003  0.002  0.002 

gpp <- h.post.predict.dmc(hsamples1.group)
tmp <- lapply(gpp, function(x){plot.pp.dmc(x, style="cdf") })

save_data (raw.data,hsamples.good,hsamples1.good,hsamples,hsamples1,
           hsamples.group,hsamples1.group,gpp,hest,gest,
           file="dmc_4_7group_random.RData")
