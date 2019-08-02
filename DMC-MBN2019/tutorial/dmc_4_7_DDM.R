##################  DMC Lesson 4: Multiple Subjects

### Lesson 4.7: Hierarchical DDM Model 2 x 2 with a rate effect on factor F.


rm(list=ls())

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("DDM","ddm.R") 

# load_data ("dmc_4_7_fixed.RData")
# load_data ("dmc_4_7_random.RData")

## Set up a DDM Model, rate effect of factor F
model <- model.dmc(
  p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
  match.map = list(M=list(s1="r1",s2="r2")),
  factors=list(S=c("s1","s2"),F=c("f1","f2")),
  constants = c(st0=0,d=0),
  responses = c("r1","r2"),
  type = "rd")

# Population distribution. 
pop.mean <- c(a=2,  v.f1=4, v.f2=3, z=0.5,sv=1, sz=0.3, t0=0.3)
pop.scale <-c(a=0.5,v.f1=.5,v.f2=.5,z=0.1,sv=.3,sz=0.1,t0=0.05)
pop.prior <- prior.p.dmc(
  dists = rep("tnorm",7),
  p1=pop.mean,                           
  p2=pop.scale,
  lower=c(0,-5, -5, 0, 0, 0, 0),
  upper=c(5, 7,  7, 1, 2, 1, 1)
)
##  Check population distributions
# par(mfcol=c(2,3)); for (i in names(pop.prior)) plot.prior(i,pop.prior)

# NOTE that we enforce bounds on all parameters to
# keep them either inside the definitional ranges of the parameters (e.g.,
# 0 < z < 1, 0 < sz < 1, t0 > 0, a > 0, sv > 0) or that avoid numerical 
# problems (from very large or very small values of v). We use these same 
# bounds on all priors througout the lesson, both at the data level and for 
# the location parameters at the hyper level.

# Simulate some data
raw.data <- h.simulate.dmc(model, p.prior = pop.prior, n = 250, ns = 40)
data.model <- data.model.dmc(raw.data, model)

# # Take a look at the first 10 subjects data
# par(mfrow=c(2,5)) # Row 1 = subjects 1..5, row2 = 6..10
# for (i in 1:40)
#   plot.cell.density(data.cell=data.model[[i]][data.model[[i]]$S=="s1",],C="r1")


# Take a look at parameters
ps <- round( attr(raw.data, "parameters"), 2); 
round(apply(ps,2,mean),2)
#    a v.f1 v.f2    z   sv   sz   t0 
# 1.93 4.06 3.09 0.50 0.95 0.30 0.31 
round(apply(ps,2,sd),2)
#    a v.f1 v.f2    z   sv   sz   t0 
# 0.48 0.50 0.56 0.12 0.30 0.10 0.04 

### FIT FIXED EFFECTS

# Specify a broader prior than the true population distribution.
# NOTE the bounds.
p.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm"),
  p1=pop.mean,                           
  p2=pop.scale*5,
  lower=c(0,-5, -5, 0, 0, 0, 0),
  upper=c(5, 7,  7, 1, 2, 1, 1)
)

# Start with thin of 5 assuming at least as bad as LNR  
samples  <- h.samples.dmc(nmc = 50, p.prior, data.model, thin = 5)

save(samples,file="4_7.RData")

samples  <- h.run.unstuck.dmc(samples, p.migrate = .05, cores = 20)
samples1  <- h.run.converge.dmc(h.samples.dmc(nmc=100, samples=samples,thin=5), 
                                cores = 40, nmc=40)


#Parameter chains look well converged and mixed
plot.dmc(samples1,density=FALSE,smooth=FALSE,subject=1,layout=c(2,4))

# All chains are well converged. 
gelman.diag.dmc(samples1)
#   28   31   12   21   25   38   35   39    3    5    8   34   10   13   23   22 
# 1.06 1.06 1.06 1.07 1.07 1.07 1.07 1.07 1.07 1.07 1.07 1.08 1.08 1.08 1.08 1.08 
#   24   16    2    4   37   40    1   19   33    9   36   11    7   14   27   32 
# 1.08 1.08 1.08 1.08 1.08 1.08 1.08 1.08 1.09 1.09 1.09 1.09 1.09 1.09 1.09 1.09 
#   18   15   26   17    6   30   20   29 
# 1.09 1.10 1.10 1.10 1.10 1.10 1.10 1.10 
# Mean
# [1] 1.08


# Effective sample size: good yield given 21 x 140 = 2940 samples
es <- effectiveSize.dmc(samples1)
round(apply(do.call(rbind,es),2,median))
#    a v.f1 v.f2    z   sz   sv   t0 
#  644  634  638  681  599  625  646 
# Subject 5 by far the worst
apply(do.call(rbind,es),2,min)
#    a v.f1 v.f2    z   sz   sv   t0 
#  451  469  407  521  430  452  484 

# Good fits?
pp <- h.post.predict.dmc(samples1) 
tmp <- lapply(pp, function(x){plot.pp.dmc(x, style="cdf") })

# Check against true values: fairly good parameter recovery.
est <- summary.dmc(samples1)

est.mean <- t(data.frame(lapply(est,function(x){
  x$statistics[,"Mean"]})))[,dimnames(ps)[[2]]]

# Pretty good
round(apply(est.mean,2,mean),3)
#     a  v.f1  v.f2     z    sv    sz    t0 
# 2.076 4.082 3.101 0.508 1.031 0.305 0.311 
round(apply(est.mean,2,mean)-apply(ps,2,mean),3)
#      a   v.f1   v.f2      z     sv     sz     t0 
#  0.064  0.071  0.064  0.000 -0.022 -0.013 -0.004 

# Scale for a, v and sv off by quite a bit
round(apply(est.mean,2,sd),3)
#     a  v.f1  v.f2     z    sv    sz    t0 
# 0.562 0.755 0.604 0.118 0.397 0.109 0.043 
round(apply(est.mean,2,sd)-apply(ps,2,sd),3)
#      a   v.f1   v.f2      z     sv     sz     t0 
#  0.101  0.183  0.097  0.000  0.069 -0.005 -0.002 

# save_data(raw.data,p.prior,samples,samples1,pp,est,file="dmc_4_7_fixed.RData")

### FIT RANDOM EFFECTS

# Use all truncated normal priors for locations
# NOTE the bounds.
mu.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm"),
  p1=pop.mean,                           
  p2=pop.scale*5,
  lower=c(0,-5, -5, 0, 0, 0, 0),
  upper=c(5, 7,  7, 1, 2, 1, 1)
)
# par(mfcol=c(2,3)); for (i in names(mu.prior)) plot.prior(i,mu.prior)

# As in the LBA example suppose we have little idea about individual variation 
# apart from being sure that standard deviations are not huge. Again, we use
# uniform priors for the hyper-scale parameters on the range 0-1. Note that
# similar results to those reported below were obtained with a range of
# different vague priors.
sigma.prior <- prior.p.dmc(
  dists = rep("beta", length(p.prior)),
  p1=c(a=1, v.f1=1,v.f2 = 1, z=1,sv=1, sz=1, t0=1),p2=c(1,1,1,1,1,1,1),
  upper=c(2,2,2,2,2,2,2)
)
# par(mfcol=c(2,3)); for (i in names(sigma.prior)) plot.prior(i,sigma.prior)

# Make a hyper-prior list
pp.prior <- list(mu.prior, sigma.prior)


hstart <- make.hstart(samples1)
theta1 <- make.theta1(samples1)

# Make a new samples object from the 40-subject data, longer and with more
# immediate thinning that LNR given likely more autocorrelation
hsamples <- h.samples.dmc(nmc = 50, p.prior = p.prior, pp.prior=pp.prior, 
                          data = data.model, thin = 5, hstart.prior=hstart,theta1=theta1)
hsamples <- h.run.unstuck.dmc(hsamples, cores=21, report=10, 
                              p.migrate = .05, h.p.migrate = .05)
hsamples1 <- h.run.converge.dmc(h.samples.dmc(nmc=100, samples=hsamples), 
                                nmc=50,cores=21)
# ...
# Try 4: N = 200, Effective N = NA, Hyper mpsrf = 1.05
# Subject mpsrf achieved (sorted):
#   15   10   36   31   34   37    6   35    2   12   20   16   28   23   30   17 
# 1.09 1.07 1.06 1.06 1.06 1.06 1.06 1.06 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 
#   33   22    3    8    7    9   39   14   32   38   21    1   13   24   40   27 
# 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.04 1.04 
#   11    4   25   19   18    5   29   26 
# 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04

# Hyper likelihoods stable 
plot.dmc(hsamples1,hyper=TRUE,pll.chain=TRUE)
# Chains look nice, although sv.h2 is very broad
plot.dmc(hsamples1,hyper=TRUE,layout=c(2,7))

gelman.diag.dmc(hsamples1,hyper=TRUE)
#        Point est. Upper C.I.
# a.h1          1.01       1.01
# v.f1.h1       1.00       1.01
# v.f2.h1       1.00       1.01
# z.h1          1.01       1.02
# sv.h1         1.05       1.07
# sz.h1         1.00       1.01
# t0.h1         1.01       1.01
# a.h2          1.01       1.01
# v.f1.h2       1.01       1.01
# v.f2.h2       1.01       1.01
# z.h2          1.01       1.02
# sv.h2         1.05       1.08
# sz.h2         1.02       1.04
# t0.h2         1.01       1.01
# 
# Multivariate psrf
# 
# 1.05


# Good yield at hyper given nominal 21*200=4200
effectiveSize.dmc(hsamples1,hyper=TRUE)
#    a.h1 v.f1.h1 v.f2.h1    z.h1   sv.h1   sz.h1   t0.h1    a.h2 v.f1.h2 v.f2.h2 
#    1972    2248    1887    1839    1156    1613    1928    1998    1858    1926 
#    z.h2   sv.h2   sz.h2   t0.h2 
#    1801     614    1304    2075 

# Less so at individual subject level  
es <- effectiveSize.dmc(hsamples1)
round(apply(data.frame(es),1,mean))        
#    a v.f1 v.f2    z   sz   sv   t0 
#  886  886  893  919  888  898  920 
round(apply(data.frame(es),1,min))        
#    a v.f1 v.f2    z   sz   sv   t0 
#  737  770  808  795  761  720  812 

# Parameter recovery
hest <- summary.dmc(hsamples1,hyper=TRUE)
# Mean recovery better than for indivdiual fits except for sv
round(hest$statistics[,"Mean"][1:7],3)
#    a.h1 v.f1.h1 v.f2.h1    z.h1   sv.h1   sz.h1   t0.h1 
#   2.007   3.948   2.994   0.508   0.955   0.312   0.312 
round(hest$statistics[,"Mean"][1:7]-apply(ps,2,mean),3)
#    a.h1 v.f1.h1 v.f2.h1    z.h1   sv.h1   sz.h1   t0.h1 
#  -0.005  -0.063  -0.043   0.000  -0.097  -0.006  -0.002 

# Reasonable updating but clearly worst for sv.
plot.dmc(hsamples1,hyper=TRUE,p.prior=pp.prior,layout=c(2,7))


# Scale improved on individual fits except for sv
round(hest$statistics[,"Mean"][8:14],3)
#    a.h2 v.f1.h2 v.f2.h2    z.h2   sv.h2   sz.h2   t0.h2 
#   0.498   0.644   0.499   0.122   0.725   0.109   0.045 
round(hest$statistics[,"Mean"][8:14]-apply(ps,2,sd),3)
#    a.h2 v.f1.h2 v.f2.h2    z.h2   sv.h2   sz.h2   t0.h2 
#   0.038   0.072  -0.008   0.004   0.397  -0.005   0.000 

# Fits good
hpp <- h.post.predict.dmc(hsamples1)
plot.pp.dmc(hpp, style="cdf")
tmp <- lapply(hpp, function(x){plot.pp.dmc(x, style="cdf") })


# save_data(raw.data,pp.prior,hsamples,hsamples1,hpp,hest,file="dmc_4_7_random.RData")

