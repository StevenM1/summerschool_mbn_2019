##################  DMC Lesson 3: Sampling


### Lesson 3.4:  Sampling and assessing a single DDM subject


rm(list=ls()) 
# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("DDM","ddm.R")

# load_data ("dmc_3_4.RData")

p.map <- list(a="1",v="1",z="1",d="1",sz="1",sv="1",t0="1",st0="1")
const <- c(st0=0,d=0)

# Design is the same simple one used in previous examples.
model <- model.dmc(p.map,constants=const,
                   match.map=list(M=list(s1="r1",s2="r2")),
                   factors=list(S=c("s1","s2")),
                   responses=c("r1","r2"),
                   type="rd")

p.vector  <- c(a=1,v=1,z=0.5,sv=1,sz=0.2,t0=.15)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)

# Accuracy around 70%
par(mfrow=c(1,2))
plot.cell.density(data.cell=data.model[data.model$S=="s1",],C="r1",xlim=c(0,2),main="s1")
plot.cell.density(data.cell=data.model[data.model$S=="s2",],C="r2",xlim=c(0,2),main="s2")


# Profiles are broad for variability parameters and a little
# off for the sz variability parameter in particular. 
# NOTE: This can be very slow!
par(mfrow=c(2,3)); ylim = c(200,475)
profile.dmc("a",.9,1.1,p.vector,data.model,ylim=ylim)
profile.dmc("v",.5,1.5,p.vector,data.model,ylim=ylim)
profile.dmc("z",.45,.55,p.vector,data.model,ylim=ylim)
profile.dmc("sv",.1,2,p.vector,data.model,ylim=ylim)
profile.dmc("sz",.1,.45,p.vector,data.model,ylim=ylim)
profile.dmc("t0",.14,.16,p.vector,data.model,ylim=ylim)

# Note if you look at ddm.R you will see a function "bad" in the 
# likelihood.dmc function (see also tutorial 1_5). It enforces these bounds by 
# setting the likelihood to zero when these constraints are breached.
# You will note that there are other bounds enforced, including some a
# little inside the theoretically allowable bounds in order to avoid numerical
# problems in the C code that performs the numerical integration to get the
# ddm likelihood. It can sometimes have problems with the extreme estimates that
# can occur during sampling. Note that the set implemented in "bad" may 
# sometimes need augmentation by editing the ddm.R script. 

# Here we use beta(1,1) (i.e., uniform) priors on the z and sz parameters
# (so they are bounded in the 0-1 range, recall both are defined RELATIVE to a), 
# the same for t0 (to keep it greater than zero, which it must be by definition, 
# and < 1 to impose prior knowledge that values larger than 1 second are very
# unlikley in most applications) and finally truncated normal priors
# on the remaining parameters (to enforce a > 0, sv > 0 and t0 > 0) with broad 
# bounds otherwise but not unbounded in order to avoid numerical issues, e.g.
# very large a or very large or small v). These choices might vary in other 
# applications, and interact with the likelihood restrictions described above. 

p.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","beta","tnorm","beta","beta"),
  p1=c(a=1,v=0,z=1,sv=1,sz=1,t0=1),                           
  p2=c(a=1,v=2,z=1,sv=1,sz=1,t0=1),
  lower=c(0,-5,NA,0,NA,NA),
  upper=c(2, 5,NA,2,NA,NA)
)
par(mfcol=c(2,3)); for (i in names(p.prior)) plot.prior(i,p.prior)

# Parameters of the DDM are more strongly correlated than those of the LNR
# hence longer burnin and more stuck chains are to be expected, so try a
# similar longer burnin as for LBA.
samples <- samples.dmc(nmc=400,p.prior,data.model)
samples <- run.dmc(samples, report = 25, cores=4,p.migrate=.05)
plot.dmc(samples,pll.chain=TRUE,layout=c(1,3))
plot.dmc(samples)
plot.dmc(samples,pll.chain=TRUE,start=300,layout=c(1,3))
plot.dmc(samples,start=300)


# Looks good, get 500 without migration. 
samples1 <- run.dmc(samples.dmc(nmc=500,samples=samples),cores=4,report=50)
plot.dmc(samples1,pll.chain=TRUE)
plot.dmc(samples1)

# R-hat shows whole series not fully converged, but it doesn't look to be moving
gelman.diag.dmc(samples1)
#    Point est. Upper C.I.
# a        1.10       1.15
# v        1.10       1.16
# z        1.10       1.15
# sz       1.22       1.42
# sv       1.12       1.18
# t0       1.12       1.18
# 
# Multivariate psrf
# 
# 1.19

# Effective size is still fairly small, so probably just need a longer series.
effectiveSize.dmc(samples1)
#   a   v   z  sz  sv  t0 
# 339 314 344 292 326 291 


samples2 <- run.dmc(samples.dmc(nmc=500,samples=samples1,add=TRUE), 
                    cores=4,report=25)
plot.dmc(samples2)

# Now close to converged, sz is a little tardy. 
gelman.diag.dmc(samples2)
#    Point est. Upper C.I.
# a        1.04       1.06
# v        1.04       1.07
# z        1.04       1.05
# sz       1.14       1.23
# sv       1.04       1.06
# t0       1.05       1.08
# 
# Multivariate psrf
# 
# 1.14

effectiveSize.dmc(samples2)
#   a   v   z  sz  sv  t0 
# 594 573 625 486 562 526 

# Quite autocorrelted
acf.dmc(samples2)

# Looks like we might need an even longer series. Lets try this with automatic 
# instead so we don’t have to figure it out manually.
samples.auto <- samples.dmc(nmc=100,p.prior,data=samples$data)

# Takes 4 cycles of 100 to get rid of stuck chains
samples1.auto <- run.unstuck.dmc(samples.auto,p.migrate=.05,verbose=TRUE)
# 10  20  30  40  50  60  70  80  90  100  
# ...
# Bad chains: 8 7 1 13 5 4 6 9 18
# ...
# Bad chains: 13 18 6 1 12 16 11 14 5
# ...
# Bad chains: 18 15
# ...
# Bad chains: None

# Then sample to convergence and a minimum effectiveSize of 500. Takes 1300 to
# get there and near the end Rhat is not monotonic decreasing.
samples2.auto <- run.converge.dmc(samples1.auto,minN=500,nmc=50,max.try=50,verbose=TRUE)
# 110  120  130  140  150  
# [1] "N = 150 Multivariate psrf achieved = 1.33899146633459"
# ...
# [1] "N = 1100 Multivariate psrf achieved = 1.10805186781507"
# 1110  1120  1130  1140  1150  
# [1] "N = 1150 Multivariate psrf achieved = 1.11406650045681"
# 1160  1170  1180  1190  1200  
# [1] "N = 1200 Multivariate psrf achieved = 1.11669803801465"
# 1210  1220  1230  1240  1250  
# [1] "N = 1250 Multivariate psrf achieved = 1.10591229677169"
# 1260  1270  1280  1290  1300  
# [1] "Final multivariate psrf = 1.09652788931662"
# Effective sample size
#   a   v   z  sz  sv  t0 
# 716 669 772 508 732 591 

# Looking in more detail we see sz is still the highest. However, this is
# probably OK, Gelman et al. (2004) recommendation in the BDA text was that 
# values < 1.2 are OK.
gelman.diag.dmc(samples2.auto)
#    Point est. Upper C.I.
# a        1.03       1.05
# v        1.03       1.05
# z        1.03       1.04
# sz       1.11       1.17
# sv       1.03       1.05
# t0       1.06       1.09
# 
# Multivariate psrf
# 
# 1.1


# Good number of samples, but strong autocorrelation
acf.dmc(samples2.auto)

# Confirm that no chains are stuck
pick.stuck.dmc(samples2.auto,cut=10,verbose=TRUE)

# Tabled estimates show good recovery, although sz in particular has a wide CI.
summary.dmc(samples2.auto)
#       2.5%    25%    50%    75%  97.5%
# a  0.99093 1.0014 1.0080 1.0137 1.0260
# v  0.95276 0.9972 1.0225 1.0480 1.0965
# z  0.49555 0.4980 0.4993 0.5007 0.5031
# sz 0.06808 0.2363 0.2806 0.3147 0.3673
# sv 0.83226 0.9890 1.0625 1.1372 1.2720
# t0 0.14657 0.1497 0.1511 0.1522 0.1539

# Good updating although sv and sz less so than the other parameters
plot.dmc(samples2.auto,p.prior=p.prior)


# Strong correlations but not as strong as LBA
pairs.dmc(samples2.auto)

# Sample posterior predictive to check fit
pp <- post.predict.dmc(samples2.auto)

# Good fit to pdf
plot.pp.dmc(pp)

# And cdf
plot.pp.dmc(pp,"cdf")

# Posterior predictive p value for robust skew (rather slow!)
ppp.dmc(samples2.auto,plot.density=TRUE,
        fun=function(data) diff(diff(quantile(data$RT,probs=c(.25,.5,.75)))) )
# [1] 0.278

# save_data (p.vector,data,samples,samples1,samples2,
#            samples1.auto,samples2.auto,pp,file="dmc_3_4.RData")





