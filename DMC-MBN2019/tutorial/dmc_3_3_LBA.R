##################  DMC Lesson 3: Sampling

### Lesson 3.3:  Sampling and assessing a single LBA subject
rm(list=ls()) 

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("LBA","lba_B.R")

# load_data ("dmc_3_3.RData")

# NB: For the LBA one of B, mean_v or sd_v must be fixed in at least one
#     cell of the design. In this case we fix sd_v for both cells
p.map <- list(A="1",B="1",mean_v="M",sd_v="1",t0="1",st0="1")
const <- c(st0=0,sd_v=1)

# Same simple design as for previous LNR examples
model <- model.dmc(p.map,constants=const,match.map=list(M=list(s1=1,s2=2)),
                   factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="norm")


# Simulate some data, with around 65% accuracy
p.vector  <- c(A=.25,B=.35,mean_v.true=1,mean_v.false=.25,t0=.2)

data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)

# Data distributions similar to LNR model
par(mfrow=c(1,2))
plot.cell.density(data.cell=data.model[data.model$S=="s1",],C="r1",xlim=c(0,2))
plot.cell.density(data.cell=data.model[data.model$S=="s2",],C="r2",xlim=c(0,2))

# Slow errors
crct <- (data.model$S=="s1" & data.model$R=="r1") |
  (data.model$S=="s2" & data.model$R=="r2")
round(tapply(data.model$RT,list(data.model$S,C=crct),mean),2)

# Give t0 a uniform prior from 0.1-1s, other priors normal, truncated
# below for A and B as they must be positive, unbounded for the v
p.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm","beta"),
  p1=c(A=.3,B=.3,mean_v.true=1,mean_v.false=0,t0=1),                           
  p2=c(1,1,3,3,1),lower=c(0,0,NA,NA,.1),upper=c(NA,NA,NA,NA,1)
)
par(mfcol=c(2,3)); for (i in names(p.prior)) plot.prior(i,p.prior)


# Parameters of the LBA are more strongly correlated than those of the LNR
# hence longer burnin and more stuck chains are to be expected, so try a longer
# burnin run than for LNR.
samples <- samples.dmc(nmc=400,p.prior,data.model)
samples <- run.dmc(samples, report = 25, cores=4,p.migrate=.05)
plot.dmc(samples,pll.chain=TRUE)
plot.dmc(samples)
plot.dmc(samples,start=300)

# Looks like burnt in, so get a longer run without migration
samples1 <- run.dmc(samples.dmc(nmc=500,samples=samples), 
                    cores=4,report=25)
plot.dmc(samples1)

# R-hat shows whole series is close to converged
gelman.diag.dmc(samples1)
#              Point est. Upper C.I.
# A                  1.10       1.16
# B                  1.11       1.17
# mean_v.true        1.09       1.14
# mean_v.false       1.09       1.14
# t0                 1.11       1.18
# 
# Multivariate psrf
# 
# 1.13

# Add 500 more to see if that settles things down (longer series often help)
samples2 <- run.dmc(samples.dmc(nmc=500,samples=samples1,add=TRUE), 
                    cores=4,report=50)
plot.dmc(samples2)

# Now reports converged
gelman.diag.dmc(samples2)
#              Point est. Upper C.I.
# A                  1.06       1.10
# B                  1.06       1.09
# mean_v.true        1.04       1.06
# mean_v.false       1.04       1.06
# t0                 1.06       1.09
# 
# Multivariate psrf
# 
# 1.07

# Nothing stuck
# Confirm that no chains are stuck
pick.stuck.dmc(samples2,verbose=TRUE)


# Now have over 500 of each type
effectiveSize.dmc(samples2)
#            A            B  mean_v.true mean_v.false           t0 
#          641          524          593          557          543 

# Looking at autocorrelations they are quite long, persisting above 20 
acf.dmc(samples2)


# Tabled estimates confirm some inaccuracy in estimates of A, B, and 
# mean_v.false with mean_v and t0 parameters relatively well recovered.
summary.dmc(samples2)

# Excellent updating
plot.dmc(samples2,p.prior=p.prior)

# Strong correlations are very evident.
pairs.dmc(samples2,thin=10)

# NB: When there are more conditions in the experiment parameter correlations
#     reduce and parameters become easier to estimate (at least as long as 
#     there are a reasonable number of observations per condition) and 
#     constraints can be imposed across conditions (e.g., a bias manipulation
#     only affecting thresholds and a stimulus manipulation only affecting
#     rates).

# Sample posterior predictive to check fit
pp <- post.predict.dmc(samples2)

# Good fit to pdf
plot.pp.dmc(pp)

# And cdf
plot.pp.dmc(pp,"cdf")

# Posterior predictive p value for robust skew
ppp.dmc(samples2,plot.density=TRUE,
        fun=function(data) diff(diff(quantile(data$RT,probs=c(.25,.5,.75)))) )
# [1] 0.35   [indicative only]


### Extension 1: Replacing bad chains
# samples.dmc has a parameter replace.bad.chains that can either be
# assigned a vector of integers to indicate which chains to remove or
# set to TRUE in which case it applied heuristics to identify parameter
# chains whose median is either far away from the median of all chains 
# for that parameter or which are less variable than other chains. The
# chosen chains are then replaced by randomly choosen (with replacement)
# copies of the remaining chains. The heristics are enacted in the function
# pick.stuck.pars.dmc based on a parameter "cut" which multiples the
# interquartie range of the median or variability. The default value is 
# cut=5, which serves to pick out chains that are very clearly stuck on visual 
# inspection. There were no such cases in the exmaples above so to illusrate
# with samples1 we set cut=2.

samples1.cut2 <- samples.dmc(nmc=500,samples=samples1,replace.bad.chains = TRUE,cut=2)
# Replacing bad chains
# [1] 2

# The chain removed has consistently low values in middle interations.

par(mfrow=c(1,2))
plot(samples1$theta[2,"sz",],type="l",ylab="sz")
plot(samples1$theta[2,"t0",],type="l",ylab="t0")

# This case probably does not warrent sing replace.bad.chains, but the option
# can be useful to remove more obviously bad cases.

### Extension 2: Automatic sampling
samples.auto <- samples.dmc(nmc=100,p.prior,data=samples$data)
# Replacing bad chains
# [1] 2

# First repeatedly get a new set of nmc samples until there are no stuck chains.
# Here it takes 3 cycles of 100
samples1.auto <- run.unstuck.dmc(samples.auto,p.migrate=.05,verbose=TRUE)
# 10  20  30  40  50  60  70  80  90  100  
# ...
# Bad chains: 2 9 14 12 7 4 1
# ...
# Bad chains: 9 15 11 1 2 12
# ...
# Deviation of mean chain log-likelihood from median of means
#    15     4     1     5    10    14     8     6     7     2    12    13    11     3     9 
#  3.11  1.20  1.01  0.78  0.15  0.14  0.01  0.00 -0.05 -0.05 -0.32 -0.39 -0.42 -0.59 -0.73 
# Bad chains: None

# Then sample to convergence and a minimum effectiveSize of 500 
samples2.auto <- run.converge.dmc(samples1.auto,minN=500,nmc=50,verbose=TRUE,max.try=10)
# [1] "N = 150 Multivariate psrf achieved = 1.23301489360127"
# 160  170  180  190  200  
# ...
# 560  570  580  590  600  
# [1] "Final multivariate psrf = 1.12487413392696"
# Effective sample size
#            A            B  mean_v.true mean_v.false           t0 
#          349          323          362          355          328 

# Not quite there, so run some more, gets there in 7 cycles 
# (by defualt max.try=100).        
samples2.auto <- run.converge.dmc(samples2.auto,minN=500,nmc=50,verbose=TRUE)
# [1] "N = 650 Multivariate psrf achieved = 1.11106747394672"
# ...
# [1] "Final multivariate psrf = 1.06470818922721"
# Effective sample size
#            A            B  mean_v.true mean_v.false           t0 
#          549          517          597          582          508 


# save_data (samples,samples1,samples2,pp,
#            samples.auto,samples1.auto,samples2.auto,file="dmc_3_3.RData")



### Extension 2: Check out dmc_3_3_LBAup.R in the Extras directory. 

# There we fit a slightly more complex model with exactly the same 
# parameterization as here but just one more (sd_v) paramerer estimated, and show 
# that although it should in theory be identifiable in practice it is not with
# even the large data sample used here. This is shown to be true when parameters
# are estimated on two different scales. Although this failure is not generally
# true (i.e., with other parameter values recovery is fine) it illustrates the
# importance of checking parameter recovery for you model and design, preferably
# in a parameter region representative of your data.

