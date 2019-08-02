##################  DMC Lesson 3: Sampling

### Lesson 3.3_LBAup:  Sampling and assessing a single LBA subject, using an  
###   unbounded parameterization (except st0), all parameters sampled on log scale
###   Model like 3.3 but estimates sd_v.true (fix sd_v.false=1). Shows that 
###   even though this model should be identified it fails in both the unbounded
###   and original parameterization (see end of file).

rm(list=ls()) 

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("LBA","lba_Bup.R")
load("dmc_3_3_LBAup.RData")


p.map <- list(A="1",B="1",mean_v="M",sd_v="M",t0="1",st0="1")
const <- c(st0=0,sd_v.false=0) # NOTE: exp(sd_v.false)=1
model <- model.dmc(p.map,constants=const,match.map=list(M=list(s1=1,s2=2)),
                   factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="norm")


# Simulate some data, with around 65% accuracy
p.vector  <- c(A=log(.25),B=log(.35),mean_v.true=log(1),mean_v.false=log(.25),
               sd_v.true=log(1),t0=log(.2))

raw.data <- simulate.dmc(p.vector,model,n=1e4)
data.model <- data.model.dmc(raw.data,model)

# Data distributions similar to LNR model
par(mfrow=c(1,2))
plot.cell.density(data.cell=data.model[data.model$S=="s1",],C="r1",xlim=c(0,3))
plot.cell.density(data.cell=data.model[data.model$S=="s2",],C="r2",xlim=c(0,3))

# Unbounded, put in untrans so can look at priors on natural scale
p.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm","tnorm","tnorm"),
  p1=c(A=log(.3),B=log(.3),mean_v.true=1,mean_v.false=0,sd_v.true=log(1),t0=log(1)),                           
  p2=c(.5,.5,3,3,.5,.5),
  untrans=c(A="exp",B="exp",t0="exp",sd_v.true="exp")  
)
par(mfcol=c(2,3)); for (i in names(p.prior)) plot.prior(i,p.prior,natural=TRUE)

# Initially converges to local minimium before breaking out! 
samples <- samples.dmc(nmc=1000,p.prior,data.model)
samples <- run.dmc(samples, report = 25, cores=18,p.migrate=.05)
plot.dmc(samples,pll.chain=TRUE)
plot.dmc(samples,start=700)

# Looks like burnt in, so get a longer run without migration and a little thinning
samples1 <- run.dmc(samples.dmc(nmc=200,samples=samples,thin=5), 
                    cores=18,report=10)
plot.dmc(samples1)
plot.dmc(samples1,p.prior=p.prior)


# samples <- samples.dmc(nmc=500,samples=samples1,add=TRUE)
# samples <- run.dmc(samples, report = 25, cores=18,p.migrate=.05)


gelman.diag.dmc(samples1)

effectiveSize.dmc(samples1)

# Recovery is fairly awful!
check.recovery.dmc(samples1,p.vector)
#                    A     B mean_v.true mean_v.false sd_v.true    t0
# True           -1.39 -1.05        0.00        -1.39      0.00 -1.61
# 2.5% Estimate  -1.25 -1.41       -0.42        -1.56      0.04 -1.60
# 50% Estimate   -0.91 -1.22       -0.16        -1.21      0.15 -1.50
# 97.5% Estimate -0.65 -1.00        0.05        -0.91      0.27 -1.42
# Median-True     0.48 -0.17       -0.16         0.18      0.15  0.11


# Correlations not terribly strong except B and t0.
pairs.dmc(samples1)

# Sample posterior predictive to check fit
pp <- post.predict.dmc(samples1)

# Good fit to pdf
plot.pp.dmc(pp)

# And cdf
plot.pp.dmc(pp,"cdf")

# save(samples,samples1,pp,file="dmc_3_3_LBAup.RData")

# BUT NORMAL PARAMETERIZAITON ALSO FAILS

load_model ("LBA","lba_B.R")


p.map <- list(A="1",B="1",mean_v="M",sd_v="M",t0="1",st0="1")
const <- c(st0=0,sd_v.false=1) 
modelN <- model.dmc(p.map,constants=const,match.map=list(M=list(s1=1,s2=2)),
                    factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="norm")


# Back to original numbers
p.vectorN  <- c(A=.25,B=.35,mean_v.true=1,mean_v.false=.25,
                sd_v.true=1,t0=.2)

data.modelN <- data.model.dmc(raw.data,modelN)


p.priorN <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm","tnorm","tnorm"),
  p1=c(A=.3,B=.3,mean_v.true=1,mean_v.false=0,sd_v.true=1,t0=1),                           
  p2=c(1,1,3,3,1,1),lower=c(0,0,NA,NA,0,0),upper=c(NA,NA,NA,NA,NA,NA)
)
par(mfcol=c(2,3)); for (i in names(p.priorN)) plot.prior(i,p.priorN)

# Initially converges to local minimium before breaking out! 
samplesN <- samples.dmc(nmc=1000,p.priorN,data.modelN)
samplesN <- run.dmc(samplesN, report = 25, cores=18,p.migrate=.05)
plot.dmc(samplesN,pll.chain=TRUE)
plot.dmc(samplesN,start=700)

# Looks like burnt in, so get a longer run without migration and a little thinning
samples1N <- run.dmc(samples.dmc(nmc=200,samples=samplesN,thin=5), 
                     cores=18,report=10)
plot.dmc(samples1N)
plot.dmc(samples1N,p.prior=p.priorN)

gelman.diag.dmc(samples1N)

effectiveSize.dmc(samples1N)

# Recovery is fairly awful!
check.recovery.dmc(samples1N,p.vectorN)
#                    A    B mean_v.true mean_v.false sd_v.true    t0
# True            0.25 0.35        1.00         0.25      1.00  0.20
# 2.5% Estimate   0.01 0.26       -0.27        -1.82      0.91  0.17
# 50% Estimate    0.20 0.37        0.02        -1.37      0.99  0.19
# 97.5% Estimate  0.38 0.54        0.25        -1.01      1.13  0.23
# Median-True    -0.05 0.02       -0.98        -1.62     -0.01 -0.01

save(samples,samples1,pp,samplesN,samples1N,file="dmc_3_3_LBAup.RData")

