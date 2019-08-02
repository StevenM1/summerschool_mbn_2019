##################  DMC Lesson 3: Sampling

### Lesson lba1v:  Sampling and assessing a single LBA subject with a
###                one-dimensional v model

rm(list=ls()) 

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("LBA","lba_Bdv.R")

# load_data ("lba_Bdv_lesson.RData")

# NB1: New parameter d is "drive" sum of match and mismatch inputs.
# NB2: For the LBA one dimensional model no need to fix any parameters for scaling
#      purposes, but in order to enforce make v.false = 1-v.true need to set all
#      v.false to a constant -Inf (for more complicated designs this is true for
#      EVERY v.false parameter)
model <- model.dmc(
  p.map <- list(A="1",B="1",d="1",mean_v="M",sd_v="S",t0="1",st0="1"),
  constants=c(st0=0,mean_v.false=-Inf,sd_v.s1=0.4),
  match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="norm")
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"           "B"           "d"           "mean_v.true" "sd_v.s2"     "t0"         
# 
# Constants are (see attr(,"constants") ):
#          st0 mean_v.false      sd_v.s1 
#            0         -Inf            1 
# 
# Model type = norm (posdrift= TRUE ) 



# Simulate some data, with around 70% accuracy
p.vector  <- c(A=.25,B=.35,d=1,mean_v.true=.75,sd_v.s2=0.6,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)

# About 70% accuracy for s2 and 80% for s1 (lower sv)
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
  dists = c("tnorm","tnorm","tnorm","tnorm","tnorm","beta"),
  p1=c(A=.3,B=.3,d=.5,mean_v.true=1,sd_v.s2=.3,t0=1),                           
  p2=c(1,1,3,3,1,1),lower=c(0,0,0,NA,0,.1),upper=c(NA,NA,NA,NA,NA,1)
)
# par(mfcol=c(2,3)); for (i in names(p.prior)) plot.prior(i,p.prior)


samples <- samples.dmc(nmc=400,p.prior,data.model)
samples <- run.dmc(samples, report = 10, cores=18,p.migrate=.05)


plot.dmc(samples,pll.chain=TRUE)
# Looks stable by 350
plot.dmc(samples,pll.chain=TRUE,start=350)
plot.dmc(samples,start=350,layout=c(2,3))

# Looks like burnt in, so get a longer run without migration and add some thinning
# as looks fairly autocorrelated
samples1 <- run.dmc(samples.dmc(nmc=100,samples=samples,thin=5),cores=18)

plot.dmc(samples1,pll.chain=TRUE)
plot.dmc(samples1,layout=c(2,3))
# Priors nicely updated
plot.dmc(samples1,layout=c(2,3),p.prior=p.prior)

# R-hat shows whole series is close to converged
gelman.diag.dmc(samples1)
#             Point est. Upper C.I.
# A                 1.02       1.03
# B                 1.02       1.04
# d                 1.02       1.03
# mean_v.true       1.02       1.03
# sd_v.s2           1.02       1.04
# t0                1.02       1.04
# 
# Multivariate psrf
# 
# 1.05

effectiveSize.dmc(samples1)
#           A           B           d mean_v.true     sd_v.s2          t0 
#         413         454         465         532         433         429 

# Could have moved thin up to 10 
acf.dmc(samples1)

# Excellent recovery
check.recovery.dmc(samples1,p.vector)
#                  A    B    d mean_v.true sd_v.s2   t0
# True           0.25 0.35 1.00        0.75    0.60 0.20
# 2.5% Estimate  0.19 0.32 0.97        0.74    0.59 0.18
# 50% Estimate   0.25 0.36 1.03        0.77    0.60 0.20
# 97.5% Estimate 0.30 0.41 1.09        0.80    0.62 0.21
# Median-True    0.00 0.01 0.03        0.02    0.00 0.00



# Strong correlations are very evident.
pairs.dmc(samples1)

# Sample posterior predictive to check fit
pp <- post.predict.dmc(samples1)

# Good fit to pdf
plot.pp.dmc(pp)

# And cdf
plot.pp.dmc(pp,"cdf")



save_data (samples,samples1,pp,file="lba_Bdv_lesson.RData")
