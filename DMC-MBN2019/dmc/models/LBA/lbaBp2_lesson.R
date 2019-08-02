#  Sampling and assessing a single LBA subject with a B mixture

rm(list=ls()) 

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("LBA","lba_Bp2.R")

# load_data ("lbaBp2_lesson.RData")

# LBA uniform/normal version
model <- model.dmc(
  p.map=list(A="1",B1="1",B2="1",pb2="1",mean_v="M",sd_v="M",t0="1",st0="1"),
  constants=c(st0=0,sd_v.false=1),match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="normBp2")
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"            "B1"           "B2"           "pb2"          "mean_v.true"  "mean_v.false" "sd_v.true"   
# [8] "t0"          
# 
# Constants are (see attr(,"constants") ):
#        st0 sd_v.false 
#          0          1 
# 
# Model type = normBp2


# Simulate some data
p.vector  <- c(A=.25,B1=.25,B2=1,pb2=.5,mean_v.true=1.5,mean_v.false=.5,sd_v.true=.5,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)

# Data distributions, double threshold effect quite evident.
par(mfrow=c(1,2))
plot.cell.density(data.cell=data.model[data.model$S=="s1",],C="r1",xlim=c(0,3))
plot.cell.density(data.cell=data.model[data.model$S=="s2",],C="r2",xlim=c(0,3))

# Check profiles, they all look OK
par(mfrow=c(2,4))
profile.dmc("A",      .1,  .5,p.vector,data.model,ylim=NA)
profile.dmc("B1",      .15,  .5,p.vector,data.model,ylim=NA)
profile.dmc("B2",      .5,  1.5,p.vector,data.model,ylim=NA)
profile.dmc("pb2",      .25,  .75,p.vector,data.model,ylim=NA)
profile.dmc("mean_v.true",  1,  2,p.vector,data.model,ylim=NA)
profile.dmc("mean_v.false",0, 1,p.vector,data.model,ylim=NA)
profile.dmc("sd_v.true",  .2,  1,p.vector,data.model,ylim=NA)
profile.dmc("t0",     .05,.35,p.vector,data.model,ylim=NA)



p.prior <- prior.p.dmc(
  dists = rep("tnorm",8),
  p1=c(A=.3,B1=.3,B2=.3,pb2=.5,mean_v.true=1,mean_v.false=0,sd_v.true=1,t0=1),                           
  p2=c(1,1,1,1,3,3,1,1),lower=c(0,0,0,0,NA,NA,0,.1),upper=c(NA,NA,NA,1,NA,NA,NA,1)
)
# par(mfcol=c(2,4)); for (i in names(p.prior)) plot.prior(i,p.prior)



# Try estimation
samples <- samples.dmc(nmc=100,p.prior,data.model)
samples1 <- run.unstuck.dmc(samples,p.migrate=.05,verbose=TRUE,cores=24)

# Looks OK, but needs thinning.
get.thin(samples1)

samples2 <- run.dmc(samples.dmc(samples=samples1,nmc=100,thin=10),cores=24)

# NB On some replicaiton runs this didnt always work smoothly e.g., A appeared
#    stable then changed, causing a chain to get stuck again.


gelman.diag.dmc(samples2)
# Potential scale reduction factors:
# 
#              Point est. Upper C.I.
# A                  1.01       1.02
# B1                 1.02       1.03
# B2                 1.02       1.04
# pb2                1.01       1.01
# mean_v.true        1.03       1.04
# mean_v.false       1.02       1.04
# sd_v.true          1.02       1.04
# t0                 1.01       1.02
# 
# Multivariate psrf
# 
# 1.03

plot.dmc(samples2,pll.chain=TRUE)
plot.dmc(samples2,layout=c(2,4))
plot.dmc(samples2,layout=c(2,4),p.prior=p.prior)

check.recovery.dmc(samples2,p.vector)
#                   A    B1   B2  pb2 mean_v.true mean_v.false sd_v.true   t0
# True           0.25  0.25 1.00 0.50        1.50         0.50      0.50 0.20
# 2.5% Estimate  0.23  0.21 0.94 0.50        1.41         0.38      0.47 0.19
# 50% Estimate   0.26  0.24 1.00 0.50        1.49         0.51      0.50 0.20
# 97.5% Estimate 0.28  0.29 1.06 0.51        1.58         0.65      0.53 0.21
# Median-True    0.01 -0.01 0.00 0.00       -0.01         0.01      0.00 0.00

save_data (samples1,samples2,file="lbaBp2_lesson.RData")
