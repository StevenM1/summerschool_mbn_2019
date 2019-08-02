##################  DMC Lesson 6: More Models

### Lesson 6.1:  Sampling and assessing a single LBA subject with go failure
# "go failure" means that on some proportion of trials participants fail to 
# respond: we estimate a parameter "gf" to account for this.

rm(list=ls()) 

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("LBA-GoFailure","lba_Bgf.R")

# load_data ("dmc_6_1.RData")

model <- model.dmc(
  p.map=list(A="1",B="1",mean_v="M",sd_v="M",t0="1",st0="1",gf="1"),
  constants=c(st0=0,sd_v.false=1),match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="normGF")
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"            "B"            "mean_v.true"  "mean_v.false" "t0"           "gf"          
# 
# Constants are (see attr(,"constants") ):
#  st0 sd_v 
#    0    1 
# 
# Model type = norm (posdrift= TRUE ) 


# Simulate some data, with around 65% accuracy
p.vector  <- c(A=.25,B=.75,mean_v.true=1.5,mean_v.false=.25,sd_v.true=.6,t0=.2,gf=-1)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)

# Data distributions similar to LNR model
par(mfrow=c(1,2))
plot.cell.density(data.cell=data.model[data.model$S=="s1",],C="r1",xlim=c(0,2))
plot.cell.density(data.cell=data.model[data.model$S=="s2",],C="r2",xlim=c(0,2))

head(likelihood.dmc(p.vector,data.model))

# Check profiles, they all look OK
par(mfrow=c(2,4))
profile.dmc("A",      .1,  .5,p.vector,data.model,ylim=NA)
profile.dmc("B",      .5,  1.5,p.vector,data.model,ylim=NA)
profile.dmc("mean_v.true",  1,  2,p.vector,data.model,ylim=NA)
profile.dmc("mean_v.false",0, .5,p.vector,data.model,ylim=NA)
profile.dmc("sd_v.true",  .2,  1,p.vector,data.model,ylim=NA)
profile.dmc("t0",     .05,.35,p.vector,data.model,ylim=NA)
profile.dmc("gf",     -1.5,-.5,p.vector,data.model,ylim=NA)



p.prior <- prior.p.dmc(
  dists = rep("tnorm",7),
  p1=c(A=.3,B=.3,mean_v.true=1,mean_v.false=0,sd_v.true=1,t0=1,gf=-1.5),                           
  p2=c(1,1,3,3,1,1,2),lower=c(0,0,NA,NA,0,.1,NA),upper=c(NA,NA,NA,NA,NA,1,NA)
)
par(mfcol=c(2,4)); for (i in names(p.prior)) plot.prior(i,p.prior)



# Try estimation
samples <- samples.dmc(nmc=100,p.prior,data.model)
samples1 <- run.unstuck.dmc(samples,p.migrate=.05,verbose=TRUE,cores=21)

samples2 <- run.converge.dmc(samples.dmc(samples=samples1,nmc=100),
                             nmc=50,verbose=TRUE,max.try=10,cores=21)
# MPSRF:  1.09314162172984 
# [1] "Final multivariate psrf = 1.09314162172984"
# Effective sample size
#            A            B  mean_v.true mean_v.false    sd_v.true           t0           gf 
#           69           54           66           37           44          127           40

# High autocorrelation!          
samples3 <- run.dmc(samples.dmc(samples=samples2,nmc=100,thin=25),report=1,cores=21)    

# Not stable until 20
samples4 <- samples.dmc(samples=samples3,remove=1:20,add=TRUE,nmc=0)

effectiveSize.dmc(samples4)
#            A            B  mean_v.true mean_v.false    sd_v.true           t0           gf 
#          712          621          748          750          875          669         1095 

plot.dmc(samples4,pll.chain=TRUE)
plot.dmc(samples4,layout=c(2,4))
plot.dmc(samples4,layout=c(2,4),p.prior=p.prior)

# Usual recovery, poorer for of A and B
check.recovery.dmc(samples4,p.vector)
#                    A    B mean_v.true mean_v.false sd_v.true    t0    gf
# True            0.25 0.75        1.50         0.25      0.60  0.20 -1.00
# 2.5% Estimate   0.01 0.64        1.46         0.15      0.59  0.16 -1.04
# 50% Estimate    0.15 0.89        1.61         0.41      0.63  0.18 -1.02
# 97.5% Estimate  0.35 1.03        1.70         0.56      0.66  0.23 -1.00
# Median-True    -0.10 0.14        0.11         0.16      0.03 -0.02 -0.02


# Usual high correlations, but gf good.
pairs.dmc(samples4)


######### SIMPLE RT = 100% accuracy case (with two stimuli)

# Model where B depends on S as s1 and s2 in different blocks
modelS <- model.dmc(
  p.map=list(A="1",B="S",mean_v="M",sd_v="M",t0="1",st0="1",gf="1"),
  constants=c(st0=0,sd_v.false=1,mean_v.false=-100),
  match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="normGF")

# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"           "B.s1"        "B.s2"        "mean_v.true" "sd_v.true"   "t0"          "gf"         
# 
# Constants are (see attr(,"constants") ):
#          st0   sd_v.false mean_v.false 
#            0            1         -100 
# 
# Model type = normGF 

# Simulate some data, with bias to s1
p.vector  <- c(A=.25,B.s1=.5,B.s2=1,mean_v.true=1.5,sd_v.true=.6,t0=.2,gf=-1)
data.model <- data.model.dmc(simulate.dmc(p.vector,modelS,n=1e4),modelS)

# Data distributions similar to LNR model
par(mfrow=c(1,2))
plot.cell.density(data.cell=data.model[data.model$S=="s1",],C="r1",xlim=c(0,2))
plot.cell.density(data.cell=data.model[data.model$S=="s2",],C="r2",xlim=c(0,2))

head(likelihood.dmc(p.vector,data.model))

# Check profiles, they all look OK
par(mfrow=c(2,4))
profile.dmc("A",      .1,  .5,p.vector,data.model,ylim=NA)
profile.dmc("B.s1",      .25,  .75,p.vector,data.model,ylim=NA)
profile.dmc("B.s2",      .75,  1.25,p.vector,data.model,ylim=NA)
profile.dmc("mean_v.true",  1,  2,p.vector,data.model,ylim=NA)
profile.dmc("sd_v.true",  .2,  1,p.vector,data.model,ylim=NA)
profile.dmc("t0",     .05,.35,p.vector,data.model,ylim=NA)
profile.dmc("gf",     -1.5,-.5,p.vector,data.model,ylim=NA)



p.prior <- prior.p.dmc(
  dists = rep("tnorm",7),
  p1=c(A=.3,B.s1=.3,B.s2=.3,mean_v.true=1,sd_v.true=1,t0=1,gf=-1.5),                           
  p2=c(1,1,1,3,1,1,2),lower=c(0,0,0,NA,0,.1,NA),upper=c(NA,NA,NA,NA,NA,1,NA)
)
par(mfcol=c(2,4)); for (i in names(p.prior)) plot.prior(i,p.prior)



# Try estimation
samplesS <- samples.dmc(nmc=100,p.prior,data.model)
samples1S <- run.unstuck.dmc(samplesS,p.migrate=.05,verbose=TRUE,cores=21)

samples2S <- run.converge.dmc(samples.dmc(samples=samples1S,nmc=100),
                              nmc=50,verbose=TRUE,max.try=10,cores=21)
# [1] "Final multivariate psrf = 1.09719385696495"
# Effective sample size
#           A        B.s1        B.s2 mean_v.true   sd_v.true          t0          gf 
#         302         250         266         289         299         287         295

# Not so high autocorrelation!          
samples3S <- run.dmc(samples.dmc(samples=samples2S,nmc=100,thin=10),report=1,cores=21)    

# Not stable until 20, small effective sample size, get a bigger sample
samples4S <- run.dmc(samples.dmc(samples=samples3S,nmc=100,thin=25),report=1,cores=21)    

effectiveSize.dmc(samples4S)
#           A        B.s1        B.s2 mean_v.true   sd_v.true          t0          gf 
#         735         499         514         536         532         610         798 

plot.dmc(samples4S,pll.chain=TRUE)
plot.dmc(samples4S,layout=c(2,4))
# Weak updating of B
plot.dmc(samples4S,layout=c(2,4),p.prior=p.prior)

# Very poor recovery of everything except t0 and gf
check.recovery.dmc(samples4S,p.vector)
#                   A B.s1 B.s2 mean_v.true sd_v.true    t0    gf
# True           0.25 0.50 1.00        1.50      0.60  0.20 -1.00
# 2.5% Estimate  0.17 0.49 0.96        1.46      0.57  0.17 -1.01
# 50% Estimate   0.40 0.99 1.93        2.83      1.13  0.19 -0.99
# 97.5% Estimate 0.76 1.63 3.09        4.53      1.83  0.22 -0.97
# Median-True    0.15 0.49 0.93        1.33      0.53 -0.01  0.01

# Extremely strong B and v correlations
pairs.dmc(samples4S)

######### SIMPLE RT = does A=0 help?

# Model where B depends on S as s1 and s2 in different blocks
modelS1 <- model.dmc(
  p.map=list(A="1",B="S",mean_v="M",sd_v="M",t0="1",st0="1",gf="1"),
  constants=c(st0=0,sd_v.false=1,mean_v.false=-100,A=0),
  match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="normGF")

# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "B.s1"        "B.s2"        "mean_v.true" "sd_v.true"   "t0"          "gf"         
# 
# Constants are (see attr(,"constants") ):
#          st0   sd_v.false mean_v.false            A 
#            0            1         -100            0 
# 
# Model type = normGF 

# Simulate some data, with bias to s1
p.vector  <- c(B.s1=.5,B.s2=1,mean_v.true=1.5,sd_v.true=.6,t0=.2,gf=-1)
data.model <- data.model.dmc(simulate.dmc(p.vector,modelS1,n=1e4),modelS1)

# Data distributions similar to LNR model
par(mfrow=c(1,2))
plot.cell.density(data.cell=data.model[data.model$S=="s1",],C="r1",xlim=c(0,2))
plot.cell.density(data.cell=data.model[data.model$S=="s2",],C="r2",xlim=c(0,2))

head(likelihood.dmc(p.vector,data.model))

# Check profiles, they all look OK
par(mfrow=c(2,3))
profile.dmc("B.s1",      .25,  .75,p.vector,data.model,ylim=NA)
profile.dmc("B.s2",      .75,  1.25,p.vector,data.model,ylim=NA)
profile.dmc("mean_v.true",  1,  2,p.vector,data.model,ylim=NA)
profile.dmc("sd_v.true",  .2,  1,p.vector,data.model,ylim=NA)
profile.dmc("t0",     .05,.35,p.vector,data.model,ylim=NA)
profile.dmc("gf",     -1.5,-.5,p.vector,data.model,ylim=NA)



p.prior <- prior.p.dmc(
  dists = rep("tnorm",6),
  p1=c(B.s1=.3,B.s2=.3,mean_v.true=1,sd_v.true=1,t0=1,gf=-1.5),                           
  p2=c(1,1,3,1,1,2),lower=c(0,0,NA,0,.1,NA),upper=c(NA,NA,NA,NA,1,NA)
)
par(mfcol=c(2,3)); for (i in names(p.prior)) plot.prior(i,p.prior)



# Try estimation
samplesSnoA <- samples.dmc(nmc=100,p.prior,data.model)
samples1SnoA <- run.unstuck.dmc(samplesSnoA,p.migrate=.05,verbose=TRUE,cores=18)

samples2SnoA <- run.converge.dmc(samples.dmc(samples=samples1SnoA,nmc=100),
                                 nmc=50,verbose=TRUE,max.try=10,cores=18)
[1] "Final multivariate psrf = 1.09289698304832"
# Effective sample size
#        B.s1        B.s2 mean_v.true   sd_v.true          t0          gf 
#         263         264         262         263         297         389 

# Not so high autocorrelation!          
samples3SnoA <- run.dmc(samples.dmc(samples=samples2SnoA,nmc=100,thin=10),report=1,cores=18)    

# Get a bigger sample
samples4SnoA <- run.dmc(samples.dmc(samples=samples3SnoA,nmc=100,thin=25),report=1,cores=18)    

effectiveSize.dmc(samples4SnoA)
#        B.s1        B.s2 mean_v.true   sd_v.true          t0          gf 
#         757         773         773         763        1054        1281 

plot.dmc(samples4SnoA,pll.chain=TRUE)
plot.dmc(samples4SnoA,layout=c(2,3))
# Still weak updating of B
plot.dmc(samples4SnoA,layout=c(2,3),p.prior=p.prior)

# Still very poor recovery of everything except t0 and gf
check.recovery.dmc(samples4SnoA,p.vector)
#                B.s1 B.s2 mean_v.true sd_v.true   t0    gf
# True           0.50 1.00        1.50      0.60 0.20 -1.00
# 2.5% Estimate  0.35 0.70        1.05      0.42 0.19 -1.01
# 50% Estimate   0.83 1.67        2.53      1.00 0.20 -0.99
# 97.5% Estimate 1.41 2.83        4.26      1.69 0.21 -0.97
# Median-True    0.33 0.67        1.03      0.40 0.00  0.01

# Now essentiall perfect B and v correlations
pairs.dmc(samples4SnoA)

# save_data (samples3,samples4S,samples4SnoA,file="dmc_6_1.RData")
