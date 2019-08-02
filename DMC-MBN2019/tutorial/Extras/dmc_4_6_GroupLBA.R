##################  DMC Lesson 4: Multiple Subjects

### Lesson 4.6groupLBA: Hierarchical LBA Model with a simple design and smaller
### sample size where full hierarchical fails (for the same design
### with twice the sample size hierarchical works). Convergence problems
### are helped a group parameters (assumed the same for all subjects) for a
### trial-to-trial variability parameter (A) but shown to not really fix the 
### problem completely and that adding more group parameters can sometimes also 
### casue more problems. Chain update mixing in DMC (default random.phi=TRUE and 
### random.theta=TRUE) solves the convergence problem, but parameter recovery 
### still difficult.

rm(list=ls())

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("LBA","lba_B.R")

# NB: This lesson does not do fixed effects.
# load_data ("dmc_4_6group.RData")

# Set up an LBA Model, exactly as 4_6_simple (i.e., like Lesson 3.3 except let 
# sd_v vary with match, fixing scale with sd_v.false).
model <- model.dmc(
  p.map     = list(A = "1", B = "1", mean_v = "M", sd_v = "M", t0 = "1", st0 = "1"),
  match.map = list(M = list(s1 = 1, s2 = 2)),
  factors   = list(S = c("s1", "s2")),
  constants = c(st0 = 0, sd_v.false=1),  
  responses=c("r1","r2"),
  type="norm")


# Setup sample as in 4_6_simple but half the size
pop.mean <- c(A=.4,B=.6,mean_v.true=1,mean_v.false=0,sd_v.true = .5,t0=.3)
pop.scale <-c(A=.1,B=.1,mean_v.true=.2,mean_v.false=.2,sd_v.true = .1,t0=.05)  
pop.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm","tnorm","tnorm"),
  p1=pop.mean,p2=pop.scale,
  lower=c(0,0,NA,NA,0,.1),upper=c(NA,NA,NA,NA,NA,1)
)
# #  Check population distributions
# par(mfcol=c(2,3)); for (i in names(pop.prior)) plot.prior(i,pop.prior)

raw.data <- h.simulate.dmc(model, p.prior = pop.prior, n = 250, ns = 40)
data.model <- data.model.dmc(raw.data, model)

# Take a look at parameters
ps <- attr(raw.data, "parameters") 
round(apply(ps,2,mean),3)
#            A            B  mean_v.true mean_v.false    sd_v.true           t0 
#        0.398        0.614        1.040       -0.032        0.485        0.271 
round(apply(ps,2,sd),3)
#            A            B  mean_v.true mean_v.false    sd_v.true           t0 
#        0.090        0.094        0.183        0.230        0.110        0.040 


### Full hierarchical fit.

p.prior <- prior.p.dmc(
  p1=c(A=1,B=1,mean_v.true=1,mean_v.false=1,sd_v.true=1,t0=1),p2=rep(1,6),
  lower=c(0,0,NA,NA,0,.1),upper=c(NA,NA,NA,NA,NA,1)
)
mu.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm","tnorm","tnorm"),
  p1=c(A=.5,B=.5,mean_v.true=1,mean_v.false=0,sd_v.true=.5,t0=0.5),                           
  p2=c(2,2,3,3,2,2),lower=c(0,0,NA,NA,0,.1),upper=c(NA,NA,NA,NA,NA,1)
)
# par(mfcol=c(2,3)); for (i in names(mu.prior)) plot.prior(i,mu.prior)
sigma.prior <- prior.p.dmc(
  dists = rep("gamma", length(p.prior)),
  p1=c(A=1,B=1,mean_v.true=1,mean_v.false=1,sd_v.true=1,t0=1),p2=rep(1,6)
)
# par(mfcol=c(2,2)); for (i in names(sigma.prior)) plot.prior(i,sigma.prior)
pp.prior <- list(mu.prior, sigma.prior)


hsamples <- h.samples.dmc(nmc = 50, p.prior=p.prior, data=data.model, 
                          pp.prior=pp.prior,thin=10)


### Fit with default setting of random.phi=TRUE, and random.theta=TRUE, which 
## mixes up links between hyper and lower level chains in update of hyper and 
### data parameters respectively. This is the default setting for DMC as it 
### fixes convergence problems identified in the next section, see also 4_7group.

hsamples.good <- hsamples

hsamples.good <- h.run.dmc(hsamples.good, cores=8, report=1,
                           p.migrate = .05, h.p.migrate = .05)
hsamples.good <- h.run.dmc(h.samples.dmc(nmc=150,samples=hsamples.good,add=TRUE),
                           cores=8,report=1,p.migrate = .05, h.p.migrate = .05)

# Gradual increase, A scale appears to be heading to zero!
plot.dmc(hsamples.good,smooth=FALSE,density=FALSE,pll.chain=TRUE,hyper=TRUE,start=20)
plot.dmc(hsamples.good,hyper=TRUE,layout=c(2,6),smooth=FALSE,density=FALSE,start=150)

# Sample without migration
hsamples1.good <- h.run.dmc(h.samples.dmc(nmc=100,samples=hsamples.good),
                            report=1,cores=8)
hsamples1.good <- h.run.dmc(h.samples.dmc(nmc=400,samples=hsamples1.good,add=TRUE),
                            report=1,cores=8)

# Initial drop due to migration over-fitting. Good from 100
plot.dmc(hsamples1.good,smooth=FALSE,density=FALSE,pll.chain=TRUE,hyper=TRUE)
# Some A chains have zero variance (also one mean_v.false)
plot.dmc(hsamples1.good,hyper=TRUE,layout=c(2,6),smooth=FALSE,density=FALSE,start=100)

# Good mixing
gelman.diag.dmc(hsamples1.good,hyper=TRUE,start=100)
#                 Point est. Upper C.I.
# A.h1                     1          1
# B.h1                     1          1
# mean_v.true.h1           1          1
# mean_v.false.h1          1          1
# sd_v.true.h1             1          1
# t0.h1                    1          1
# A.h2                     1          1
# B.h2                     1          1
# mean_v.true.h2           1          1
# mean_v.false.h2          1          1
# sd_v.true.h2             1          1
# t0.h2                    1          1
# 
# Multivariate psrf
# 
# 1

hest.good <- summary.dmc(hsamples1.good,hyper=TRUE)

# Mean estimates not to bad, but mean_v false underestimated 
round(hest.good$statistics[,"Mean"][1:6],3)
#            A.h1            B.h1  mean_v.true.h1 mean_v.false.h1    sd_v.true.h1           t0.h1 
#           0.406           0.644           1.067           0.149           0.529           0.299 
round(hest.good$statistics[,"Mean"][1:6]-apply(ps,2,mean),3)
#            A.h1            B.h1  mean_v.true.h1 mean_v.false.h1    sd_v.true.h1           t0.h1 
#           0.008           0.032           0.042           0.102           0.019          -0.006 

# Scale OK 
round(hest.good$statistics[,"Mean"][7:12],3)
#            A.h2            B.h2  mean_v.true.h2 mean_v.false.h2    sd_v.true.h2           t0.h2 
#           0.105           0.122           0.229           0.206           0.140           0.054 
round(hest.good$statistics[,"Mean"][7:12]-apply(ps,2,sd),3)
#            A.h2            B.h2  mean_v.true.h2 mean_v.false.h2    sd_v.true.h2           t0.h2 
#           0.016           0.030           0.013           0.002           0.019           0.001 

### run with chain mixing turned off.

hsamples <- h.run.dmc(hsamples, cores=8, report=1,
                      p.migrate = .05, h.p.migrate = .05,random.phi=FALSE,random.theta=FALSE)
hsamples <- h.run.dmc(h.samples.dmc(nmc=150,samples=hsamples,add=TRUE),
                      cores=8,report=1,p.migrate = .05, h.p.migrate = .05,
                      random.phi=FALSE,random.theta=FALSE)

# Gradual increase, A scale appears to be heading to zero!
plot.dmc(hsamples,smooth=FALSE,density=FALSE,pll.chain=TRUE,hyper=TRUE,start=20)
plot.dmc(hsamples,hyper=TRUE,layout=c(2,6),smooth=FALSE,density=FALSE,start=150)

# Sample without migration
hsamples1 <- h.run.dmc(h.samples.dmc(nmc=100,samples=hsamples),
                       report=1,cores=8,random.phi=FALSE)
hsamples1 <- h.run.dmc(h.samples.dmc(nmc=400,samples=hsamples1,add=TRUE),
                       report=1,cores=8,random.phi=FALSE)

# Initial drop due to migration over-fitting. After that chains not mixing!
plot.dmc(hsamples1,smooth=FALSE,density=FALSE,pll.chain=TRUE,hyper=TRUE)
# Some A chains have zero variance (also one mean_v.false)
plot.dmc(hsamples1,hyper=TRUE,layout=c(2,6),smooth=FALSE,density=FALSE)

# Poor mixing particularly for A, B and mean_v.false scale
gelman.diag.dmc(hsamples1,hyper=TRUE)
#                 Point est. Upper C.I.
# A.h1                  1.23       1.42
# B.h1                  1.27       1.47
# mean_v.true.h1        1.13       1.22
# mean_v.false.h1       1.21       1.37
# sd_v.true.h1          1.08       1.13
# t0.h1                 1.19       1.32
# A.h2                  1.89       2.47
# B.h2                  1.11       1.20
# mean_v.true.h2        1.04       1.07
# mean_v.false.h2       1.59       1.99
# sd_v.true.h2          1.00       1.01
# t0.h2                 1.00       1.01
# 
# Multivariate psrf
# 
# 1.95

hest <- summary.dmc(hsamples1,hyper=TRUE)

# Mean estimates not to bad, but mean_v false underestimated 
round(hest$statistics[,"Mean"][1:6],3)
#            A.h1            B.h1  mean_v.true.h1 mean_v.false.h1    sd_v.true.h1           t0.h1 
#           0.394           0.652           1.066           0.148           0.529           0.297 
round(hest$statistics[,"Mean"][1:6]-apply(ps,2,mean),3)
#            A.h1            B.h1  mean_v.true.h1 mean_v.false.h1    sd_v.true.h1           t0.h1 
#          -0.005           0.040           0.041           0.101           0.019          -0.008 

# But scale well off 
round(hest$statistics[,"Mean"][7:12],3)
#            A.h2            B.h2  mean_v.true.h2 mean_v.false.h2    sd_v.true.h2           t0.h2 
#           0.056           0.133           0.232           0.173           0.139           0.054 
round(hest$statistics[,"Mean"][7:12]-apply(ps,2,sd),3)
#            A.h2            B.h2  mean_v.true.h2 mean_v.false.h2    sd_v.true.h2           t0.h2 
#          -0.034           0.041           0.016          -0.031           0.018           0.001 



### FIT RANDOM EFFECTS WITH A AS A GROUP PARAMETER

# Need a new model in which sets A to a constant in model (value doesn’t matter)
model.groupA <- model.dmc(
  p.map     = list(A = "1", B = "1", mean_v = "M", sd_v = "M", t0 = "1", st0 = "1"),
  match.map = list(M = list(s1 = 1, s2 = 2)),
  factors   = list(S = c("s1", "s2")),
  constants = c(st0 = 0, sd_v.false = 1, A = .1),  
  responses=c("r1","r2"),
  type="norm")

# Must halso drop A in data level prior
p.prior.groupA <- prior.p.dmc(
  p1=c(B=1,mean_v.true=1,mean_v.false=1,sd_v.true=1,t0=1),p2=rep(1,5),
  lower=c(0,NA,NA,0,.1),upper=c(NA,NA,NA,NA,1)
)

# Also drop A in sigma prior
sigma.prior.groupA <- prior.p.dmc(
  dists = rep("gamma", length(p.prior.groupA)),
  p1=c(B=1,mean_v.true=1,mean_v.false=1,sd_v.true=1,t0=1),p2=rep(1,5)
)
par(mfcol=c(2,3))
for (i in names(sigma.prior.groupA)) plot.prior(i,sigma.prior.groupA)

# Make new data.model and pp.prior
data.model.groupA <- data.model.dmc(raw.data, model.groupA)
pp.prior.groupA <- list(mu.prior, sigma.prior.groupA)

# Run samples with migration
hsamplesA.group <- h.samples.dmc(nmc = 50, p.prior=p.prior.groupA, data=data.model.groupA, 
                                 pp.prior=pp.prior.groupA,thin=10)
hsamplesA.group <- h.run.dmc(hsamplesA.group, cores=8, report=1,
                             p.migrate = .05, h.p.migrate = .05,random.phi=FALSE)
hsamplesA.group <- h.run.dmc(h.samples.dmc(nmc=150,samples=hsamplesA.group,add=TRUE),
                             cores=8,report=1,p.migrate = .05, h.p.migrate = .05,random.phi=FALSE)

# Much faster convergence than before.
plot.dmc(hsamplesA.group,smooth=FALSE,density=FALSE,pll.chain=TRUE,hyper=TRUE,start=20)
# Note that there is plot of scale for A as not estimated. A parameter looks less variable.
plot.dmc(hsamplesA.group,hyper=TRUE,layout=c(2,6),smooth=FALSE,density=FALSE,start=100)

# Some quite large autocorreltion in group parameter
acf.dmc(hsamplesA.group,hyper=TRUE,start=100,par=A)

# Turn off migration. Turns out need more samples as effective sample size small.
hsamplesA1.group <- h.run.dmc(h.samples.dmc(nmc=1000,samples=hsamplesA.group),
                              report=1,cores=8,random.phi=FALSE)
hsamplesA1.group <- h.run.dmc(h.samples.dmc(nmc=1000,samples=hsamplesA1.group,add=TRUE),
                              report=1,cores=8,random.phi=FALSE)


# Much better mixing for first 500 but then some chains wander off.
plot.dmc(hsamplesA1.group,smooth=FALSE,density=FALSE,pll.chain=TRUE,hyper=TRUE)

# However some A and B problems, and mean_v.false.h2 displays zero variance trap
plot.dmc(hsamplesA1.group,hyper=TRUE,layout=c(2,6),smooth=FALSE,density=FALSE)
gelman.diag.dmc(hsamplesA1.group,hyper=TRUE)
#                 Point est. Upper C.I.
# A.h1                  1.12       1.22
# B.h1                  1.11       1.20
# mean_v.true.h1        1.05       1.10
# mean_v.false.h1       1.08       1.15
# sd_v.true.h1          1.03       1.05
# t0.h1                 1.08       1.15
# B.h2                  1.03       1.06
# mean_v.true.h2        1.01       1.02
# mean_v.false.h2       1.18       1.34
# sd_v.true.h2          1.00       1.00
# t0.h2                 1.00       1.00
# 
# Multivariate psrf
# 
# 1.15

# Relatively short series due to high autocorrelation
effectiveSize.dmc(hsamplesA1.group,hyper=TRUE)
#            A.h1            B.h1  mean_v.true.h1 mean_v.false.h1 
#             311             297             478             249 
#    sd_v.true.h1           t0.h1            B.h2  mean_v.true.h2 
#             951             507            1399            4540 
# mean_v.false.h2    sd_v.true.h2           t0.h2 
#             495            6508            4469 


hestA.group <- summary.dmc(hsamplesA1.group,hyper=TRUE)
# Similar mean recovery, mean_v.false still off
round(hestA.group$statistics[,"Mean"][1:6],3)
#            A.h1            B.h1  mean_v.true.h1 mean_v.false.h1    sd_v.true.h1           t0.h1 
#           0.400           0.639           1.056           0.128           0.526           0.300 
round(hestA.group$statistics[,"Mean"][1:6]-apply(ps,2,mean),3)
#            A.h1            B.h1  mean_v.true.h1 mean_v.false.h1    sd_v.true.h1           t0.h1 
#           0.001           0.026           0.031           0.081           0.015          -0.005 

# Slightly better SD recovery, B worst.
round(hestA.group$statistics[,"Mean"][7:11],3)
#            B.h2  mean_v.true.h2 mean_v.false.h2    sd_v.true.h2           t0.h2 
#           0.140           0.232           0.170           0.139           0.054 
round(hestA.group$statistics[,"Mean"][7:11]-apply(ps,2,sd)[-1],3)
#            B.h2  mean_v.true.h2 mean_v.false.h2    sd_v.true.h2           t0.h2 
#           0.048           0.016          -0.034           0.018           0.001 


### FIT RANDOM EFFECTS WITH A and mean_v.false AS GROUP PARAMETERS

# DMC can do several group parameters at once. Here we explore whether making
# mean_v.false a group parameter to see if that fixes estimation problems, 
# turns out to make things worse. Group parameters may only be suitable for
# trial to trial variability parameters that have weak influence on data. 

# Set both A and mean_v.false to a constant in model 
model.groupAvf <- model.dmc(
  p.map     = list(A = "1", B = "1", mean_v = "M", sd_v = "M", t0 = "1", st0 = "1"),
  match.map = list(M = list(s1 = 1, s2 = 2)),
  factors   = list(S = c("s1", "s2")),
  constants = c(st0 = 0, sd_v.false = 1, A = .1, mean_v.false = 1),  
  responses=c("r1","r2"),
  type="norm")


# Data level prior drops A and mean_v.false
p.prior.groupAvf <- prior.p.dmc(
  p1=c(B=1,mean_v.true=1,sd_v.true=1,t0=1),p2=rep(1,4),
  lower=c(0,NA,0,.1),upper=c(NA,NA,NA,1)
)

# Drop A and mean_v.false in sigma prior
sigma.prior.groupAvf <- prior.p.dmc(
  dists = rep("gamma", length(p.prior.groupAvf)),
  p1=c(B=1,mean_v.true=1,sd_v.true=1,t0=1),p2=rep(1,4)
)
par(mfcol=c(2,2))
for (i in names(sigma.prior.groupAvf)) plot.prior(i,sigma.prior.groupAvf)

data.model.groupAvf <- data.model.dmc(raw.data, model.groupAvf)
pp.prior.groupAvf <- list(mu.prior, sigma.prior.groupAvf)


hsamplesAvf.group <- h.samples.dmc(nmc = 50, p.prior=p.prior.groupAvf, 
                                   data=data.model.groupAvf, pp.prior=pp.prior.groupAvf,thin=10)

# Run with migration
hsamplesAvf.group <- h.run.dmc(hsamplesAvf.group, cores=8, report=1,
                               p.migrate = .05, h.p.migrate = .05,random.phi=FALSE)
hsamplesAvf.group <- h.run.dmc(h.samples.dmc(nmc=150,samples=hsamplesAvf.group,add=TRUE),
                               cores=8,report=1,p.migrate = .05, h.p.migrate = .05,random.phi=FALSE)

# # Much faster convergence than before.
plot.dmc(hsamplesAvf.group,smooth=FALSE,density=FALSE,pll.chain=TRUE,hyper=TRUE,start=20)
plot.dmc(hsamplesAvf.group,hyper=TRUE,layout=c(2,6),smooth=FALSE,density=FALSE,start=100)

# Run with migration off
hsamplesAvf1.group <- h.run.dmc(h.samples.dmc(nmc=100,samples=hsamplesAvf.group),
                                report=1,cores=8,random.phi=FALSE)
hsamplesAvf1.group <- h.run.dmc(h.samples.dmc(nmc=400,samples=hsamplesAvf1.group,add=TRUE),
                                report=1,cores=8,random.phi=FALSE)

# Strong decrees in likelihoods then chains not mixing well, 
plot.dmc(hsamplesAvf1.group,smooth=FALSE,density=FALSE,pll.chain=TRUE,hyper=TRUE)

# Several parameters become unstable
plot.dmc(hsamplesAvf1.group,hyper=TRUE,layout=c(2,6),smooth=FALSE,density=FALSE)
gelman.diag.dmc(hsamplesAvf1.group,hyper=TRUE,start=100)
#                 Point est. Upper C.I.
# A.h1                  1.19       1.38
# B.h1                  1.15       1.30
# mean_v.true.h1        1.11       1.21
# mean_v.false.h1       1.19       1.39
# sd_v.true.h1          1.01       1.03
# t0.h1                 1.02       1.05
# B.h2                  1.11       1.22
# mean_v.true.h2        1.15       1.31
# sd_v.true.h2          1.18       1.36
# t0.h2                 1.05       1.10
# 
# Multivariate psrf
# 
# 1.25


hestAvf.group <- summary.dmc(hsamplesAvf1.group,hyper=TRUE,start=100)
# Mean recovery worse
round(hestAvf.group$statistics[,"Mean"][1:6],3)
#            A.h1            B.h1  mean_v.true.h1 mean_v.false.h1    sd_v.true.h1           t0.h1 
#           0.428           0.841          -2.193           0.210           0.380           0.223 
round(hestAvf.group$statistics[,"Mean"][1:6]-apply(ps,2,mean),3)
#            A.h1            B.h1  mean_v.true.h1 mean_v.false.h1    sd_v.true.h1           t0.h1 
#           0.029           0.229          -3.218           0.163          -0.130          -0.082 
# SD recovery, terrible
round(hestAvf.group$statistics[,"Mean"][7:10],3)
#           B.h2 mean_v.true.h2   sd_v.true.h2          t0.h2 
#          0.234         11.467          2.967          0.144 
round(hestAvf.group$statistics[,"Mean"][7:10]-apply(ps,2,sd)[-c(1,4)],3)
#           B.h2 mean_v.true.h2   sd_v.true.h2          t0.h2 
#          0.142         11.252          2.846          0.090 


### CONCLUSION: Nothing worked very well

save_data(raw.data,hsamples.good,hsamples1.good,hsamples,hsamples1,hsamplesA.group,hsamplesA1.group,
          hsamplesAvf.group,hsamplesAvf1.group,file="dmc_4_6group.RData")

