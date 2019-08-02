rm(list=ls())
source ("dmc/dmc.R")
load_model ("EXG-SS","exgSSprobitN.R")

is.tf <- TRUE
is.gf <- TRUE
use.staircase <- TRUE

  model <- model.dmc(
    # SS stands for trial type (GO or Stop-signal [SS]):
    factors=list(S=c("s1","s2","s3"),SS=c("GO","SS"),NC=c("a2","a3")),
    # NR stands for "No response", i.e., go omission & succesfull inhibitions:
    responses=c("NR","r1","r2","r3"),
    constants=c(N.a2=3,N.a3=4),
    # Match scores correct responses for each GO stimulus as usual, and also scores 
    # the "correct" stimulus corresponding to an NR response, but the latter has
    # no effect (it just avoids a standard check making sure that each response is scored): 
    match.map=list(M=list(s1="r1",s2="r2",s3="r3",s1="NR")),
    p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1",N="NC"),
    type="exgss")

# Parameter vector names are: ( see attr(,"p.vector") )
#  [1] "mu.true"     "mu.false"    "sigma.true"  "sigma.false" "tau.true"   
#  [6] "tau.false"   "tf"          "muS"         "sigmaS"      "tauS"       
# [11] "gf"         
# 
# Constants are (see attr(,"constants") ):
# N.a2 N.a3 
#    3    4 
# 
# Model type = exgss

#     p.vector  <- c(mu.true=.5,mu.false=.60,muS=.2,
#                  sigma.true=.05,sigma.false=.030,sigmaS=.03,
#                  tau.true=.08,tau.false=.04,tauS=.05,
#                  tf=-6,gf=-6)

    # This parameter setting will produce about 25% (relatively slow) errors:
    p.vector  <- c(  mu.true=.5,    mu.false=.60,    muS=.2,
                 sigma.true=.05,sigma.false=.030,sigmaS=.03,
                   tau.true=.08,   tau.false=.04,  tauS=.05,
                 tf=qnorm(.1),gf=qnorm(.1))

    # check.p.vector(p.vector,model)

# # This's how (toy) stop-signal data look like with fixed and staircase SSDs
# # Full factorial 12 cells, but remove s3 in a2 for both GO and SS
# n <- c(6,6,0,6,6,0,6,6,6,6,6,6)  
# SSD <- c(rep(Inf,3),rep(.32,3),rep(Inf,3),rep(.32,3)) # Per cell, two ignored
# # SSD <- c(rep(Inf,2*6),rep(.32,2*6),rep(Inf,3*6),rep(.32,3*6)) # Per value
# 
# data_fixed <- data.model.dmc(simulate.dmc(p.vector,model,n=n,SSD=SSD),model)
# data_fixed
# likelihood.dmc(p.vector,data_fixed)
# 
# data_stair <- data.model.dmc(simulate.dmc(p.vector,model,n=n,staircase=.05, SSD=SSD),model)
# data_stair
# likelihood.dmc(p.vector,data_stair)

# Make more realistic data:
n <- c(375,375,0,125,125,0,375,375,375,125,125,125)*5
SSD <- c(rep(Inf,3),rep(.25,3),rep(Inf,3),rep(.25,3)) # Per cell, two ignored
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,n=n,SSD=SSD),model)


# Check the profiles
    p.vector  <- c(  mu.true=.5,    mu.false=.60,    muS=.2,
                 sigma.true=.05,sigma.false=.030,sigmaS=.03,
                   tau.true=.08,   tau.false=.04,  tauS=.05,
                 tf=qnorm(.1),gf=qnorm(.1))

par(mfrow=c(2,6))
profile.dmc("mu.true",.25,.75,p.vector,data)
profile.dmc("mu.false",.25,.75,p.vector,data)
profile.dmc("muS",.1,.3,p.vector,data)
profile.dmc("sigma.true",.025,.075,p.vector,data)
profile.dmc("sigma.false",.02,.04,p.vector,data)
profile.dmc("sigmaS",.02,.04,p.vector,data)
profile.dmc("tau.true",.06,.1,p.vector,data)
profile.dmc("tau.false",.02,.06,p.vector,data)
profile.dmc("tauS",.03,.07,p.vector,data)
profile.dmc("gf",qnorm(.05),qnorm(.15),p.vector,data)
profile.dmc("tf",qnorm(.05),qnorm(.15),p.vector,data)



# Plot go RT distribution;
# P(NA) = Probability of go omission;
# Accuracy is computed as correct go/(all go - go omissions):
correct <- as.numeric(data$S)==(as.numeric(data$R)-1)
layout(1)
plot.cell.density(data[data$SS=="GO",],
                  C=correct[data$SS=="GO"],
                  xlim=c(0,5),ymax=5,main="Go RTs")

# Overall accuracy (i.e., all trials-(errors + omission)
tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]
#  GO 
# 0.653  

# Show the different SSDs:
sort(tapply(as.character(data$SSD),data[,c("SS")],unique)$SS)
# "0.15" "0.2"  "0.25" "0.3"  "0.35" "0.4"  "0.45"

# Show the number of trials for each SSD:
Ns = tapply(data$RT,data$SSD,length)
Ns
# 0.15  0.2 0.25  0.3 0.35  0.4 0.45  Inf 
#    1   10   43   79   77   35    5  750 

# Show response rate:
tapply(!is.na(data$RT),data[,c("SS")],mean)
#   GO    SS 
# 0.896 0.500 

# Response rate broken down by SSD & corresponding inhibition function:
tapply(!is.na(data$RT),data[,c("SS","SSD")],mean)
layout(1)
plot_SS_if.dmc(data)  # P(Respond) increases as a function of SSD, as it should
Ns

# Median signal-respond RT per SSD:
tapply(data$RT,data$SSD,median,na.rm=TRUE)
# 0.15  0.2      0.25       0.3      0.35       0.4      0.45        Inf 
# NA  0.4169316 0.5367189 0.5176544 0.5497371 0.5468397 0.5735215 0.5588533

# Plot median signal-respond RT per SSD:
plot_SS_srrt.dmc(data) # Median SRRT increases as a function of SSD, as it should

# Show number of signal-respond RT per SSD:
Nr = tapply(!is.na(data$RT),data[,c("SS","SSD")],sum)[2,]
Nr
# 0.15  0.2 0.25  0.3 0.35  0.4 0.45  Inf 
#   0    1   10   32   47   30    5   NA 

# Signal-respond RTs should be faster than go RTs:
hist(data$RT[data$SS=="GO"],breaks="fd",main="SRRT vs. Go RT",freq=F,ylim=c(0,8))
lines(density(data$RT[data$SS=="SS"],na.rm=T),col="red",lwd=2)

# Uniform (scaled beta) priors:
p1 <- p.vector; p1[1:length(p1)] <- 1; p1[10:11] <- 3
p.prior <- prior.p.dmc(
  dists = rep("tnorm",length(p1)),p1=p.vector,p2=p1, 
  lower=c(rep(0,length(p1)-2),-6,-6),upper=c(rep(2,3),rep(.5,6),rep(6,2)) 
)
par(mfcol=c(2,6)); for (i in names(p.prior)) plot.prior(i,p.prior)

# Start sampling
samples <- samples.dmc(nmc=100,p.prior,data)
samples <- run.unstuck.dmc(samples,report = 10,cores=33,p.migrate=.05,verbose=TRUE)
layout(1)
plot.dmc(samples,pll.chain=TRUE)
save.image()

samples2 <- run.converge.dmc(samples.dmc(samples=samples,nmc=50,thin=5),
  report=10, cores=33,cut=1.1,verbose=TRUE,nmc=50,minN=500,max.try=20)
save.image()

# Posterior loglikelihood looks nice:
layout(1)
plot.dmc(samples2,pll.chain=TRUE,density=FALSE,smooth=FALSE)
# Parameter chains look like fat hairy catapillars:
plot.dmc(samples2,layout=c(2,6))

# Rhat looks good:
gelman.diag.dmc(samples2)
#             Point est. Upper C.I.
# mu.true           1.01       1.01
# mu.false          1.01       1.01
# sigma.true        1.01       1.01
# sigma.false       1.01       1.01
# tau.true          1.01       1.01
# tau.false         1.01       1.01
# tf                1.10       1.14
# muS               1.03       1.04
# sigmaS            1.03       1.04
# tauS              1.03       1.04
# gf                1.01       1.02
# 
# Multivariate psrf
# 
# 1.1

# Good parameter recovery except tf
check.recovery.dmc(samples2,p.vector)
#                mu.true mu.false sigma.true sigma.false tau.true
# True               0.5      0.6       0.05        0.03     0.08
# 2.5% Estimate      0.5      0.6       0.05        0.03     0.08
# 50% Estimate       0.5      0.6       0.05        0.03     0.08
# 97.5% Estimate     0.5      0.6       0.05        0.03     0.09
# Median-True        0.0      0.0       0.00        0.00     0.00
#                tau.false    tf  muS sigmaS tauS    gf
# True                0.04 -1.28 0.20   0.03 0.05 -1.28
# 2.5% Estimate       0.03 -4.70 0.18   0.01 0.00 -1.35
# 50% Estimate        0.04 -1.41 0.20   0.04 0.05 -1.31
# 97.5% Estimate      0.04 -1.08 0.24   0.06 0.09 -1.28
# Median-True         0.00 -0.12 0.00   0.01 0.00 -0.03

# Note as tf is on the probit scale the miss is not so bad 
round(pnorm(c(-4.7,-1.41,-1.08)),4)
# [1] 0.0000 0.0793 0.1401

# Fits 
pp <- post.predict.dmc(samples2,n.post=200,save.simulation=TRUE)

layout(1)
plot_SS_if.dmc(data=samples2$data,sim=pp)
#    SSD   n     p
# 1 0.10   2 1.000
# 2 0.15  29 0.860
# 3 0.20 186 0.715
# 4 0.25 600 0.555
# 5 0.30 986 0.425
# 6 0.35 851 0.385
# 7 0.40 379 0.640
# 8 0.45  85 0.215
# 9 0.50   7 0.495

layout(1)
plot_SS_srrt.dmc(data=samples2$data,sim=pp)
#    SSD nrt     p n.sim
# 1 0.15   2 0.922   192
# 2 0.20  27 0.505   200
# 3 0.25 159 0.970   200
# 4 0.30 436 0.275   200
# 5 0.35 548 0.340   200
# 6 0.40 301 0.925   200
# 7 0.45  78 0.515   200
# 8 0.50   7 0.950   200

save_data(data,samples,samples2,pp,file="dmc_6_5_EXG34probit.RData")


