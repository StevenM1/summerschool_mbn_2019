##################  DMC Lesson 6: More Models

### Lesson 6.3: Go-NoGo Models

rm(list=ls()) 
# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")

# load_data ("dmc_6_3.RData")

# LNR version, two choice, same model as Lessons 3_1

load_model ("LNR-GoNoGo","lnrgng.R")

# Setup is like a standard LNR except for type="lnrgng"
model <- model.dmc(p.map=list(meanlog="M",sdlog="M",t0="1",st0="1"),
                   match.map=list(M=list(s1=1,s2=2)),constants=c(st0=0),
                   factors=list(S=c("s1","s2")),responses=c("r1","r2"),
                   type="lnrgng")
p.vector  <- c(meanlog.true=-1,meanlog.false=0,
               sdlog.true=1,sdlog.false=1,t0=.2)

# Note that the FIRST accumulator is assumed to be the stop accumulator. Data 
# looks much as usual except now NA for RT where the first (no-go) "response"
# occurs. Note you can both fail to no-go for no-go stimulus and fail to go
# for go stimulus.
simulate.dmc(p.vector,model,n=5)
#     S  R        RT
# 1  s1 r1        NA
# 2  s1 r2 0.6634085
# 3  s1 r1        NA
# 4  s1 r1        NA
# 5  s1 r2 1.9323661
# 6  s2 r2 1.2065659
# 7  s2 r2 0.6394280
# 8  s2 r1        NA
# 9  s2 r2 0.5671596
# 10 s2 r1        NA


# Lets get a large sample
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)

# Prior as usual
p.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","gamma","gamma","beta"),
  p1=c(meanlog.true=-1,meanlog.false=0, sdlog.true=2,sdlog.false=2,t0=1),                           
  p2=c(1,1,1,1,1),lower=c(NA,NA,NA,NA,.1),upper=c(NA,NA,NA,NA,1)
)
par(mfcol=c(2,3)); for (i in names(p.prior)) plot.prior(i,p.prior)

# Do sampling
samples <- samples.dmc(nmc=400,p.prior,data.model)
samples <- run.dmc(samples, report = 10, cores=15,p.migrate=.05)

# Converged by around 250
plot.dmc(samples,pll.chain=TRUE,start=250)

samples1 <- run.dmc(samples.dmc(samples=samples,nmc=200),cores=15)
samples1 <- run.dmc(samples.dmc(samples=samples1,add=TRUE,nmc=100),cores=15)

# Looks good
plot.dmc(samples1,pll.chain=TRUE)
gelman.diag.dmc(samples1)
#               Point est. Upper C.I.
# meanlog.true        1.02       1.03
# meanlog.false       1.02       1.04
# sdlog.true          1.02       1.03
# sdlog.false         1.02       1.03
# t0                  1.01       1.02
# 
# Multivariate psrf
# 
# 1.03

# OK Size
effectiveSize.dmc(samples1)
#  meanlog.true meanlog.false    sdlog.true   sdlog.false            t0 
#           322           291           272           296           304 
# Good recovery
check.recovery.dmc(samples1,p.vector)          
#               meanlog.true meanlog.false sdlog.true sdlog.false   t0
# True                  -1.00          0.00       1.00        1.00 0.20
# 2.5% Estimate         -1.01         -0.03       0.97        0.97 0.19
# 50% Estimate          -0.99          0.01       0.99        1.00 0.20
# 97.5% Estimate        -0.97          0.04       1.02        1.03 0.20
# Median-True            0.01          0.01      -0.01        0.00 0.00

# Excellent fit
pp <- post.predict.dmc(samples1)
plot.pp.dmc(pp,style = "cdf")


### LBA version, two choice, same as lesson 3_3

load_model ("LBA-GoNoGo","lba_Bgng.R")

# Same simple design as for previous LNR examples
model <- model.dmc(p.map=list(A="1",B="1",mean_v="M",sd_v="1",t0="1",st0="1"),
                   constants=c(st0=0,sd_v=1),match.map=list(M=list(s1=1,s2=2)),
                   factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="normgng")


# Simulate some data, with around 65% accuracy
p.vector  <- c(A=.25,B=.35,mean_v.true=1,mean_v.false=.25,t0=.2)
data.model1 <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)

p.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm","beta"),
  p1=c(A=.3,B=.3,mean_v.true=1,mean_v.false=0,t0=1),                           
  p2=c(1,1,3,3,1),lower=c(0,0,NA,NA,.1),upper=c(NA,NA,NA,NA,1)
)
par(mfcol=c(2,3)); for (i in names(p.prior)) plot.prior(i,p.prior)


# Do sampling
samples2 <- samples.dmc(nmc=400,p.prior,data.model1)
samples2 <- run.dmc(samples2, report = 10, cores=15,p.migrate=.05)

# Converged by around 200
plot.dmc(samples2,pll.chain=TRUE,start=200)

samples3 <- run.dmc(samples.dmc(samples=samples2,nmc=200),cores=15)
samples3 <- run.dmc(samples.dmc(samples=samples3,add=TRUE,nmc=200),cores=15)

# Looks good
plot.dmc(samples3,pll.chain=TRUE)
gelman.diag.dmc(samples3)

# OK Size
effectiveSize.dmc(samples3)
#            A            B  mean_v.true mean_v.false           t0 
#          174          165          199          194          171 

# Good recovery
check.recovery.dmc(samples3,p.vector)          
#                   A     B mean_v.true mean_v.false   t0
# True           0.25  0.35        1.00         0.25 0.20
# 2.5% Estimate  0.07  0.25        0.81         0.02 0.17
# 50% Estimate   0.25  0.34        0.97         0.23 0.20
# 97.5% Estimate 0.35  0.50        1.14         0.46 0.23
# Median-True    0.00 -0.01       -0.03        -0.02 0.00


# Excellent fit
pp1 <- post.predict.dmc(samples3)
plot.pp.dmc(pp1,style = "cdf")



# save_data (data.model,samples,samples1,pp,
#            data.model1,samples2,samples3,pp1,file="dmc_6_3.RData")


