##################  DMC Lesson 2: Priors, Posteriors and Adding New Models
#
# Note: Before running these scripts, the current working directory must be set 
# to the top-level folder containing the dmc and tutorial subfolders 
# 
# NB: Lessons 2.2 uses the LNR model, either "lnr.R" or "lnrPP.R". 
#       The former uses the natural parameterization (see also lesson 1.4 and 
#       2.1), the latter a postive parameterization ("PP"), estimating the zero 
#       bounded parameters (t0 and sdlog) on the log scale (see transform.dmc 
#       in lnrPP, which exponentiates these parameters before they are passed 
#       to n1PDF.lnr)

###Lesson 2.2:  Advanced Priors

rm(list=ls())
source ("dmc/dmc.R")

### ---------------------------We start with the natural parametrization

load_model ("LNR","lnr.R")

# Define a simple design with one stimulus factor and a 5 parameter model
# Only a M(atch) effect, just like in Lesson 1.4
# NB: Remember LNR does not have a scaling parameter, so only st0 is fixed in 
#     constants (see also Lesson 2.1).
model <- model.dmc(type="lnr",constants=c(st0=0),
                   p.map=list(meanlog="M",sdlog="M",t0="1",st0="1"),
                   match.map=list(M=list(left="LEFT",right="RIGHT")),
                   factors=list(S=c("left","right")),
                   responses=c("LEFT","RIGHT"))

# The following parameter setting provides a reasonable RT distribution, 
# with approximately 25% errors:
p.vector  <- c(meanlog.true=-1,meanlog.false=0,sdlog.true=1,sdlog.false=1,t0=.2)

# This generates data: 
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1000),model)
head(data.model)

# This is how the data look like:
layout(1)
plot.score.dmc(data.model)

### Example 1: Lognormal priors on the natural scale

# Remember, meanlog.true and meanlog.false are unbounded;sdlog.true, 
# sdlog.false, and t0 must be positive
p.prior <- prior.p.dmc(
  # tnorm for unbounded meanlog
  # lnorm for sdlog and t0 (sdlog>0 and t0>0)
  dists = c("tnorm","tnorm","lnorm","lnorm","lnorm"),
  p1=c(meanlog.true=-1,meanlog.false=0,     # Sets tnorm mean = true value; prior mode = true value     
       sdlog.true=1,sdlog.false=1,          # Sets lnorm meanlog = 1
       t0=log(.2)+1),                       # Sets log(t0) = log(.2)+1 = -0.61           
  p2=c(meanlog.true=1,meanlog.false=1,      # Sets tnorm SD = 1
       sdlog.true=1,sdlog.false=1,          # Sets lnorm sdlog = 1; prior mode 
       #  = exp(meanlog-sdlog^2) 
       #  = exp(1-1^2) = 1 = true value
       t0=1)                                # Sets lnorm sdlog = 1; prior mode 
  # = exp(meanlog-sdlog^2)  
  # = exp((log(.2)+1)-1^2) 
  # = 0.2 = true value
  # No lower or upper bounds, use default values
) 

# NB. To make mode of lnrom = x: meanlog = log(x) + sdlog^2   
# That's why we set p1 = log(0.2) + 1

unlist(p.prior$meanlog.true)  # Default lower and upper are added
unlist(p.prior$sdlog.true)    # Default lower is added; no upper bound for lnorm
unlist(p.prior$t0)            # Default upper is added; no upper bound for lnorm

# Plot p.prior: 
par(mfcol=c(2,3))
for (i in names(p.prior)) plot.prior(i,p.prior,ylim=c(0,1.2))

### Example 2: Shifted gamma and shifted lognormal priors on the natural scale

# Use shifted gamma and shifted lognormal to make sdlog<.05 and t0<.1 impossible
p.prior2 <- prior.p.dmc(
  # tnorm for unbounded meanlog
  # gamma for sdlog (sdlog>0)
  # lnorm for t0 (t0>0)
  dists = c("tnorm","tnorm","gamma","gamma","lnorm"),
  p1=c(meanlog.true=-1,meanlog.false=0,      # Sets tnorm mean = true value; 
       #  prior mode = true value     
       sdlog.true=2,sdlog.false=2,           # Sets gamma shape = 2  
       t0=log(.1)+1),                        # Sets lnorm meanlog = log(.1)+1 = -1.3     
  p2=c(meanlog.true=1,meanlog.false=1,       # Sets tnorm SD = 1
       sdlog.true=.95,sdlog.false=.95,       # Sets gamma scale = .95;
       # gamma mode = (shape-1)*scale
       t0=1),                                # Sets lnorm sdlog = 1
  lower=c(meanlog.true=NA,meanlog.false=NA,  
          sdlog.true=.05,sdlog.false=.05,    # Shifts with .05; shifted gamma mode 
          # = lower + (shape-1)*scale 
          # = .05 + (2-1)*.95 = 1 = true value
          t0=.1)                             # Shift with .1;  
  # shifted lnorm mode 
  # = lower + exp(meanlog-sdlog^2) 
  # = .1 + exp((log(.1)+1)-1^2) 
  # = 0.2 = true value
)

p.prior2$sdlog.true   # lower is set to .05; no upper bound for gamma
p.prior2$t0           # lower is set to .1; no upper bound for lnorm

# Plot p.prior2: 
par(mfcol=c(2,3))
for (i in names(p.prior2)) plot.prior(i,p.prior2,ylim=c(0,1.2))

# Zoom in:
par(mfcol=c(1,2))
plot.prior("sdlog.true",p.prior2,ylim=c(0,0.4),xlim=c(0,1))
plot.prior("t0",p.prior2,ylim=c(0,0.6),xlim=c(0,1)) 

# Compare lnorm from p.prior and shifted gamma from p.prior2:
par(mfrow=c(1,2))
# lnorm with mode=1; lower bound=1
plot.prior("sdlog.true",p.prior,n.point=1e4,xlim=c(0,5))  
# shifted gamma with mode=1; lower bound=.5, steeper leading edge and skinnier tail
plot.prior("sdlog.true",p.prior2,n.point=1e4,xlim=c(0,5)) 

# Compare lnorm from p.prior and shifted lnorm from p.prior2:
par(mfrow=c(1,2))
# lnorm with mode=.2; lower bound=1
plot.prior("t0",p.prior,n.point=1e4,xlim=c(0,4))  
# shifted lnorm with mode=.2; lower bound=.1
# variance is smaller (meanlog smaller), skew is same (sdlog same)
plot.prior("t0",p.prior2,n.point=1e4,xlim=c(0,4)) 


### --------------------------- We continue with the positive parametrization

# Here we estimate the zero bounded parameters (t0 and sdlog) on the log scale 
# (see transform.dmc in lnrPP, which exponentiates these parameters 
# before they are passed to n1PDF.lnr)

load_model ("LNR","lnrPP.R")

p.vector
# meanlog.true meanlog.false    sdlog.true   sdlog.false    t0 
# -1.0           0.0           1.0           1.0           0.2
logp.vector = c(p.vector[1:2],log(p.vector[3:5]))
logp.vector
# meanlog.true meanlog.false    sdlog.true   sdlog.false    t0 
# -1.000000      0.000000      0.000000      0.000000     -1.609438 

### Example 3: Normal priors on the log scale
# OK to use normals on sdlog and t0, because they are now on the log scale 
# and not bounded by 0:
p.prior3 <- prior.p.dmc(
  # tnorm is the default prior distribution if a vector of "dists" isn't specified
  p1=logp.vector,        # sets tnorm mode = true value for all parameters    
  p2=c(1,1,1,1,1),       # Sets tnorm SD = 1 for all parameters
  lower=rep(0,5),        # Sets lower = 0 for all parameters
  upper=rep(1,5))        # Sets upper = 1 for all parameters

p.prior3$meanlog.true
p.prior3$sdlog.true
p.prior3$t0

# Plot p.prior3: 
par(mfcol=c(2,3))
for (i in names(p.prior3)) plot.prior(i,p.prior3,xlim=c(-1,2),ylim=c(0,2))

# This prior makes absolutely no sense:
logp.vector
# meanlog.true meanlog.false    sdlog.true   sdlog.false      t0 
#  -1.000000      0.000000      0.000000      0.000000     -1.609438

# You are free to chose your priors so that they reflect your (subjective) beliefs, 
# but they shouldn't exclude probable parameter values (here the true values).

# Let's try again:
p.prior4 <- prior.p.dmc(
  # tnorm is the default prior distribution if a vector of "dists" isn't specified
  p1=logp.vector,       # sets tnorm mode = true value for all parameters    
  p2=c(1,1,1,1,1)       # Sets tnorm SD = 1 for all parameters     
  # No crazy bounds enforced
)

p.prior4$meanlog.true
p.prior4$sdlog.true
p.prior4$t0

# Plot p.prior4: 
par(mfcol=c(2,3))
for (i in names(p.prior4)) plot.prior(i,p.prior4,xlim=c(-5,5),ylim=c(0,0.5))

# sdlog and t0 are transformed parameters, on a log scale. 
# To see what the prior looks like on the natural scale, set the "trans argument 
# to the name of the inverse (exp) transformation:
par(mfrow=c(1,4))
plot.prior("sdlog.true",p.prior4)             # log scale
plot.prior("sdlog.true",p.prior4,trans="exp") # natural scale
plot.prior("t0",p.prior4)                     # log scale
plot.prior("t0",p.prior4,trans="exp")         # natural scale

# Compare p.prior (lognormal prior on natural scale) and p.prior4 (normal prior 
# on log scale). Gives the same shape (check Wikipedia to see how to tinker with 
# the parameters to match them exactly):
par(mfrow=c(1,4))
plot.prior("t0",p.prior)                      # lognormal prior on natural scale
plot.prior("t0",p.prior,trans="log")          # ...and then log transformed
plot.prior("t0",p.prior4)                     # normal prior on log scale with 
plot.prior("t0",p.prior4,trans="exp")         # ...and then exponentiated

# NB: As this example illustrates, because of the prior, scale matters in Bayes! 
#     Make sure that your priors make sense!

# You can also indicate the appropriate transformation to the natural scale when
# specifying the prior. You can do this for all parameters or only those 
# that are not on the natural scale (i.e., an "identity" transform) using names.
p.prior5 <- prior.p.dmc(
  p1=logp.vector,                     
  p2=c(1,1,1,1,1), 
  # trasnform sdlog and t0 back back to natural scale
  untrans=c(sdlog.true="exp",sdlog.false="exp",t0="exp")  
)

p.prior5$meanlog.true
# attr(,"untrans")
# [1] "identity"

p.prior5$sdlog.true
# attr(,"untrans")
# sdlog.true 
# "exp"

p.prior5$t0
# attr(,"untrans")
# t0 
# "exp"

# By default plotting will apply these transformation:
par(mfcol=c(2,3)); for (i in names(p.prior5)) plot.prior(i,p.prior5)

# To turn off the automatic transformation 
par(mfcol=c(2,3)); for (i in names(p.prior5)) plot.prior(i,p.prior5,natural=FALSE)


### --------------------------------------- Choose your priors wisely! 

###--- Prior bounds (see also p.prior3)

# "log.prior.dmc" calculates the log of the value of the prior density of a set 
# of parameter values 

# We use p.prior2 as an example: 
# p.prior2 specifies tnorm for meanlogs, shifted (.05) gamma for sdlogs and 
# shifted (.1) lnorm for t0:
par(mfcol=c(2,3))
for (i in names(p.prior2)) plot.prior(i,p.prior2,ylim=c(0,1.2))

# As our setup ensures that the modes of the prior densities are at the 
# p.vector values, the maximum density values are:
exp(log.prior.dmc(p.vector,p.prior2))
# meanlog.true meanlog.false    sdlog.true   sdlog.false      t0 
# 0.3989423     0.3989423     0.3872415     0.3872415     2.4197072 

# Moving below these values decreases the density:
exp(log.prior.dmc(p.vector/.5,p.prior2))
# meanlog.true meanlog.false    sdlog.true   sdlog.false      t0 
# 0.2419707     0.3989423     0.2774220     0.2774220     1.3233575 

# The same is true for moving above:
exp(log.prior.dmc(p.vector*1.5,p.prior2))
# meanlog.true meanlog.false    sdlog.true   sdlog.false      t0 
# 0.3520653     0.3989423     0.3491807     0.3491807     1.9029780 

# Given the lower bounds on sdlog and t0,
# sufficiently small values have zero density:
exp(log.prior.dmc(p.vector-c(0,0,.99,.99,.15),p.prior2))
#  meanlog.true meanlog.false    sdlog.true   sdlog.false       t0 
#   0.3989423     0.3989423     0.0000000     0.0000000     0.0000000 

# and hence log(density) of -Inf
log.prior.dmc(p.vector-c(0,0,.99,.99,.15),p.prior2)
#  meanlog.true meanlog.false    sdlog.true   sdlog.false     t0 
#    -0.9189385    -0.9189385          -Inf          -Inf    -Inf 

# As posterior likelihood = data likelihood * prior,
# this implies that there is no chance of accepting such values,
# so posterior likelihood is also zero.

###--- Finally, let's revisit prior sensitivity

# posterior likelihood = data likelihood * prior (not the same as posterior 
#                                         distribution, needs to be normalized!)
# log posterior likelihood = log data likelihood + log prior
# "log.posterior.dmc" calculates the summed log posterior likelihood, which is
#  the sum of the log likelihood and the sum of the log prior given specific 
#  values of the parameters:
log.posterior.dmc(p.vector,p.prior2,data.model) 
#  -43128.77  (this changes for each new data set)

# The value of the summed log likelihood is: 
sum(log(likelihood.dmc(p.vector,data.model)))
#  -43125.92  (this chages for each new data set)

# The value of the summed log prior is:
sum(log.prior.dmc(p.vector,p.prior2))
#  -2.851644   (this does not depend on the data)

# In this case the prior is overwhelmed by the data likelihood due to the very
# large number of data points (1,000*2 RTs), meaning that the
# data likelihood contributes more to log posterior likelihood than the prior,
# except where parameters cross bounds or are in the far tail of the prior: 

log.posterior.dmc(p.vector*1e4,p.prior2,data.model) 
# -50057136    
sum(log(likelihood.dmc(p.vector*1e4,data.model)))
# -46051.7     
sum(log.prior.dmc(p.vector*1e4,p.prior))
# -49990132    

# And this is what happens with scarce data
# Prior is more influential even in the high likelihood region:
data.model2 <- data.model.dmc(simulate.dmc(p.vector,model,n=5),model)

log.posterior.dmc(p.vector,p.prior2,data.model2) 
# -9.041214
sum(log(likelihood.dmc(p.vector,data.model2)))
# -6.18957
sum(log.prior.dmc(p.vector,p.prior2))
# -2.851644
