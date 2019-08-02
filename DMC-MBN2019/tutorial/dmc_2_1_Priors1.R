##################  DMC Lesson 2: Priors, Posteriors and Adding New Models
#
# Note: Before running these scripts, the current working directory must be set 
# to the top-level folder containing the dmc and tutorial subfolders 
#
### Lesson 2.1:  Basic Priors and Posteriors

# For our Bayesian analysis of the fit of the model (LNR) to the data (as 
# simulated in our call to simulate.dmc), we firstly need to specify the prior 
# probability distribution for each of the parameters in our model (i.e., 
# meanlog.true, meanlog.false, sdlog.true, sdlog.false, t0). That is, our 
# parameters are treated as random variables having their own associated 
# probability distributions. But which distributions? Firstly, we'll specify 
# the beta distribution for all distributions - uninformative priors for which, 
# as we'll see, each value is as likely as any other within the allowable range. 
# Then we'll be more infromative specifying the normal distribution for the 
# meanlogs, gamma for the SDs, and beta only for t0.


rm(list=ls())
source ("dmc/dmc.R")

# We use the lognormal race as an example:
load_model ("LNR","lnr.R")

# Define a simple design with one stimulus factor and a 5-parameter model
# Only a M(atch) effect, just like in Lesson 1.4
# NB: Remember LNR does not have a scaling parameter, so only st0 is fixed 
#     in constants.

model <- model.dmc(type="lnr",constants=c(st0=0),
                   p.map=list(meanlog="M",sdlog="M",t0="1",st0="1"),
                   match.map=list(M=list(left="LEFT",right="RIGHT")),
                   factors=list(S=c("left","right")),
                   responses=c("LEFT","RIGHT"))

# The following parameter setting provides a reasonable RT distribution, with 
# approximately 25% errors:
p.vector  <- c(meanlog.true=-1,meanlog.false=0,sdlog.true=1,sdlog.false=1,t0=.2)

# This generates data: 
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1000),model)
head(data.model)

# This is what the data look like:
layout(1)
plot.score.dmc(data.model)

# Errors are a bit slower than correct responses:
correct = as.numeric(data.model$S)==as.numeric(data.model$R)
round(tapply(data.model$RT,list(correct),mean),2)  

###---------------------------------------------------------- Create the prior

# "prior.p.dmc" makes a list of prior distributions 
# Supported priors:
# a) Truncated normal (tnorm), where: p1 = mean, p2 = sd (relies on the msm 
#    package, see ?ptnorm)
# b) Beta, where: p1 = shape1 and p2 = shape2 (see ?pbeta, note that that the
#          uniform distribution is a special case of the beta with p1=p2=1).
# c) Gamma, where p1 = shape and p2 = scale (see ?pgamma) (NB. Watch out 
#    WinBUGS folks: p2 is not the rate!)
# d) Lognormal, where p1 = meanlog and p2 = sdlog (see ?plnorm)
#
# Use lower and upper to set the distribution support: 
# a) for tnorm these define the lower and upper bounds; NAs are filled with 
#    default values of -Inf and Inf, i.e., normal distribution (see ?pnorm)
# b) for beta these define the lower and upper bounds (i.e., scaled beta 
#    distribution); NAs are filled with default values of 0 and 1
#    --> p1=1 & p2=1 & lower=0 (default) & upper=1 (default) creates Uniform(0,1)
#    --> p1=1 & p2=1 & lower=l & upper=u creates Uniform(l,u)
# c) for gamma lower shifts the distribution to exclude small values
# d) for lognormal lower shifts the distribution to exclude small values

### Example 1: Beta (and uniform) prior

# Remember, meanlog.true and meanlog.false are unbounded; sdlog.true, 
# sdlog.false, and t0 must be positive
p.prior <- prior.p.dmc(
  # Define the beta priors
  dists = c("beta","beta","beta","beta","beta"),
  # p1=1 and p2=1 for all parameters-->Beta(1,1)
  p1=c(meanlog.true=1,meanlog.false=1,sdlog.true=1,sdlog.false=1,t0=1),                    
  p2=c(meanlog.true=1,meanlog.false=1,sdlog.true=1,sdlog.false=1,t0=1),      
  # Defines the support
  lower=c(-4,-4,0,0,0.1),     # Scales the Beta(1,1) distributions 
  upper=c(4,4,4,4,1)    
)

# Before moving on, have a look at what we just mode. p.prior is a list with one 
# entry for each element of p1 with names taken from p1:
p.prior$meanlog.true
p.prior$sdlog.true
p.prior$t0

# The output is hard to read; let's flatten the list structure and see all the 
# values in a single line by using unlist():
unlist(p.prior$meanlog.true)
# shape1 shape2  lower  upper    log 
#    1      1     -4      4      1  

# Here is how we can access the type of prior distribution:
attr(p.prior$meanlog.true,"dist")
# "beta_lu"

# Now let's just see how prior.p.dmc has used default values for our parameters
# when, instead, we just give "NA" for some of the upper/lower bounds.
# We'll see that the following gives the exact same prior specification as we 
# got before: 
p.prior2 <- prior.p.dmc(
  # Define the beta priors
  dists = c("beta","beta","beta","beta","beta"),
  # p1=1 and p2=1 for all parameters-->Beta(1,1)
  p1=c(meanlog.true=1,meanlog.false=1,sdlog.true=1,sdlog.false=1,t0=1),                    
  p2=c(meanlog.true=1,meanlog.false=1,sdlog.true=1,sdlog.false=1,t0=1),      
  # Defines the support    
  lower=c(-4,-4,NA,NA,0.1),     # Scales the Beta(1,1) distributions 
  upper=c(4,4,4,4,NA)           # NA causes default values to be used
)

# Now let's compare the p.prior and p.prior2 values for sdlog.true and t0:
unlist(p.prior$sdlog.true)
unlist(p.prior2$sdlog.true) # Default lower is added
# shape1 shape2  lower  upper    log 
#    1      1      0      4      1
unlist(p.prior$t0)
unlist(p.prior2$t0)         # Default upper is added
# shape1 shape2  lower  upper    log 
#    1.0    1.0    0.1    1.0    1.0

# Plot the prior using the convenience function "plot.prior"
# Support can be set using xlim, by default it is over bounds for beta, 
# for other distributions it uses lower and upper where defined, or
# guesses a reasonable range (but that might need tuning):
par(mfcol=c(1,2))
plot.prior("meanlog.true",p.prior2,ylim = c(0,1.2))
plot.prior("meanlog.true",p.prior2,ylim = c(0,1.2),xlim=c(0,1)) # Zooms in

# But rather than plotting the distribution for each of meanlog.true, 
# meanlog.false, etc., individually, we can plot them in a loop. 
par(mfcol=c(2,3))
for (i in names(p.prior2)) plot.prior(i,p.prior2,ylim = c(0,1.2))


### Example 2: Combination of different priors for different parameters

# We just used the beta-distribution to generate priors/as the prior 
# probability distribution for ALL of meanlog.true, meanlog.false, ... t0.
# Now we're going to mix it up, using
# (a) the truncated normal distribution to generate priors for the meanlogs,
# (b) the gamma distribution for the SDlogs, and 
# (c) the beta distribution only for t0.
# Remember: in giving arguments to prior.p.dmc() to generate values from each 
# of these distribution functions,  p1 and p2 mean different things (as needed
# by each probability function) -- specifically: mean and SD for normal, 
# shape and scale for gamma, and shape1 and shape2 for beta

p.prior3 <- prior.p.dmc(
  # tnorm for unbounded meanlog
  # gamma for sdlog (sdlog>0)
  # beta for t0 
  dists = c("tnorm","tnorm","gamma","gamma","beta"),
  p1=c(meanlog.true=-1,meanlog.false=0, # sets tnorm mean = true value; prior 
       #  mode = true value
       sdlog.true=2,sdlog.false=2,      # Sets gamma shape = 2  
       t0=1),                           
  p2=c(meanlog.true=1,meanlog.false=1,  # Sets tnorm SD = 1
       sdlog.true=1,sdlog.false=1,      # Sets gamma scale = 1; prior mode = 
       #  (shape-1)*scale = (2-1)*1=1=true value
       t0=1),                           # Creates Beta(1,1); lower & upper 
  #  creates Uniform(.1,1)
  lower=c(NA,NA,NA,NA,.1),              # NA causes default values to be used
  upper=c(NA,NA,NA,NA,NA)               # upper is ignored by gamma
)

unlist(p.prior3$meanlog.true) # Default lower and upper are added
unlist(p.prior3$sdlog.true)   # Default lower is added; upper is ignored
unlist(p.prior3$t0)           # Default upper is added

# Plot p.prior3: 
par(mfcol=c(2,3))
for (i in names(p.prior3)) plot.prior(i,p.prior3)

###------------------------------------- Does the prior matter? That depends...

# Compare two prior settings for the same data

# (1) p.prior3 as defined above; these priors are weakly informative 
# and with the exception of t0, peak at true value:
par(mfcol=c(2,3))
for (i in names(p.prior3)) plot.prior(i,p.prior3,ylim = c(0,1.2))

# (2) "Misspecified" diffuse prior:
p.prior4 <- prior.p.dmc(
  # tnorm for unbounded meanlog parameters
  # gamma sdlog parameters (sdlog>0)
  # beta for t0
  dists = c("tnorm","tnorm","gamma","gamma","beta"),
  p1=c(meanlog.true=0,meanlog.false=0,  # sets prior mode to 0
       sdlog.true=4,sdlog.false=4,      # Sets gamma shape = 4  
       t0=1),                           
  p2=c(meanlog.true=15,meanlog.false=15, # Sets tnorm SD = 15; very diffuse tnorm prior
       sdlog.true=1,sdlog.false=1,       # Sets gamma scale to 1; prior mode = 
       #  (shape-1)*scale = (4-1)*1=3!=true value
       t0=1),                            
  lower=c(NA,NA,NA,NA,NA),   
  upper=c(NA,NA,NA,NA,NA)                
)

par(mfcol=c(2,3))
for (i in names(p.prior4)) plot.prior(i,p.prior4,ylim = c(0,1.2))

# Now we're going to shift from just generating priors to seeing the 
# posterior distributions that come from combining the priors with their 
# likelihood. You will look at how to estimates posterior distribution in a
# later lesson. For now you will use dmc_2.1.RData, which contains the results  
# fitting the LNR model to two data sets, data.big features 2x1000 RTs and 
# of data.small contains 2x15 RTs. Both data sets were fit using a weakly 
# informative prior (p.prior3) and a diffuse prior (p.prior4):
load_data("dmc_2.1.RData")

### Prior doesn't really matter with sufficient data (as long as it's not 
#    crazy...see lesson 2.2)

data.big <- data.model.dmc(simulate.dmc(p.vector,model,n=1000),model)
head(data.big)

# Compare results obtained with p.prior3 (samples.big3) and p.prior4 
# (samples.big4) for data.big

# Priors and posteriors for data.big with p.prior3
# Red is prior, black is posterior
# Posteriors are nicely updated: 
plot.dmc(samples.big3,layout=c(2,3),p.prior=p.prior3)
p.vector
# meanlog.true meanlog.false    sdlog.true   sdlog.false   t0 
# -1.0           0.0           1.0           1.0           0.2 

# Note that you can control the attributes of the prior and posterior lines 
# seperately. We can also turn off the "rug" plot of individual posterior
# observations shown as tick marks on the x-axis by setting the line width
# to zero.
plot.dmc(samples.big3,layout=c(2,3),p.prior=p.prior3,
  prior.lty=2,prior.col="green",prior.lwd=2,
  post.lty=3,post.col="red",post.lwd=3,rug.lwd=0)

# We're now going to look at the credible intervals for our estimates to 
# see how well the model did. We can use two methods to retrieve these 
# estimates: summary.dmc() or check.recovery.dmc().

# True values are well within 95% CI:
summary.dmc(samples.big3)
# ...
# 2. Quantiles for each variable:
#  
#                 2.5%      25%       50%     75%    97.5%
# meanlog.true  -1.06824 -1.03440 -1.016088 -0.9961 -0.95915
# meanlog.false -0.08316 -0.03564 -0.008615  0.0188  0.06787
# sdlog.true     0.95699  0.98806  1.004447  1.0220  1.05399
# sdlog.false    0.94723  0.98585  1.006268  1.0295  1.07159
# t0             0.19469  0.19947  0.201747  0.2037  0.20718

# Parameter recovery can be summarized as follows:
check.recovery.dmc(samples.big3,p.vector)
#                   meanlog.true meanlog.false sdlog.true sdlog.false   t0
# True                  -1.00          0.00       1.00        1.00     0.20
# 2.5% Estimate         -1.07         -0.08       0.96        0.95     0.19
# 50% Estimate          -1.02         -0.01       1.00        1.01     0.20
# 97.5% Estimate        -0.96          0.07       1.05        1.07     0.21
# Median-True           -0.02         -0.01       0.00        0.01     0.00

# So, the true values are well within the 95% CI. Specifically, for example, 
# having set meanlog.true to equal -1, and taking the 95% credible interval 
# from the 2.5 and 97.5 percentiles, we can say that we are 95% sure that the 
# value of meanlog.true lies in the interval [-1.07, -0.96]. 

# Priors and posteriors for data.big with p.prior4
# Essentially the same as with p.prior3:
plot.dmc(samples.big4,layout=c(2,3),p.prior=p.prior4)

# Compare posteriors obtained with 
# p.prior3 (samples.big3) and p.prior4 (samples.big4) 
# for data.big (note that this code is looping over the arrays that store
# the "theta" parameter estimates in the samples objects):
par(mfcol=c(2,3))
for(i in names(p.vector)){
  #Posteriors obtained with weakly informative p.prior3
  hist(samples.big3$theta[,i,],main=i,freq=F,breaks="fd") 
  #Posteriors obtained with misspecified & diffuse p.prior4
  lines(density(samples.big4$theta[,i,]),col="red",lwd=2) 
}

### Prior matters with scarce data

#data.small <- data.model.dmc(simulate.dmc(p.vector,model,n=15),model)
head(data.small)

# First compare results obtained with p.prior3 (samples.small3) and p.prior4 
# (samples.small4) for data.small

# Priors and posteriors for data.small with p.prior3
# Priors are barely updated:
plot.dmc(samples.small3,layout=c(2,3),p.prior=p.prior3)
p.vector
# meanlog.true meanlog.false  sdlog.true   sdlog.false   t0 
# -1.0           0.0           1.0           1.0         0.2 

# True values are well within 95% CI, but posteriors are now much wider:
summary.dmc(samples.small3)
#...
# 2. Quantiles for each variable:
#
#                 2.5%     25%      50%     75%   97.5%
# meanlog.true  -1.8642 -1.4746 -1.27844 -1.0861 -0.6948
# meanlog.false -0.5696 -0.1749  0.04214  0.3289  1.0244
# sdlog.true     0.9009  1.2072  1.40163  1.6451  2.1863
# sdlog.false    0.7528  1.0285  1.22035  1.4756  2.1347
# t0             0.1609  0.2133  0.22780  0.2365  0.2438

# First compare samples.big3 and samples.small3 (same prior, but different n)
# Less data = more uncertainty:
par(mfcol=c(2,3))
for(i in names(p.vector)){
  # for data.small with p.prior3
  hist(samples.small3$theta[,i,],main=i,freq=F,breaks="fd",ylim=c(0,8))
  # for data.big with p.prior3
  lines(density(samples.big3$theta[,i,]),col="red",lwd=2)               
}

# Priors and posteriors for data.small with p.prior4
# Priors seem updated (but they were very wide to begin with...), 
# but the posteriors are very wide...
plot.dmc(samples.small4,layout=c(2,3),p.prior=p.prior4)

# Now compare samples.big4 and samples.small4 (same prior, but different n)
# Less data = more uncertainty:
par(mfcol=c(2,3))
for(i in names(p.vector)){
  # for data.small with p.prior4
  hist(samples.small4$theta[,i,],main=i,freq=F,breaks="fd",ylim=c(0,8)) 
  # for data.big with p.prior4
  lines(density(samples.big4$theta[,i,]),col="red",lwd=2)               
}

# Finally let's compare posteriors obtained with 
# p.prior3 (samples.small3) and p.prior4 (samples.small4) 
# for data.small:
par(mfcol=c(2,3))
for(i in names(p.vector)){
  # Posteriors obtained with weakly informative p.prior3
  hist(samples.small3$theta[,i,],main=i,freq=F,breaks="fd") 
  # Posteriors obtained with misspecified & diffuse p.prior4
  lines(density(samples.small4$theta[,i,]),col="red",lwd=2) 
}

# Can you explain this particular pattern of differences in posteriors between 
# samples.small3 and samples.small4?
