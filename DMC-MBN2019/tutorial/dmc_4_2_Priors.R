##################  DMC Lesson 4: Multiple Subjects


### Lesson 4.2:  Hierarchical models, priors 

# IN a hierarchical model we assume that each subject's parameters are sampled 
# from population distributions. For an introduction to the benefits of 
# such models see Shiffrin, R., Lee, M., Kim, W., & Wagenmakers, E.-J. (2008). 
# A Survey of Model Evaluation Approaches With a Tutorial on Hierarchical 
# Bayesian Methods. Cognitive Science: a Multidisciplinary Journal, 32(8), 
# 1248–1284. http://doi.org/10.1080/03640210802414826

# Now we have two types of posterior likelihood to consider:
# a) The posterior likelihood of each subject's parameters, conditional 
#    on the data 
# b) The posterior likelihood of the population parameters, conditional on 
#    the estimates of each subject's parameters.
#
# (a) is the same as in the single subject case (but with one instance for each 
# subject), except that the prior is calculated relative to the current 
# estimates of the population parameters (i.e., it changes over MC iterations)
#
# For (b) we need to specify a prior for the population parameters (which does
# not change over iterations). The relevant likelihood is of the current 
# parameters for each subject (which do change over iterations) under the 
# population model.
#
# 
# To illustrate this process, we will first make a model and sample some data.
# So we do not have to worry about bounds we will use the positive 
# parameterization LNR model.

rm(list=ls())
source ("dmc/dmc.R")
load_model ("LNR","lnrPP.R")

# load_data ("dmc_4_2.RData")

# Use the same model as the single-subject LNR example
model <- model.dmc(p.map=list(meanlog="M",sdlog="M",t0="1",st0="1"),
                   match.map=list(M=list(s1=1,s2=2)),constants=c(st0=0),
                   factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="lnr")

# Specify Gaussian population distributions 
p.mu  <- c(meanlog.true=-1,meanlog.false=0,                 # natural scale
           sdlog.true=log(1),sdlog.false=log(1),t0=log(.2)) # log scale

# Use reasonably tight population standard deviations
p.sigma <- c(.2,.2,.2,.2,.1); names(p.sigma) <- names(p.mu)

# Make the population distributions (because they act in the same role as the
# prior in the single subject case we give them the same name)
p.prior <- prior.p.dmc(p1=p.mu,p2=p.sigma,
                       untrans=c(sdlog.true="exp",sdlog.false="exp",t0="exp"))

# Plot the population distributions on the natural scale. Because the standard
# deviations are small the positive skew of the log-scale parameters is hardly
# evident.
par(mfcol=c(2,3)); for (i in names(p.prior)) plot.prior(i,p.prior)

# Simulate 40 subjects with 250 observations each per cell. Note that 
# the parameters for each subject differ as they are sampled from p.prior
ns <- 40
raw.data <- h.simulate.dmc(model,p.prior=p.prior,ns=ns,n=250)
data.model <- data.model.dmc(raw.data,model)

# Plot the s1 data for 10 subjects (s2 will be the same up to sampling error)
par(mfrow=c(2,5)) # Row 1 = subjects 1..5, row2 = 6..10
for (i in 1:10)
  plot.cell.density(data.cell=data.model[[i]][data.model[[i]]$S=="s1",],C="r1")

# Note the variation in accuracy and how that maps to corresponding parameters
ps <- attr(raw.data,"parameters")
print (ps)
# Note that after the model is added parameters are stored with each subject. 
# The following will extract and make them into a data frame similar to the one
# in the raw data.
ps.natural <- ps; ps.natural[,3:5] <- exp(ps[,3:5])
round(ps.natural,2)
# Accuracy is roughly proportional to a measure like d' 
round(apply(ps.natural[1:10,1:2],1,diff)/
        apply(ps.natural[1:10,3:4],1,function(x){sqrt(sum(x^2))}),2)


# Now make a prior for the population parameters. We won’t make use of it at
# this stage, but it will be necessary later if we were to try to sample
# a hierarchical model based on the data we just simulated.

# We assume a normal prior for p.mu and gamma prior for p.sigma.
#
# The p.mu prior standard deviations should be fairly large if we don’t
# know much a priori about population values, and means should be at roughly
# plausible values. 
p.mu.mu <- c(0,0,log(1),log(1),log(0.2)); names(p.mu.mu) <- names(p.mu)
p.mu.sigma <- c(3,3,1,1,1); names(p.mu.sigma) <- names(p.mu)
mu.prior <- prior.p.dmc(p1=p.mu.mu,p2=p.mu.sigma,
                        untrans=c(sdlog.true="exp",sdlog.false="exp",t0="exp"))

par(mfcol=c(2,3)); for (i in names(mu.prior)) plot.prior(i,mu.prior)
#
# For p.sigma we assume shape = 1 (i.e., an exponential distribution) and set 
# scales that are fairly broad. The exponential and, more generally the
# gamma, are positively skewed distributions. You can play with the setting of
# shape and use the plots below to understand what they can do. An important 
# point is that they are always positive, which is necessary for these
# sigma parameters.
p.sigma.shape <- rep(1,5); names(p.sigma.shape) <- names(p.sigma)
p.sigma.scale <- c(1,1,.5,.5,.2)
sigma.prior <- prior.p.dmc(p1=p.sigma.shape,p2=p.sigma.scale,
                           dists=rep("gamma",length(p.sigma)))
# N.B. Does not make sense to untrans these scale parameters

# Again these are quite vague priors
par(mfcol=c(2,3)); for (i in names(sigma.prior)) plot.prior(i,sigma.prior)

# save_data(model, data.model, ps, p.prior, p.mu, p.sigma, mu.prior, sigma.prior,
#           file="dmc_4_2.RData")
