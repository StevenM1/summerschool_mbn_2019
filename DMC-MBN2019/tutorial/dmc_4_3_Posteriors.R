##################  DMC Lesson 4: Multiple Subjects


### Lesson 4.3:  Hierarchical models, posteriors

rm(list=ls())
source ("dmc/dmc.R")
load_model ("LNR","lnrPP.R")
load_data ("dmc_4_2.RData")

### Calculate the posterior log-likelihoods at the two levels 

# (a) The posterior likelihood of the data conditional on the estimates of each 
# subject's parameters is calculated per subject in the same way as the single 
# subject case.

# For example, for the first subject at the true parameters, the likelihood is
sum(log(likelihood.dmc(ps[1,],data.model[[1]])))
# [1] -218.7203

# and the corresponding prior is
sum(log.prior.dmc(ps[1,],p.prior))
# [1] -1.443087

# Which sum to:
log.posterior.dmc(ps[1,],p.prior,data=data.model[[1]])
# [1] -220.1634


# For (b) The posterior likelihood of the population parameters conditional on 
# the estimates of each subject's parameters consist of:

# The sum of the log-likelihoods of the subject's true parameters (ps) 
# under the (true) population distributions (p.prior) is given by
tmp <- h.log.likelihood.dmc(ps=ps,p.prior=p.prior,pp=list(p.mu,p.sigma)) 
# pp specifies the parameters of the population prior

# Values are calculated for each subject, the relevant value is the one summed 
# over subjects
sum(tmp)   
# [1] 55.22315

# The prior summed over population parameters (here the true values) is 
sum(log.prior.dmc(p.mu,mu.prior)) + sum(log.prior.dmc(p.sigma,sigma.prior))
# [1] -5.551741

# Which sum to:
h.log.posterior.dmc(ps,p.prior,
                    pp=list(p.mu,p.sigma),
                    pp.prior=list(mu.prior,sigma.prior))
# [1] 49.67141
#
# NB1: The two types of population parameters (pp) and their priors are 
#      combined into lists for the pp and pp.prior arguments to this function.
#      The entries in these two lists must be in the same order.


# In the h.log.posterior function the parameter values in p.prior have no 
# influence; they are replaced by the values in pp, but p.prior still needs to
# be passed to the function to act as a container for the values in pp. 
# For example:
p1=rep(NA,5); names(p1) <- names(attr(model,"p.vector"))
na.prior <- prior.p.dmc(p1,p2=rep(NA,5))
h.log.posterior.dmc(ps,p.prior=na.prior,
                    pp=list(p.mu,p.sigma),
                    pp.prior=list(mu.prior,sigma.prior))
# [1] 49.67141


# Just as for the data level, profile likelihoods can be plotted at the 
# population level, but now the parameter estimates (ps) act as the data.
# The following plots the profiles for the two population parameters (indexed 
# by the p.num argument) for each parameter type (indexed by the p.name
# argument) given p.prior. 

# Here the convenience function "assign.pp" is used to first assign the values 
# of pp to p.prior. 
p.prior <- assign.pp(pp=list(p.mu,p.sigma),p.prior)

# Use large sample of subject parameters to check that the maxima are near the
# true values (i.e, p.mu and p.sigma)
tmp <- h.simulate.dmc(model,p.prior=p.prior,ns=1000,n=2)
ps.big <- attr(tmp,"parameters")

ylim=c(0,250)
par(mfrow=c(2,5)) # Top row is the location hyper-parameter, bottom row the scale
h.profile.dmc(p.name="meanlog.true",p.num=1,min.p=-2,max.p=1,ps.big,p.prior,ylim=ylim)
h.profile.dmc(p.name="meanlog.false",p.num=1,min.p=-2,max.p=1,ps.big,p.prior,ylim=ylim)
h.profile.dmc(p.name="sdlog.true",p.num=1,min.p=log(0.25),max.p=log(2),ps.big,p.prior,ylim=ylim)
h.profile.dmc(p.name="sdlog.false",p.num=1,min.p=log(0.25),max.p=log(2),ps.big,p.prior,ylim=ylim)
h.profile.dmc(p.name="t0",p.num=1,min.p=log(0.1),max.p=log(0.4),ps.big,p.prior,digits=3,ylim=c(0,1000))
h.profile.dmc(p.name="meanlog.true",p.num=2,min.p=.05,max.p=1,ps.big,p.prior,ylim=ylim)
h.profile.dmc(p.name="meanlog.false",p.num=2,min.p=.05,max.p=1,ps.big,p.prior,ylim=ylim)
h.profile.dmc(p.name="sdlog.true",p.num=2,min.p=.05,max.p=1,ps.big,p.prior,ylim=ylim)
h.profile.dmc(p.name="sdlog.false",p.num=2,min.p=.05,max.p=1,ps.big,p.prior,ylim=ylim)
h.profile.dmc(p.name="t0",p.num=2,min.p=.05,max.p=.4,ps.big,p.prior,ylim=c(0,1000))


