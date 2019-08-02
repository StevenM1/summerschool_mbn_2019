##################  DMC Lesson 3: Sampling

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# THIS LESSON ALSO REQUIRES      # 
#        PACKAGE loo             #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

### Lesson 3.5:  Model Selection

# In this lesson we will examine various "information criteria" that are used
# to perform model selection, favouring the model that provides a good fit while
# also being reasonably simple. The advanced lesson dmc_5_3_ModelSelection.R
# does the same for hierarchical models (best looked at after completing 
# lesson 4).

rm(list=ls()) 
# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")

# Let's compare two DDM models for a design with a stimulus manipulation
# (S) and a second manipulation (A). In model 0 neither factor has an effect,
# In model 1 A has an effect on v. 

# All of the fitting was done beforehand (see commented out code at the end of 
# this file). Modify that code if you want to redo it yourself, but use a 
# smaller n.trials.per.cell, here a large value is used to show (effectively)
# asymptotic results
load_model ("DDM","ddm.R") 
load_data ("dmc_3_5.RData")

# LOOK AT THE DATA
par(mfrow=c(2,4))
plot.cell.density(data.cell=data0[data0$S=="left" & data0$A=="a1",],
                  C="LEFT",xlim=c(0,2),ymax=3,main="0: a1.left")
plot.cell.density(data.cell=data0[data0$S=="right" & data0$A=="a1",],
                  C="RIGHT",xlim=c(0,2),ymax=3,main="0: a1.right")
plot.cell.density(data.cell=data0[data0$S=="left" & data0$A=="a2",],
                  C="LEFT",xlim=c(0,2),ymax=3,main="0: a2.left")
plot.cell.density(data.cell=data0[data0$S=="right" & data0$A=="a2",],
                  C="RIGHT",xlim=c(0,2),ymax=3,main="0: a2.right")
plot.cell.density(data.cell=data1[data1$S=="left" & data1$A=="a1",],
                  C="LEFT",xlim=c(0,2),ymax=3,main="1: a1.left")
plot.cell.density(data.cell=data1[data1$S=="right" & data1$A=="a1",],
                  C="RIGHT",xlim=c(0,2),ymax=3,main="1: a1.right")
plot.cell.density(data.cell=data1[data1$S=="left" & data1$A=="a2",],
                  C="LEFT",xlim=c(0,2),ymax=3,main="1: a2.left")
plot.cell.density(data.cell=data1[data1$S=="right" & data1$A=="a2",],
                  C="RIGHT",xlim=c(0,2),ymax=3,main="1: a2.right")

round(tapply(data0$RT,cbind(data1[,c("A","S")],
                            C=data0$R==toupper(data0$S)),mean),3)
round(tapply(data1$RT,cbind(data1[,c("A","S")],
                            C=data0$R==toupper(data0$S)),mean),3)

# In data1 a2 is ~10% more accurate and ~10ms slower than a1, 
# data0 falls in the middle

# We will now examine model selection statistics to compare the true and false
# models. Note that the effectiveSize (~500, despite being based on 18-21,000
# actual samples!) is probably too small, with several thousand being a better
# minimal amount given the importance of min and var statistics, see below.

# First lets look at the fits
style="cdf"

# For data1 the misfit of the false model is clear
plot.pp.dmc(pp1.true,style=style,layout=c(2,2))
plot.pp.dmc(pp1.false,style=style,layout=c(2,2))

# However for data0 both fit very well
plot.pp.dmc(pp0.true,style=style,layout=c(2,2))
plot.pp.dmc(pp0.false,style=style,layout=c(2,2))

# Model selection should take model complexity into account as well as fit 
# (preferring models that are simpler when they fit "well enough"),
# but getting the balance between fit and simplicity can be tricky.

# Dstats.dmc returns a list of statistics about the posterior deviance, a 
# measure of MISFIT to the data (so SMALLER is BETTER) that can be calculated for 
# each posterior parameter sample. It can return the entire set of samples ($D) 
# with save=TRUE, but by default only returns summaries.

# Get statistics saving the deviances 
true0 <- Dstats.dmc(samples0.true,save=TRUE)
true1 <- Dstats.dmc(samples1.true,save=TRUE)
false0 <- Dstats.dmc(samples0.false,save=TRUE)
false1 <- Dstats.dmc(samples1.false,save=TRUE)

round(unlist(true0[-6]),2)
#       np    meanD     varD     minD    Dmean 
#     6.00 -2180.25    13.60 -2186.21 -2184.79 
round(unlist(false0[-6]),2)
#       np    meanD     varD     minD    Dmean 
#     7.00 -2179.65    12.69 -2186.08 -2185.15 

round(unlist(true1[-6]),2)
#       np    meanD     varD     minD    Dmean 
#     7.00 -2986.08    13.38 -2992.06 -2991.43 
round(unlist(false1[-6]),2)
#       np    meanD     varD     minD    Dmean 
#     6.00 -2614.92    10.79 -2619.97 -2619.21 

# np is number of parameters
# minD is the minimum deviance (best fit)
# Dmean is the deviance of the mean of parameters, it is also an estimate
#   of the best fit, working best if the distribution of D is gaussian. Note 
#   that the more complex model tends to fit better (smaller D). The difference
#   is much bigger for data1 where there is extra structure to fit.
# meanD is the mean of the posterior deviances. By itself it can be used as a 
#   model selection measure that (somewhat) takes account of model complexity
#   (see Liu, C. C., Aitkin, M. (2008). Bayes factors: Prior sensitivity and model 
#   generalizability, 52(6), 362-375.), although a method based on the full
#   distribution of D is preferred by Aitkin (see below). Note that unlike minD and
#   Dmean it is the true models which win on this measure. 
# varD is the variance of the posterior deviance, half of this value was 
#   suggested as a measure of the number of effective parameters (pd) by
#   Sturtz, S., Ligges, U., & Gelman, A. (2005). R2WinBUGS: A package for running
#   WinBUGS from R. Journal of Statistical Software, 12, 1-16.

# pd.dmc provides three estimates of pd (the number of parameters taking 
# account of correlations that modify the degree to which they act 
# independently). 

round(unlist(pd.dmc(true0)),2)
# Pmean  Pmin  Pvar 
#  4.54  5.96  6.80 
round(unlist(pd.dmc(false0)),2)
# Pmean  Pmin  Pvar 
#  5.50  6.43  6.35 

round(unlist(pd.dmc(true1)),2)
# Pmean  Pmin  Pvar 
#  5.35  5.99  6.69 
round(unlist(pd.dmc(false1)),2)
# Pmean  Pmin  Pvar 
#  4.29  5.05  5.40 

# Two estimates of pd are based on the distance between meanD and an estimate
# of the best fit, either Dmean or minD
# Pmean = meanD - Dmean
# Pmin  = meanD - minD
# All three measures roughly correspond to the idea that greater variability in
# D indicates a more flexible model. plot.deviance.dmc shows this graphically.
# The function can be used without first calling Dstats using samples=argument.

par(mfrow=c(2,2))
plot.deviance.dmc(true0,main="true0",xlim=c(-2190,-2150))
plot.deviance.dmc(true1,main="true1") 
plot.deviance.dmc(false0,main="false0",xlim=c(-2190,-2150))
plot.deviance.dmc(false1,main="false1") 

# Note that the distributions are positively skewed, so it is worth thinking about
# whether to use minD or Dmean. Also the long tail can make estimates of the
# mean and especially variance rather unstable. All of these facts underline
# the need for a fairly large effective sample size to do model selection.

# A posterior likelihood ratio test can be calculated by taking the difference
# in the posterior deviances from two models D1-D2 (NB: the deviance is 
# proportional to the log-likelihood, so the difference is the log of the 
# likelihood ratio). The test produces the probability that D1 wins (D1 < D2).
# Aitkin, M., Boys, R. J., & Chadwick, T. (2005). Bayesian point null hypothesis 
#   testing via the posterior likelihood ratio. Statistics and Computing, 15(3), 
#   217-230.

# The test can be performed as follows, with graphical output, using either the
# list output of Dstats.dmc or raw deviances. If the two models have different
# numbers of deviances the longer one is truncated.

# Unsurprisingly the true model wins clearly for data1
posterior.lr.dmc(true1,false1,plot=TRUE)
# pD1 
#   1 

# However the LR test produces only weak evidence for the true model for data0,
# reflecting the fact that the complexity correction is not sufficient.
posterior.lr.dmc(true0,false0,plot=TRUE)
#       pD1 
# 0.5567778 

# The Deviance Information Criterion (DIC) is the most popular model selection
# method, where: DIC = meanD + pd, as proposed by:
# Spiegelhalter, D. J., Best, N. G., Carlin, B. P., & van der Linde, A. (2002). 
#   Bayesian measures of model complexity and fit. Journal of the Royal 
#   Statistical Society: B, 64, 583-639.
# A little algebra shows DIC = Dmean + 2pd (or minD + 2pd), which is why it is
# sometimes compared to the Akaike Information Criteria (AIC), D + 2p, used with
# maximum likelihood estimation. Like AIC it tends to pick models that are too
# complex asymptotically.

# The function IC.dmc calculates DIC in various ways. use.pd can specify which
# pd calculated by pd.dmc to use (e.g.,Pmean,Pmin,Pvar) but by default takes
# the best of Pmin and Pmean (i.e,. whichever of minD and Dmean, calculated by
# Dstats.dmc, is smaller and therefore a better estimate of the minimum). Here, 
# this is Pmin.

# Again unsurprisingly DIC works well with data1
IC.dmc(true1,DIC=TRUE) 
# [1] -2980.093
IC.dmc(false1,DIC=TRUE) 
# [1] -2609.874


# However data0 is more challenging
IC.dmc(true0, DIC=TRUE) # True model
# [1] -2174.287
IC.dmc(false0,DIC=TRUE) # Model too complex
# [1] -2173.219

# There are various guidelines to what DIC difference is substantial. A good
# guide is provided by DIC weights, which are analogous to AIC weights, see
# Wagenmakers, E. J., & Farrell, S. (2004). AIC model selection using Akaike 
#  weights. Psychonomic Bulletin & Review, 11(1), 192-196.
# These can be interpreted as probabilities that a model is best. They are of 
# the following form, where d is the difference between ICs. For example:
d=c(0,6) # the first model is the best, the other is worse by 6
exp(-d/2)/sum(exp(-d/2))
# So a difference of 6 means you are >95% in favour of a model.
d=c(0,10) # the first model is the best, the other is worse by 10
exp(-d/2)/sum(exp(-d/2))
# is about >99%. The following function is set up to compare any number of 
# models so pass results as a list. Again the result is strong for data1
wIC.dmc(list(true1,false1),DIC=TRUE)
#        M1           M2
# IC-min  0 3.702186e+02
# w       1 4.055573e-81
# Here "IC-min" is d above. Not so clear for data0 but clearer than the LR test
wIC.dmc(list(true0,false0),DIC=TRUE,digits=3)
#          M1   M2
# IC-min 0.00 1.07
# w      0.63 0.37

# The following paper suggested a correction to DIC which is simple for 
# non-hierarchical models, BPIC = meanD + 2pd (note the bigger penalty, recall
# that DIC = Dmean + 2pd, which with some algebra becomes DIC = meanD + pd)
# Ando, T. (2007). Bayesian predictive information criterion for the evaluation 
#   of hierarchical Bayesian and empirical Bayes models. Biometrika, 94(2), 
#   443-458. 
# A later paper showed it is approximately correct even for hierarchical models.
# Ando, T. (2011). Predictive Bayesian Model Selection. American Journal of 
#   Mathematical and Management Sciences, 31(1-2), 13-38. 
# dmc uses this as its default measure. For the difficult data0 case
wIC.dmc(list(true0,false0))
#               M1        M2
# IC-min 0.0000000 1.5393735
# w      0.6834531 0.3165469
# Some conventions consider IC differences of more than 3 (or 3.2) as moderate
# evidence, and whatever the convention in this case at least this seems best.


# The model selection approaches examined so far are not "fully Bayesian", in
# that they rely on point estimates. Two alternative information criteria that 
# avoid this problem are provided by the "loo" package. loo stands for 
# "leave-out one", a type of cross-validation which involves fitting to all but
# one data point and then assessing predictive performance for the left out data
# point. Literally doing this for every data point is computationally 
# prohibitive, so loo provides two approximations. The following paper describes
# and tests them relative to full Bayesian leave-out-one cross-validation as the
# gold standard.
# Vehtari, A., Gelman, A., and Gabry, J. (2015). Efficient implementation of 
# leave-one-out crossvalidation and WAIC for evaluating fitted Bayesian models.

# Both methods require the "pointwise" log-likelihoods for each data point. By 
# default DMC only stores their sum over data points (see lesson 3_1) as the 
# full iterations x chains x data points array can be very large (e.g., ~5.5 GB
# for each of the models examined here!). 

# The size can be reduced by "thinning"; using the log-likelihood for only each 
# "thin-pointwise" iteration with little loss of information due to the high 
# chain autocorrelation. To chose a value of thinning it is useful to plot the
# autocorrelation function for the summed posterior likelihood. This can be done
# with the same function used for posterior parameter samples but setting 
# par=NA. 
acf.dmc(samples0.true,par=NA,chain=1,lag.max=150)
acf.dmc(samples0.false,par=NA,chain=1,lag.max=150)
acf.dmc(samples1.true,par=NA,chain=1,lag.max=150)
acf.dmc(samples0.false,par=NA,chain=1,lag.max=150)

# Here we set a large value of the maximum lag to show that that the 
# auto-correlation is extreme, remaining significant beyond lag 100, which
# reinforces the point that there are not nearly enough samples here for 
# stable model selection.

# Given these considerations we did not store the pointwise log-likelihoods
# but rather calculate them based on the stored parameter estimates and data
# using the "trial_log_likes" function, and set thin-pointwise=100 (a value of
# 10 produced similar results but with much larger objects). Because
# the thinning is quite extreme (yielding only 10 values per chain, so 180 
# values for each of the 40,000 data points), this takes only ~15sec each to 
# run, but can be slow for more realistic sample sizes. 
#   NB: These objects are not stored in "dmc_3_5.RData" as even with extreme 
#       thinning they total ~250Mb!)

samples0.true.pll <- trial_log_likes(samples0.true,thin_pointwise = 100)
samples0.false.pll <- trial_log_likes(samples0.false,thin_pointwise = 100)
samples1.true.pll <- trial_log_likes(samples1.true,thin_pointwise = 100)
samples1.false.pll <- trial_log_likes(samples1.false,thin_pointwise = 100)

# The first method, WAIC (the widely applicable or Watanabe-Akaike information 
# criterion) is easy to compute but can be distorted by outlier data values,
# although that is mainly a problem only in small data samples. 

# The "waic.dmc" function can be used to calculate WAIC and an associated 
# estimate of number of effective parameters (p). For both it calculates an
# associated standard error, with asymptotic normality holding as the number
# of data points grows large. The following se_waic estimates highlight the
# substantial uncertainty about the waic values due to the small (effective)
# sample size.

samples0.true.waic <- waic.dmc(samples0.true.pll,digits=2,save=TRUE)
#        p     se_p     waic  se_waic 
#     5.60     0.21 -2174.73   419.24 
samples0.false.waic <- waic.dmc(samples0.false.pll,digits=2,save=TRUE)
#        p     se_p     waic  se_waic 
#     6.91     0.22 -2172.48   419.57 
samples1.true.waic <- waic.dmc(samples1.true.pll,digits=2,save=TRUE)
#        p     se_p     waic  se_waic 
#     5.96     0.13 -2980.33   421.09 
samples1.false.waic <- waic.dmc(samples1.false.pll,digits=2,save=TRUE)
#        p     se_p     waic  se_waic 
#     5.37     0.12 -2609.43   419.84 

# waic.dmc can also store the list created by the underlying waic 
# function from the loo package. This includes a $pointwise matrix which
# includes the contributions of each data point to waic.

# The "loocompare.dmc" function can compare the outputs of waic.dmc for two
# models, taking the difference (second - first) and calculating a standard
# error that takes advantage of the correlation in the two waic estimates due
# to each being calculated on the same data. See ?compare for the loo package
# discussion of the standard errors

# Despite the large amount of uncertainty in the waic values the correlation
# results in a strong indication that the waic_diff(erence) strongly favors
# true model 1 (the first model is favored by a positive waic_diff)

loocompare.dmc(samples1.true.waic,samples1.false.waic,digits=3)
# waic_diff        se 
#     370.9      38.6 

# For true model 0 the support is more modest
loocompare.dmc(samples0.true.waic,samples0.false.waic,digits=3)
# waic_diff        se 
#      2.26      1.07 


# loocompare can also be used to calculate model weights by passing the
# outputs of waic.dmc as a list (this will work with more than two objects).

# As in earlier approaches the support for the true model 1 is strong
loocompare.dmc(list(true1=samples1.true.waic,false1=samples1.false.waic),digits=3)
#        true1   false1
# IC-min     0 3.71e+02
# w          1 2.88e-81

# It is less strong for true model 0 (although better than any other approach)
loocompare.dmc(list(true0=samples0.true.waic,false0=samples0.false.waic),digits=2)
#        true0 false0
# IC-min  0.00   2.26
# w       0.76   0.24


# The second information criterion, looic, is a better approximation but takes
# a little more computing. It uses "Pareto smoothed importance sampling", which
# requires estimation of a Pareto distribution shape parameter (k) for each 
# data point. Values of k < 0.5 are desirable with 0.5<k<1 indicating potential
# problems and k>1 more severe problems, so warning messages are printed when
# this occurs. The "looic.dmc" function can be used to calculate looic and 
# store the object created by the underlying loo function from the loo package.

samples0.true.looic <- looic.dmc(samples0.true.pll,digits=2,save=TRUE)
#        p     se_p    looic se_looic 
#     5.57     0.21 -2174.80   419.24 
# Warning message:
# In looic.dmc(samples0.true.pll, digits = 2, save = TRUE) :
#   129 (0%) Pareto k estimates between 0.5 and 1
# See PSIS-LOO description (?'loo-package') for more information
samples0.false.looic <- looic.dmc(samples0.false.pll,digits=2,save=TRUE)
# All Pareto k estimates OK (k < 0.5)
#        p     se_p    looic se_looic 
#     6.87     0.22 -2172.56   419.56 
samples1.true.looic <- looic.dmc(samples1.true.pll,digits=2,save=TRUE)
# All Pareto k estimates OK (k < 0.5)
#        p     se_p    looic se_looic 
#     5.93     0.13 -2980.40   421.09 
samples1.false.looic <- looic.dmc(samples1.false.pll,digits=2,save=TRUE)
# All Pareto k estimates OK (k < 0.5)
#        p     se_p    looic se_looic 
#     5.33     0.12 -2609.50   419.84 

# samples0.true indicates 129/40000 data points produced mildly problematic k
# values, likely because of the very small sample size (180).

# loocompare.dmc can operate on the saved looic.dmc objects in just the same
# way as for waic.dmc saved objects.


# Another way of testing if it is necessary to have a difference between v for
# a1 and a2 is to do posterior predictive tests to see if the credible interval
# for their sampled difference contains zero. Recall that in model 1 v.a1=.75 
# and v.a2=1.25, so we would expect v.a2-v.a1 = 0.5, whereas in model 0 they are
# the same so we would expect no difference.

vdiff <- function(p) {p["v.a2"]-p["v.a1"]}

vdiffs1 <- p.fun.dmc(samples=samples1.true,fun=vdiff)
round(c(mean(vdiffs1),quantile(vdiffs1,probs=c(.025,.975))),2)
#        2.5% 97.5% 
#  0.47  0.42  0.52 

vdiffs0 <- p.fun.dmc(samples=samples0.false,fun=vdiff)
round(c(mean(vdiffs0),quantile(vdiffs0,probs=c(.025,.975))),2)
#        2.5% 97.5% 
#  0.01 -0.04  0.06 


# ########################## DATA AND FIT GENERATION
# 
# n.trials.per.cell = 1e4
# 
# model0 <- model.dmc(constants=c(st0=0,d=0),type="rd",
#   p.map=list(a="1",v="1",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
#   match.map=list(M=list(left="LEFT",right="RIGHT")),responses=c("LEFT","RIGHT"),
#   factors=list(S=c("left","right"),A=c("a1","a2")))
# p.vector0  <- c(a=1,v=1,z=0.5,sv=1,sz=0.2,t0=.15)
# data0 <- data.model.dmc(simulate.dmc(p.vector0,model0,n=n.trials.per.cell),model0)
# CHECK THE LIKELIHOODS, SZ NOT WELL DEFINED
# par(mfrow=c(2,3)); ylim = c(0,1000)
# profile.dmc("a",.1,2,p.vector0,data0,ylim=ylim)
# profile.dmc("v",.1,2,p.vector0,data0,ylim=ylim)
# profile.dmc("z",.2,.8,p.vector0,data0,ylim=ylim)
# profile.dmc("sv",.1,2,p.vector0,data0,ylim=ylim)
# profile.dmc("sz",.1,.9,p.vector0,data0,ylim=ylim)
# profile.dmc("t0",.01,.9,p.vector0,data0,ylim=ylim)
# p.prior0 <- prior.p.dmc(
#   dists = c("tnorm","tnorm","beta","tnorm","beta","beta"),
#   p1=c(a=1,v=0,z=1,sv=1,sz=1,t0=1),                           
#   p2=c(a=1,v=2,z=1,sv=1,sz=1,t0=1),
#   lower=c(0,-5,NA,0,NA,NA),
#   upper=c(2, 5,NA,2,NA,NA)
# )
# model1 <- model.dmc(constants=c(st0=0,d=0),type="rd",
#   p.map=list(a="1",v="A",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
#   match.map=list(M=list(left="LEFT",right="RIGHT")),responses=c("LEFT","RIGHT"),
#   factors=list(S=c("left","right"),A=c("a1","a2")))
# 
# p.vector1  <- c(a=1,v.a1=.75,v.a2=1.25,z=0.5,sv=1,sz=0.2,t0=.15)
# data1 <- data.model.dmc(simulate.dmc(p.vector1,model1,n=n.trials.per.cell),model1)
# # LOOK AT THE DATA
# # AGAIN, SZ NOT WELL DEFINED
# par(mfrow=c(2,4)); ylim = c(0,2000)
# profile.dmc("a",.1,2,p.vector1,data1,ylim=ylim)
# profile.dmc("v.a1",.1,2,p.vector1,data1,ylim=ylim)
# profile.dmc("v.a2",.1,2,p.vector1,data1,ylim=ylim)
# profile.dmc("z",.2,.8,p.vector1,data1,ylim=ylim)
# profile.dmc("sv",.1,2,p.vector1,data1,ylim=ylim)
# profile.dmc("sz",.1,.9,p.vector1,data1,ylim=ylim)
# profile.dmc("t0",.01,.9,p.vector1,data1,ylim=ylim)
# p.prior1 <- prior.p.dmc(
#   dists = c("tnorm","tnorm","tnorm","beta","tnorm","beta","beta"),
#   p1=c(a=1,v.a1=0,v.a2=0,z=1,sv=1,sz=1,t0=1),                           
#   p2=c(a=1,v.a1=2,v.a2=2,z=1,sv=1,sz=1,t0=1),
#   lower=c(0,-5,-5,NA,0,NA,NA),
#   upper=c(2, 5,5,NA,2,NA,NA)
# )
# 
# # TRUE MODEL FITS
# 
# # FITTING MODEL0 TO MODEL0 DATA 
# samples0.1 <- samples.dmc(nmc=500,p.prior0,data0)
# samples0.1 <- run.dmc(samples0.1, report = 25, cores=4,p.migrate=.05)
# plot.dmc(samples0.1,layout=c(3,4))
# samples0.2 <- run.dmc(samples.dmc(nmc=500,samples=samples0.1), 
#                       cores=4,report=25)
# plot.dmc(samples0.2,layout=c(3,4))
# samples0.3 <- run.dmc(samples.dmc(nmc=500,samples=samples0.2,add=TRUE), 
#                       cores=4,report=25)
# plot.dmc(samples0.3,layout=c(3,4))
# gelman.diag(theta.as.mcmc.list(samples0.3),transform=TRUE)
# effectiveSize(theta.as.mcmc.list(samples0.3))
# summary.dmc(samples0.3)
# pairs.dmc(samples0.3)
# 
# pp0 <- post.predict.dmc(samples=samples0.3)
# plot.pp.dmc(pp0) 
# plot.pp.dmc(pp0,"cdf")
# # add in pointwise likelihood
# samples0.4 <- run.dmc(
#   tmp=samples.dmc(nmc=10,samples=samples0.3,pointwise_likelihood = TRUE)
#   ,cores=4,report=1)
# 
# 
# # FITTING MODEL1 TO MODEL1 DATA 
# samples1.1 <- samples.dmc(nmc=500,p.prior1,data1)
# samples1.1 <- run.dmc(samples1.1, report = 10, cores=4,p.migrate=.05)
# plot.dmc(samples1.1,layout=c(4,4))
# samples1.2 <- run.dmc(samples.dmc(nmc=500,samples=samples1.1), 
#                     cores=4,report=25)
# plot.dmc(samples1.2,layout=c(4,4))
# samples1.3 <- run.dmc(samples.dmc(nmc=500,samples=samples1.2,add=TRUE), 
#                       cores=4,report=25)
# summary.dmc(samples1.3)
# effectiveSize(theta.as.mcmc.list(samples1.3))
# gelman.diag(theta.as.mcmc.list(samples1.3),transform=TRUE)
# pp1 <- post.predict.dmc(samples1.3)
# plot.pp.dmc(pp1)
# plot.pp.dmc(pp1,"cdf")
# pairs.dmc(samples1.3)
# 
# # FALSE MODEL FITS
# 
# # FITTING MODEL1 TO MODEL0 DATA 
# 
# data10 <- data.model.dmc(data0,model1)
# samples10.1 <- samples.dmc(nmc=500,p.prior=p.prior1,data=data10)
# samples10.1 <- run.dmc(samples10.1, report = 10, cores=4,p.migrate=.05)
# plot.dmc(samples10.1,layout=c(4,4))
# samples10.2 <- run.dmc(samples.dmc(nmc=500,samples=samples10.1), 
#                     cores=4,report=25)
# plot.dmc(samples10.2,layout=c(4,4))
# gelman.diag(theta.as.mcmc.list(samples10.2),transform=TRUE)
# samples10.3 <- run.dmc(samples.dmc(nmc=500,samples=samples10.2,add=TRUE), 
#                       cores=4,report=25)
# # sz not good, try a fresh batch to see if it can be improved.
# gelman.diag(theta.as.mcmc.list(samples10.3),transform=TRUE)
# samples10.4 <- run.dmc(samples.dmc(nmc=1000,samples=samples10.3), 
#                        cores=4,report=25)
# # Improved a bit
# gelman.diag(theta.as.mcmc.list(samples10.4),transform=TRUE)
# plot.dmc(samples10.4,layout=c(4,4))
# summary.dmc(samples10.4)
# effectiveSize(theta.as.mcmc.list(samples10.4))
# pp10 <- post.predict.dmc(samples10.4)
# plot.pp.dmc(pp10)
# plot.pp.dmc(pp10,"cdf")
# pairs.dmc(samples10.4)
# 
# # FITTING MODEL0 TO MODEL1 DATA
# data01 <- data.model.dmc(data1,model0)
# samples01.1 <- samples.dmc(nmc=500,p.prior0,data01)
# samples01.1 <- run.dmc(samples01.1, report = 25, cores=4,p.migrate=.05)
# plot.dmc(samples01.1,layout=c(3,4))
# samples01.2 <- run.dmc(samples.dmc(nmc=500,samples=samples01.1), 
#                       cores=4,report=25)
# plot.dmc(samples01.2,layout=c(3,4))
# samples01.3 <- run.dmc(samples.dmc(nmc=500,samples=samples01.2,add=TRUE), 
#                       cores=4,report=25)
# # sz not good, try a fresh batch to see if it can be improved.
# gelman.diag(theta.as.mcmc.list(samples01.3),transform=TRUE)
# samples01.4 <- run.dmc(samples.dmc(nmc=1000,samples=samples01.3), 
#                        cores=4,report=25)
# # Improved a bit
# gelman.diag(theta.as.mcmc.list(samples01.4),transform=TRUE)
# plot.dmc(samples01.4,layout=c(3,4))
# effectiveSize(theta.as.mcmc.list(samples01.4))
# summary.dmc(samples01.4)
# pp01 <- post.predict.dmc(samples01.4)
# plot.pp.dmc(pp01)
# plot.pp.dmc(pp01,"cdf")
# pairs.dmc(samples01.4)
# 
# save_data (p.vector0,data0,p.prior0,samples0.1,samples0.2,samples0.3,pp0,
#      p.vector1,data1,p.prior1,samples1.1,samples1.2,samples1.3,pp1,
#      samples10.1,samples10.2,samples10.3,pp10,samples10.4,pp10,
#      samples01.1,samples01.2,samples01.3,pp01,samples01.4,pp01,
#      file="tutorial/data/dmc_3_5_setup.RData")
# 
# # FINAL SET OF SAMPLES
# samples00 <- samples0.3
# samples11 <- samples1.3
# samples10 <- samples10.4
# samples01 <- samples01.4
# pp00 <- pp0
# pp11 <- pp1
# samples0.true <- samples00
# samples0.false <- samples10
# samples1.true <- samples11
# samples1.false <- samples01
# pp0.true <- pp00
# pp0.false <- pp10
# pp1.true <- pp11
# pp1.false <- pp01
# 
# save_data (data0,data1,samples0.true,samples0.false,samples1.true,
#            samples1.false,pp0.true,pp0.false,pp1.true,pp1.false,
#            file="dmc_3_5.RData")
# 


