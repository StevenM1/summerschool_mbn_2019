##################  DMC Lesson 4: Multiple Subjects


### Lesson 4.5:  Multiple subjects as random effects, LNR model

# Same setup as lesson 4.2, 40 subjects

rm(list=ls())

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("LNR","lnrPP.R")

load_data ("dmc_4_4.RData")
# load_data ("dmc_4_5.RData")

# From lesson 4.2 recall that we can use prior.p.dmc to make priors for the 
# population level hyper-parameters. These settings are quite vague.

# Location prior
p.mu.mu <- c(0,0,log(1),log(1),log(0.2)); names(p.mu.mu) <- names(p.mu)
p.mu.sigma <- c(3,3,1,1,1); names(p.mu.sigma) <- names(p.mu)
mu.prior <- prior.p.dmc(p1=p.mu.mu,p2=p.mu.sigma,
                        untrans=c(sdlog.true="exp",sdlog.false="exp",t0="exp"))
# Scale prior
p.sigma.shape <- rep(1,5); names(p.sigma.shape) <- names(p.sigma)
p.sigma.scale <- c(1,1,.5,.5,.2)
sigma.prior <- prior.p.dmc(p1=p.sigma.shape,p2=p.sigma.scale,
                           dists=rep("gamma",length(p.sigma)))

# Make a hyper-prior list
pp.prior=list(mu.prior,sigma.prior)


# Make a new samples object from the same 40 subject data as in the last lesson
# Specifying pp.prior tells h.samples.dmc to add an attribute "hyper" to samples 
# object to enable hierarchical sampling.
hsamples <- h.samples.dmc(nmc=100,p.prior,data.model,thin=10,pp.prior=pp.prior)


# Run sampling allowing migration at both levels.
# NB: You can separately tune gamma for hierarchical (h.gamma.mult) and subject 
# (gamma.mult) levels. By default h.gamma.mult = NA, so gamma is set as a
# uniform sample on 0.5-1 (this may not always be suitable).

# On the second set of samples all stuck chains are gone.
hsamples <- h.run.unstuck.dmc(hsamples,cores=4,p.migrate=0.05,h.p.migrate=0.05)

# Now run with migration turned off in order to get proper samples from the 
# posterior. Convergence is immediate (verbose=TRUE, the default, causes the 
# report below to be printed).
hsamples1 <- h.run.converge.dmc(h.samples.dmc(nmc=100,samples=hsamples),
                                nmc=50,cores=4,verbose=TRUE)
# Iterations = 100, Effective N = NA, Hyper mpsrf = 1.02
# Subject mpsrf achieved (sorted):
#   28    9   22   14   24   13    7    1   30   27   31   37   33   39    2   17   26   35 
# 1.05 1.05 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.03 
#   21   16   23   40   10    4   15   25   19   12    8    5   36   18   20   34   38   11 
# 1.03 1.03 1.03 1.03 1.03 1.03 1.03 1.03 1.03 1.03 1.03 1.03 1.03 1.03 1.03 1.03 1.03 1.03 
#    3   32    6   29 
# 1.02 1.02 1.02 1.02 

# Normally it can take a number of trys to get convergence (the default maximum
# is 20). On each try nmc iterations are added, and corresponding hyper multivariate 
# gelman.diag computed. Thye hyper multivarite gelman.diag is then computed for the new 
# series but with the first nmc iterations removed, and whichever is best is kept. 
# In the initial stages the earliest iterations are typically repeatedly discarded,
# then towards the end the series length grows untill there are sufficient 
# samples (recall gelman.diag can be high simply because the series is short).

# The following convenience function returns the Multivariate psrfs at the 
# subject level and the hyper level as a vector.
h.gelman.diag.dmc(hsamples1)
#  1.02  1.01  1.01  1.01  1.01  1.01  1.01  1.01  1.01  1.01  1.01  1.01 
#    36    16    40     6    18    27    34    11    12    31    15    35 
#  1.01  1.01  1.01  1.01  1.01  1.02  1.02  1.02  1.02  1.02  1.02  1.02 
#    14    25    17    24    22    21     4    30    37    13    23    33 
#  1.02  1.02  1.02  1.02  1.02  1.02  1.02  1.02  1.02  1.02  1.02  1.02 
#    28     1     7     9    26 
#  1.03  1.03  1.03  1.03  1.03 

# More detail on the hyper level convergence can be got as follows 
gelman.diag.dmc(hsamples1,hyper=TRUE)
#                  Point est. Upper C.I.
# meanlog.true.h1       1.003       1.01
# meanlog.false.h1      1.004       1.01
# sdlog.true.h1         0.998       1.00
# sdlog.false.h1        1.008       1.02
# t0.h1                 1.008       1.02
# meanlog.true.h2       1.004       1.01
# meanlog.false.h2      1.012       1.02
# sdlog.true.h2         1.004       1.01
# sdlog.false.h2        1.003       1.01
# t0.h2                 1.007       1.02
# 
# Multivariate psrf
# 
# 1.02

# We can check the effective number of hyper samples, which are close to the
# nominal values due to thinning.
effectiveSize.dmc(hsamples1,hyper=TRUE)
#  meanlog.true.h1 meanlog.false.h1    sdlog.true.h1   sdlog.false.h1 
#             1140             1172             1417             1343 
#            t0.h1  meanlog.true.h2 meanlog.false.h2    sdlog.true.h2 
#             1395             1154             1152             1232 
#   sdlog.false.h2            t0.h2 
#             1699             1494 

# This is reflected in their being no autocorrelation at the hyper.
acf.dmc(hsamples1,hyper=TRUE)
# or subject level (here for subject 1)
acf.dmc(hsamples1[[1]])

# Both the posterior log-likelihoods and the parameter chains are very nice.
plot.dmc(hsamples1,hyper=TRUE,pll.chain=TRUE)
plot.dmc(hsamples1,hyper=TRUE,layout=c(2,5))

# The credible intervals on the hyper estimates are given by
summary.dmc(hsamples1,hyper=TRUE)
#                     2.5%      25%       50%      75%     97.5%
# meanlog.true.h1  -1.08658 -1.03811 -1.011827 -0.98607 -0.933893
# meanlog.false.h1 -0.11381 -0.07303 -0.052249 -0.03197  0.005648
# sdlog.true.h1    -0.02432  0.02201  0.045373  0.06791  0.112641
# sdlog.false.h1   -0.06654 -0.01838  0.004257  0.02864  0.068884
# t0.h1            -1.61607 -1.59464 -1.584406 -1.57411 -1.553911
# meanlog.true.h2   0.19428  0.22366  0.243516  0.26407  0.312439
# meanlog.false.h2  0.13914  0.16385  0.178030  0.19523  0.234573
# sdlog.true.h2     0.16938  0.19308  0.208564  0.22693  0.274546
# sdlog.false.h2    0.16323  0.18847  0.204424  0.22182  0.260719
# t0.h2             0.07444  0.08633  0.093450  0.10223  0.124327

# Updating of the hyper-priors is excellent, although clearly the scale 
# parameters are less updated.
plot.dmc(hsamples1,hyper=TRUE,p.prior=pp.prior,layout=c(2,5))

# The median estimates closely match the sample values
ps <- attr(raw.data,"parameters")
round(apply(ps,2,mean),2)
#  meanlog.true meanlog.false    sdlog.true   sdlog.false            t0 
#         -1.02         -0.05          0.04          0.00         -1.59 
round(apply(ps,2,sd),2)
#  meanlog.true meanlog.false    sdlog.true   sdlog.false            t0 
#          0.24          0.17          0.21          0.21          0.09 

# As in the fixed effect case the average fit is excellent.
pp <- h.post.predict.dmc(hsamples1)


### FIX ME FIX ME
plot.pp.dmc(pp)
# and individual subject fits good considering sample size.
plot.pp.dmc(pp[[1]])
plot.pp.dmc(pp[[1]],style="cdf")

# You can check out posterior correlations for a single subject
pairs.dmc(hsamples1[[1]],thin=5)

# Or for all subjects collapsed (note the larger thin, so in there are only 
# 15 samples for each subject)
pairs.dmc(hsamples1,collapse.subject=TRUE,thin=100)

# By default each subjects parameters are standardized before they are collapsed
# (i.e., divided by SD after mean is subtracted). If you dont want this:
pairs.dmc(hsamples1,collapse.subject=TRUE,scale.subjects=FALSE,thin=100)

# Or you can look at the hyper parameters
pairs.dmc(hsamples1,hyper=TRUE,thin=5)

# You can plot only the hyper location parameters using
pairs.dmc(hsamples1,location=TRUE,thin=5)

# or only the hyper scale paraemters
pairs.dmc(hsamples1,scale=TRUE,thin=5)


# save_data (pp.prior,hsamples,hsamples1,pp,file="dmc_4_5.RData")

#### EXTRAS (NB None of this is saved.)

# In this case all of the sampling went smoothly, but that may not always be the
# case. Sometimes h.run.converge.dmc will appear to converge and the gelman.diag
# will look good but visual inspection will show some movement in the early part
# of the run. In that case you might use h.samples.dmc to selectivley cull samples
# using the "remove" argument, which specifies a range of samples that will be
# omitted, just as was the case for samples.dmc. For example: 

h.samples.short <- h.samples.dmc(nmc=0,samples=hsamples1,add=TRUE,remove=1:25)

# Note that you can even remove non-contiguous sections, so, for example, you 
# could use remove=c(1:10,21:30). 

# Often you might want to combine some removal with making the final result
# longer, h.run.converge.dmc can be used to do this with its finalrun=TRUE 
# argument. For example, you could specify that you want a final series of 
# length 110 as follows after removing the first 10. 

hsamples2 <- h.run.converge.dmc(samples=hsamples1,cores=4,verbose=TRUE,report=1,
                                finalrun=TRUE,addtofinal=TRUE,finalI=110,removefromfinal=1:10)
# As this was previously converged with the same settings the convergence test
# is immediately passed and sampling commences from iteration 91.

# Suppose you want a particular effective number of hyper samples. In that case
# you can use either the finalminN argument (so the worst case parameter achieves),
# the minimum, or finalmeanN so the average is good. For example

effectiveSize.dmc(hsamples2,hyper=TRUE)
#  meanlog.true.h1 meanlog.false.h1    sdlog.true.h1   sdlog.false.h1 
#             1197             1341             1370             1376 
#            t0.h1  meanlog.true.h2 meanlog.false.h2    sdlog.true.h2 
#             1420             1316             1355             1317 
#   sdlog.false.h2            t0.h2 
#             1342             1475 

# Suppose you want a minimum of 1200. Note you have specify nmc here to give
# the number of iterations added for each try (here just 2).
hsamples2 <- h.run.converge.dmc(samples=hsamples2,cores=15,verbose=TRUE,report=1,
                                finalrun=TRUE,addtofinal=TRUE,finalminN=1200,nmc=2)
# ...
# 
# Doing final run
# 111  112  
# Iterations = 112, Effective N = 1210

# Or say you wanted an aveage of at least 1400. Given
# mean(effectiveSize.dmc(hsamples2,hyper=TRUE))
# [1] 1371.5

hsamples2 <- h.run.converge.dmc(samples=hsamples2,cores=4,verbose=TRUE,report=1,
                                finalrun=TRUE,addtofinal=TRUE,finalmeanN=1400,nmc=5)

# All of this could be put together. For example, you might want to do the inital
# converence run then get 200 fresh samples (so addtofinal=FALSE, the default):

hsamples1 <- h.run.converge.dmc(h.samples.dmc(nmc=100,samples=hsamples),
                                nmc=50,cores=4,verbose=TRUE,finalrun=TRUE,finalI=200)

# In the initial convergence process there is also "thorough" option that does
# checks on the basis of every possible univariate gelman.diag, both for the hyper
gelman.diag.dmc(hsamples1,hyper=TRUE)$psrf[,1]
# but also for each subject
unlist(lapply(gelman.diag.dmc(hsamples1),function(x){x$psrf[,1]}))

# In this case it turns out the thorough test has already been satisfied. 
hsamples1.thorough <- h.run.converge.dmc(samples=hsamples1,cores=15,verbose=TRUE,report=1)

# Note, however, that evnen this test is satistified visual inspection may still
# indicate some early non-stationarity which is to brief to effect the gelman.diag
# but which clearly should be removed using the methods shown above.


# The following code follows up the 100 participant case with small sample
# size per participant. Again convergence is fast and the fit excellent, and
# parameter recovery is only slightly worse. Running this is left as an exercise. 
hsamples.100 <- h.samples.dmc(nmc=100,p.prior,data.model1,thin=10,pp.prior=pp.prior)
hsamples.100 <- h.run.unstuck.dmc(hsamples.100,cores=4,
                                  p.migrate=0.05,h.p.migrate=0.05)
hsamples1.100 <- h.run.converge.dmc(h.samples.dmc(nmc=100,samples=hsamples.100),
                                    cores=4,nmc=50,save="tmp100")
h.gelman.diag.dmc(hsamples1.100)
summary.dmc(hsamples1.100,hyper=TRUE)
pp.100 <- h.post.predict.dmc(hsamples1.100)
plot.pp.dmc(pp.100)


# Note that the h.run.converge.dmc argument save="tmp100" saves a file 
# tmp100.RData containing an object "samples" after each try corresponding to 
# the intermediate versions of the returned object (here hsamples1.100). This can 
# be useful because the runs can be long so you do not lose too much work in 
# case of a crash.



