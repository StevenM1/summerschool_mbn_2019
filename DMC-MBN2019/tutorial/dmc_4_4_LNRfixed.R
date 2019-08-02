##################  DMC Lesson 4: Multiple Subjects


### Lesson 4.4:  Multiple subjects as fixed effects, LNR Model

# Note that methods to test differences between parameters for experiments with
# multiple subjects (applicable to both the fixed effects estimation described
# here and the hierarchical estimation described in the next lesson, with some 
# caveats) are given in the advanced lesson dmc_5_6_ParTests.R. Methods of 
# testing correlations between subject parameters and covariates are given in 
# the advanced lesson dmc_5_4_Plausible.R.

# Same setup as lesson 4.2, 40 subjects

rm(list=ls())

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("LNR","lnrPP.R")

# load_data ("dmc_4_4.RData")

# Data model
model <- model.dmc(p.map=list(meanlog="M",sdlog="M",t0="1",st0="1"),
                   match.map=list(M=list(s1=1,s2=2)),constants=c(st0=0),
                   factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="lnr")

# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "meanlog.true"  "meanlog.false" "sdlog.true"    "sdlog.false"   "t0"           

# Population distribution (used to simulate data not as a prior)
p.mu  <- c(meanlog.true=-1,meanlog.false=0,        # natural scale
           sdlog.true=log(1),sdlog.false=log(1),t0=log(.2)) # log scale

# Fairly tight distributions
p.sigma <- c(.2,.2,.2,.2,.1); names(p.sigma) <- names(p.mu)
p.prior <- prior.p.dmc(p1=p.mu,p2=p.sigma,
                       untrans=c(sdlog.true="exp",sdlog.false="exp",t0="exp"))

p.priorCP <- p.prior[-c(1:2)]
# plot population distributions
par(mfcol=c(2,3)); for (i in names(p.prior)) plot.prior(i,p.prior)

# 2e4 data points in total
raw.data <- h.simulate.dmc(model,p.prior=p.prior,n=250,ns=40)

# parameters
ps <- attr(raw.data,"parameters")
ps[,3:5] <- exp(ps[,3:5])
round(ps,2)
data.model <- data.model.dmc(raw.data,model)


# h.samples.dmc produces a list with one entry for each subject
samples <- h.samples.dmc(nmc=200,p.prior,data.model)

# Note that this is slower than the equivalent single subject example  
# even though the number of data points (2e4) is the same. By default for
# one core the nmc is reported after each subject is done.
system.time({
  samples <- h.run.dmc(samples,p.migrate=.05,cores=1)
})
#    user  system elapsed 
# 431.136   2.998 437.295 

# As with the single subject case this can be run multi-core, but now it is
# subjects rather than chains that are spread across cores. 
# N.B. There is no reporting of progress, but parallelization is quite efficient
#      here as it is per subject.
system.time({
  samples <- h.run.dmc(samples,p.migrate=.05,cores=4)
})
#    user  system elapsed 
#   5.348   1.190 163.734 

# You can specify which subject to plot (by default the first) 
plot.dmc(samples,subject=10,layout=c(2,3))
# NB: Explicitly pick a subject in e.g., theta.as.mcmc.list(samples[[10]])  

# This shows that most subjects are close to converged by 100 interations. Also,
# no chains are stuck.
lapply(samples,pick.stuck.dmc,start=100)

# This can also be checked with the following convenience function.
h.pick.stuck.dmc(samples,start=100)

# This process can be automated as follows to get a full set of samples that
# has not stuck chains.
samples <- h.run.unstuck.dmc(samples,p.migrate=.05,cores=4)

# You can conveniently look at many subjects at once to check for stuck chains
# as follows (not the assignment here just stops the return of an empty list
# from the lapply statement).
par(mfrow=c(4,5))
tmp=lapply(samples,plot.dmc,pll.chain=TRUE)


# Note that h.run.unstuck.dmc has the same final.no.migration argument 
# run.unstuck.dmc that can be used to get a final sample with no migration.

### Now, turn migration and get a bigger sample. 

# Because samples objects can get large when dealing with lots of subjects it 
# can be useful to "thin on the fly" (e.g., keep only one in every 10 
# iterations), which (given autocorrelation in chains) does not sacrifice any
# information. 

# Look at how much the effectiveSize differs from the actual number of samples 
# (200 x 15 = 3000) in the best (least correlated) case to get a rough idea of
# how much thinning to do.
es <- effectiveSize.dmc(samples)
max(unlist(lapply(es,max)))
# [1] 522

# As this around 5 times less than the actual number of samples we use thin=5.

# Note that nmc refers to the number of samples AFTER thinning, so the following
# performs 500 (nmc x thin) iternations.
samples1 <- h.run.dmc(h.samples.dmc(nmc=100,samples=samples,thin=5),cores=4)

# Note that the thin argument is also available for single subject analyses

# Check nothing has become stuck
h.pick.stuck.dmc(samples1)

# Check convergence using a convenience function that applies gelman.diag to
# each subject and returns subjects sorted from lowest to highest (it also works
# for a single subject samples object). The defaults are autoburnin=FALSE and 
# transform=TRUE, and split=TRUE. When applied to a list of subjects it returns 
# the multivariate R hat for each.
gelman.diag.dmc(samples1)
#   10   13   20    3   26   34   35   39    6   21    1    7   19   29   37   14 
# 1.05 1.05 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.07 1.07 1.07 1.07 1.07 1.07 
#   22   15    5    2   16   17   36   30   38    9   24   40   31   18   32   23 
# 1.07 1.07 1.07 1.07 1.07 1.07 1.07 1.07 1.07 1.07 1.07 1.08 1.08 1.08 1.08 1.08 
#   27   11    4   33    8   25   28   12 
# 1.08 1.08 1.08 1.09 1.09 1.10 1.10 1.10 
# Mean
# [1] 1.07

# Note that it invisibly returns the list of gelman.diag objects for each
# subject, e.g.
rhats <- gelman.diag.dmc(samples1)
rhats[[1]]
# Potential scale reduction factors:
# 
#               Point est. Upper C.I.
# meanlog.true        1.03       1.05
# meanlog.false       1.04       1.06
# sdlog.true          1.05       1.08
# sdlog.false         1.04       1.07
# t0                  1.03       1.06
# 
# Multivariate psrf
# 
# 1.07

# Suppose we wanted to get all subjects under a cutoff of Rhat = 1.05. Analogous 
# to the function for a single subject, run.converge.dmc, we can apply the 
# following function to try to achieve this aim for every subject. Note that as
# in the fitting there is no reporting of progress with multicore. 

samples2 <- h.run.converge.dmc(samples1,nmc=50,cut=1.05,cores=4)
gelman.diag.dmc(samples2)
#   21   26   11   39   36   22    9   29   25   15   35   13   32   23   37    3 
# 1.03 1.03 1.03 1.03 1.03 1.03 1.03 1.03 1.03 1.03 1.03 1.04 1.04 1.04 1.04 1.04 
#   30    2   17   33   19    1    7   12   24   40   28    8    5   10   14    4 
# 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.05 
#   20    6   27   18   34   31   38   16 
# 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 
# Mean
# [1] 1.04


# effectiveSize.dmc applies the coda effectiveSize function to each subject . 
es <- effectiveSize.dmc(samples2)

# We can assess the average and worst case as follows
round(sort(unlist(lapply(es,mean))))
#   6   8  19  38  20  13  30  14  15  26  37   4  31  18  34  27   5   2  40  35 
# 597 619 629 631 644 654 659 664 667 675 677 677 679 681 689 689 689 702 704 713 
#  10  17  22  12  11   1  28   9  25  16   7  29  21   3  24  39  33  23  36  32 
# 717 722 748 762 771 774 791 834 834 860 864 864 870 877 896 898 900 906 915 976 
round(sort(unlist(lapply(es,min))))
#  19   6  38  13  14  20   8  10  26  18  34  31  12  27  30  35  15   4   2  40 
# 531 557 576 585 587 589 592 593 598 600 601 604 607 607 608 615 616 620 621 621 
#   5  37  22  25  17  28   1  11  21  29  39   9  24   7  16   3  36  33  23  32 
# 648 650 657 666 678 712 716 725 730 766 807 814 826 834 839 847 858 859 864 918 

# Suppose we wanted a mean over parameters of at least 700 effective samples. 
samples3 <- h.run.converge.dmc(samples2,nmc=50,meanN=700,cores=4)

min(unlist(lapply(effectiveSize.dmc(samples3),mean)))
# [1] 796.2

# summary.dmc applies the coda summary function to each subject (it also 
# works in the single subject case). In the multiple-subject case it returns
# a list of coda summaries invisibly and prints a tabulation of the means. The
# mean results are close to the true values.
summary.dmc(samples3)
#      meanlog.true meanlog.false sdlog.true sdlog.false    t0
# 1           -0.77          0.05       0.19        0.28 -1.56
# 2           -1.03         -0.07      -0.17        0.10 -1.54
# ...
# 40          -1.00          0.01       0.21        0.22 -1.65
# Mean        -1.01         -0.05       0.04        0.01 -1.59

# The small samples size per particiapnt means that updating of the priors is
# not very strong. For subject 1, for example:
plot.dmc(samples3[[1]],p.prior=p.prior)

# h.post.predict.dmc calculates posterior predictives for each subject. 
pp <- h.post.predict.dmc(samples3)

# Note, this can be run in a parallel mode, with subjects run on separate cores, 
# e.g.,pp <- h.post.predict.dmc(samples3,cores=4), in which case there is no 
# progress indication.

# Goodness of fit can be assessed subject-by-subject
plot.pp.dmc(pp[[1]])
plot.pp.dmc(pp[[1]],style="cdf")

# An average over subjects can also be done for the cdf (but not the pdf)
# This is done automatically if just the list of subject posterior predictives
# is supplied as an argument. It is calculated by averaging data from all
# subjects, so weights subjects with more data more heavily. The same is true
# for the grey lines representing the predictions for each posterior sample. 
plot.pp.dmc(pp)

# We can look at particular aspects of the posterior predictions on AVERAGE as
# follows (following the example from lesson 3.2 )

qskew <- function(data) {
  out <- diff(diff(quantile(data$RT,probs=c(.25,.5,.75)))) 
  names(out) <- "qskew"
  out
}

qskew.save <- h.ppp.dmc(samples3,fun=qskew)

# In this case the observed value flls just outside the 95% CI
round(quantile(qskew.save,probs=c(.025,.975)),3)
#  2.5% 97.5% 
# 0.096 0.113 
attr(qskew.save,"observed")
# [1] 0.09303379

# Model selection can be performed by extracting a fixed effect group WAIC. 
# There are two ways to do this. One concatenates trial from each subject to 
# create interactions by trials pointwise log-likelihood matrix (i.e., points
# are trials). 

# First extract and concatenate. This can create a very big matrix so some  
# thinning is recommended, here without thinning the size is 3e7, so thinning 
# of at least 10 is needed. 
group_trial_ll <- group_trial_log_likes(samples3,thin_pointwise=10,cores=4)
# Warning message:
# In group_trial_log_likes(samples3, thin_pointwise = 10, cores = 4) :
#   Subjects do not all have the same number of interations, using first 20 for all.

# Here a warning is given because h.run.converge.dmc resulted in different
# numbers of interactions for some participants. The same number are used for 
# each subject in maintain equal weighting for each.

# The usual waic function can then be used to calculate the group waic
waic.dmc(group_trial_ll)
#            p         se_p         waic      se_waic 
#   167.671840     3.218172 23500.668086   354.402369 

# The second method sums log-likelihoods over trials and operates on an
# iterations by subjects pointwise log-likelihods matrix (i.e., points are 
# subjects). It only make sense when the number of subjects is fairly large.
# No thinning or check for size is done as that is unlikely to be an issue.
group_subject_ll <- group_subject_log_likes(samples3)

# Again the usual waic function can then be used to calculate the group waic
waic.dmc(group_subject_ll)
#            p         se_p         waic      se_waic 
#   102.011055     2.641135 23468.982077  1429.508786 

# Auto-convergence can result in different lenght series for each participant.
unlist(lapply(samples3,function(x){x$nmc}))
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 200 200 250 200 200 200 250 200 250 200 250 200 200 200 200 250 200 200 200 200 
#  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40 
# 250 200 250 250 250 200 200 250 250 200 200 250 250 200 200 250 200 200 250 200 

# Suppose you want to make them all the same length. To do this specify remove
# in h.samples.dmc. Note that cases where the range is outside the number of 
# samples will be ignored (you can also make the upper end as long as you want,
# here it is larger than any series, 201:250 will do the same thing)
tmp <- h.samples.dmc(samples=samples3,add=TRUE,nmc=0,remove=201:300)
unlist(lapply(tmp,function(x){x$nmc}))

#### Extensions 1: Replacing bad chains

# Just as for a samples.dmc for a single subject, h.samples.dmc has a
# replace.bad.chains argument. The argument must either be a list of the same 
# length as samples whose entries specify which chains to replace or TRUE in 
# which case it chooses which chains to remove by the same heuristics as the 
# single subject case: 
# 1) get mean deviation of each chain from median and identify any chains that 
#    have values more than cut*iqr away (i.e., chains consistently high or low).
# 2) Same for log(IQR) (targeting unchanging chains)
# When invoked prints out a list of what has been identified.

# For example, with the intermediate samples1 object
samples1.cut2 <- h.samples.dmc(nmc=100,samples=samples1,replace.bad.chains = TRUE)
# Replacing bad chains
# $`6`
# [1] 13
# 
# $`31`
# [1] 10

# The chains chosen in this example are only weakly problematic and so probably
# not worth removing, but with more difficult models it can be useful when a 
# small subjet of subjects prove difficult to converge.

#### Extensions 2: Example with 100 subjects

# ### Create a 100 subject samples, again 2e4 data points in total
raw.data1 <- h.simulate.dmc(model,p.prior=p.prior,n=100,ns=100)
# parameters
ps1 <- attr(raw.data1,"parameters")
ps1[,3:5] <- exp(ps1[,3:5])
round(ps1,2)
apply(ps1,2,mean)
data.model1 <- data.model.dmc(raw.data1,model)


# h.samples.dmc produces a list with one entry for each subject
samples100 <- h.samples.dmc(nmc=200,p.prior,data.model1)
samples100 <- h.run.unstuck.dmc(samples100,p.migrate=.05,cores=4)
samples100.1 <- h.run.converge.dmc(
  h.samples.dmc(nmc=100,samples=samples100,thin=5),nmc=40,cores=4)
# Check convergence
gelman.diag.dmc(samples100.1)
#   39   61   17   74   33   49   36   98   83   96   32   85   12   93    3   42   97 
# 1.03 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.04 
#   53   11   86   64   20   15   95   59   22   75   25   90   47   27   13   80   56 
# 1.04 1.04 1.04 1.04 1.04 1.04 1.04 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 
#   34    6   41   60   52   73   71   68   62   30   57   19   76   79    7   35   70 
# 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 
#   77   69   82   44   38    9   66   58   14   43   54   23   16   46   55  100   40 
# 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.06 1.06 
#   91   78   26   24   18   67   65    1   37   81    2   45   99   51   29   87   50 
# 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 
#   28   63   89    5   84   48    8   21    4   92   94   31   88   72   10 
# 1.06 1.07 1.07 1.07 1.07 1.07 1.07 1.07 1.08 1.08 1.08 1.08 1.09 1.09 1.10 
# Mean
# [1] 1.05
summary.dmc(samples100.1)
#      meanlog.true meanlog.false sdlog.true sdlog.false    t0
# 1           -0.71         -0.22       0.09        0.12 -1.73
# 2           -1.18          0.17       0.13       -0.01 -1.57
# ...
# 100         -0.80          0.19       0.01        0.18 -1.55
# Mean        -1.02          0.00       0.03        0.00 -1.60

# !!! This produces a very large object so is not saved below !!!
pp.100 <- h.post.predict.dmc(samples100.1)
# The small sample nature of individual plots is evident in irregular cdfs
plot.pp.dmc(pp.100[[1]],style="cdf")

# But overall the fit looks good.
plot.pp.dmc(pp.100)


# save_data (model,p.mu,p.sigma,p.prior,ps,raw.data,raw.data1,data.model,data.model1,samples,
#            samples1,samples2,samples3,pp,samples100,samples100.1,file="dmc_4_4.RData")

### Extensions 3: Setting different constant priors for each subject.

# Different constant priors for each subject can be imposed by the following 
# procedure. Suppose we want to use the true values for meanlog.true and 
# meanlog.false as stord in ps with a fairly tight standard deviation (0.01).
# First make a list of models, one for each subject.

ps <- attr(raw.data,"parameters")
models <- vector(mode="list",length=dim(ps)[1])
for (i in 1:dim(ps)[1]) {
  constant.prior <- prior.p.dmc(
    dists = c("tnorm","tnorm"),
    p1=c(meanlog.true=ps[i,1],meanlog.false=ps[i,2]),                           
    p2=c(.01,.01),lower=c(NA,NA),upper=c(NA,NA)
  )
  models[[i]] <- model.dmc(p.map=list(meanlog="M",sdlog="M",t0="1",st0="1"),
                           match.map=list(M=list(s1=1,s2=2)),constants=c(st0=0),
                           constant.prior=constant.prior,
                           factors=list(S=c("s1","s2")),responses=c("r1","r2"),
                           type="lnr",verbose=FALSE)
}

data.modelCP <- data.model.dmc(raw.data,models)

# Looks to have set things as required
ps[,1]-unlist(lapply(data.modelCP,function(x){attr(attr(x,"model"),"constant.prior")[[1]]$mean}))
ps[,2]-unlist(lapply(data.modelCP,function(x){attr(attr(x,"model"),"constant.prior")[[2]]$mean}))

# Run a fit to check

p.priorCP <- p.prior[-c(1:2)]
samplesCP <- h.samples.dmc(nmc=100,p.priorCP,data.modelCP)
samplesCP <- h.run.unstuck.dmc(samplesCP,p.migrate=.05,cores=40)
samplesCP.1 <- h.run.converge.dmc(
  h.samples.dmc(nmc=200,samples=samplesCP,thin=5),nmc=40,cores=40)
gelman.diag.dmc(samplesCP.1) # All good

# Recall that this was good with all estiamted
ps <- attr(raw.data,"parameters")
h.check.recovery.dmc(samples3,ps)
#                meanlog.true meanlog.false sdlog.true sdlog.false     t0
# True                 -1.020        -0.046      0.037       0.003 -1.589
# 2.5% Estimate        -1.119        -0.191     -0.060      -0.113 -1.658
# 50% Estimate         -1.011        -0.049      0.039       0.007 -1.584
# 97.5% Estimate       -0.902         0.108      0.136       0.132 -1.538
# Median-True           0.009        -0.003      0.002       0.004  0.005

# As good if not better here!
h.check.recovery.dmc(samplesCP.1,ps[,-c(1:2)])
#                sdlog.true sdlog.false     t0
# True                0.037       0.003 -1.589
# 2.5% Estimate      -0.055      -0.096 -1.642
# 50% Estimate        0.037       0.004 -1.583
# 97.5% Estimate      0.133       0.105 -1.541
# Median-True         0.000       0.001  0.006
