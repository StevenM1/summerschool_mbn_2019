##################  DMC Lesson 3: Sampling


### Lesson 3.2:  Assessing a fit

rm(list=ls())

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("LNR","lnr.R")

load_data ("dmc_3_1.RData")
# load_data ("dmc_3_2.RData")

# Occasionally a chain can get "stuck", taking a long time to join the
# other chains in the posterior mode and having a much lower posterior
# likelihood. Migration usually fixes this and stuck chains can be easily 
# detected graphically, but the following function tries to 
# automatically pick bad chains based on the difference between the mean
# chain log-likelihood (m_i). In particular it takes the median of those means 
# (median = med(m_1,m_2, ...)) then calculates the differences (median - m_i) 
# classifying as bad all chains whose average fit is worse by a value greater 
# than argument cut (default = 10, but might need to be tuned ...)

pick.stuck.dmc(samples3,cut=10,verbose=TRUE)
# Deviation of mean chain log-likelihood from median of means
#    12     5     6     7     8     9    14    11    15     1    10     3 
#  0.66  0.65  0.54  0.28  0.23  0.23  0.15  0.00 -0.04 -0.10 -0.12 -0.14 
#     2    13     4 
# -0.23 -0.28 -0.33 
# Bad chains:numeric(0)

# We can check convergence formally using "R hat", Gelman's "multivariate
# proportional scale reduction factor". See ?gelman.diag for details. The dmc
# version turns off autoburnin (which discards the first half of the series)
# and sets transform=TRUE. gelman.diag works best with normally distributed
# chains, this option applies automatic transformations that attempt to achieve
# normality.
# The same result is produced by the convenience function
gelman.diag.dmc(samples3,split=FALSE)
# Potential scale reduction factors:
# 
#               Point est. Upper C.I.
# meanlog.true        1.02       1.03
# meanlog.false       1.03       1.05
# sdlog.true          1.02       1.03
# sdlog.false         1.03       1.06
# t0                  1.01       1.03
# 
# Multivariate psrf
# 
# 1.05

# Usually a result < 1.1 is taken to indicate convergence. This is a measure of
# the ratio of variability between and within chains, so value near 1 indicates
# that the chains are "mixing" (i.e., between chain variability is similar to
# within chain variability). 

# Another characteristic of convergence is that the chains are stable (not moving
# up or down). This can be checked with gelman.diag by splitting the chains into
# first and second halves and treating each half like a chain. If the chains are
# changing systematically that inflates between chain variability. This is the 
# default check applied by the dmc version of gelman.diag
gelman.diag.dmc(samples3)
# Potential scale reduction factors:
# 
#               Point est. Upper C.I.
# meanlog.true        1.04       1.06
# meanlog.false       1.04       1.06
# sdlog.true          1.03       1.05
# sdlog.false         1.05       1.07
# t0                  1.03       1.05
# 
# Multivariate psrf
# 
# 1.07

# Note that gelman.diag can be sensitive to series length, reporting lack of 
# convergence if the effective number of samples in chains are to small. The 
# small elevation here is probably due to this effect, as splitting halves the
# chain length. Note also that if your chain lengths are uneven one sample is
# dropped in so that all split chains have the same length.

# Because MCMC samples are autocorrelated the effective number of independent
# samples obtained is often much less than nominal. The coda effectiveSize 
# function can be used to provide an estimate. More is always better but a 
# minimum of 200 is sometimes recommended.
effectiveSize.dmc(samples3)
#  meanlog.true meanlog.false    sdlog.true   sdlog.false            t0 
#           435           476           446           442           537 

# Here we see that the nominal number of samples, 500*15 = 7500 is effectively
# only around 500 independent samples (i.e., less by a factor of about 15).

# You can look at the autocorrelation function of a chain to more directly
# see the redundancy. Autocorrelation is the correlation between the current
# value and the value L iterations ago, where L is called "lag". The plot shows
# these correlations as a function of lag.
acf.dmc(samples3,chain=1,par="meanlog.true")


# You can also look at the autocorrelation of the posterior log-likelihood. Here
# we look at all 15 chains
par(mfrow=c(3,5))
for (i in 1:samples3$n.chains) acf.dmc(samples3,chain=i)

# The default number of lags used by the acf function can be changed by the 
# specifying the max.lag argument. Plotting can be turned off and the calculated 
# values for 10 lags can be saved to an object as follows.
acfs <- acf.dmc(samples3,plot=FALSE,max.lag=10)

# Auto-correlation can be reduced by thinning by a factor of k (i.e., only using
# one in every k samples from each chain). If k is chosen to be around the lag where 
# autocorrelation approaches zero then the samples become close to independent.
# Note that thinning is not necessary, and always throws away some information,
# so some papers recommend it never be done. However, if k is chosen 
# appropriately there is little loss and the number of samples that need to be
# stored is reduced.

# Thinning can be applied directly during sampling. Here we thin at k=10, which
# seems a conservative value to remove most autocorrelation. Note that the real
# number of samples taken below is 200*10 = 2000, but only 200 are stored.

samples4 <- run.dmc(samples.dmc(nmc=200,samples=samples3,thin=10), 
                    cores=4,report=10)

# The reduced autocorrelation is evident in a comparison of the posterior
# log likelihood plot for the last 200 samples without thinning
plot.dmc(samples3,start=301)
# and the 200 samples with thinning (periods of no change are now rare)
plot.dmc(samples4)

# There is now little difference between the split and non-split gelman.diag
gelman.diag.dmc(samples4)
gelman.diag.dmc(samples4,split=FALSE)

# The effective number of samples is now around half of nominal (200*15=3000).
effectiveSize.dmc(samples4)
#  meanlog.true meanlog.false    sdlog.true   sdlog.false            t0 
#          1407          1631          1489          2297          1712 

# This is confirmed by autocorrelation plots.
par(mfrow=c(3,5))
for (i in 1:samples4$n.chains) acf.dmc(samples4,chain=i)


# Assess goodness of fit by first generating posterior predictive samples.
# "post.predict" by default randomly draws (without replacement) n.post posterior 
# parameter vectors and samples the same number of data points as in data
# from them. It stores quantiles (percentiles by default) from the sample as
# well as a smoothed estimate of the density. Argument probs sets the quantiles
# (by default c(1:99)/100) and the level of smoothing for the density is 
# controlled by argument "bw" (see ?density). If n.post=NA all posterior samples
# are used. If random=FALSE then n.post samples are drawn at regular intervals 
# (this option is good for minimizing autocorrelation).
# Running this function can be a little slow.

pp <- post.predict.dmc(samples4)

# There are two styles of plot, both showing the data in black and the 
# posterior predictive in red.

# The first (default), "pdf", plots defective probability density functions,
# with each panel having lines for each response
# defective means their integral = response probability
plot.pp.dmc(pp)

# The layout parameter to plot.pp.dmc sets how many panels appear on each plot. 
# By default  all of the levels of the first two factors are displayed with 
# extra cells on further plots. The position of panel labels (drawn with 
# "legend") can be altered with the pos argument (using the keyword method,
# see ?legend).  

# The second style, "cdf", plots defective cumulative density functions, where
# defective means the cdf asymptotes at response probability. This style of plot
# also includes an indication of uncertainty in the model, grey lines 
# corresponding to the predictions of each of the n.post posterior samples. 
plot.pp.dmc(pp,"cdf")

# By default circles are drawn at the 10th, 30th, 50th, 70th and 90th 
# percentiles for both data and the model. This can be changed as follows (here
# providing the semi-deciles).
plot.pp.dmc(pp,"cdf",percentiles=seq(5,95,5))

# You can choose to show only a subset of responses with the show.response 
# argument, where the response selected is indicated by an integer (in the same
# order as responses are listed in the legend. 
plot.pp.dmc(pp,"cdf",show.response=c(1))

# This is mostly useful when there are many responses, making the plot cluttered.
# For example if there were 6 responses show.response=c(1,3,5) would pick out 
# odd ones. When there are lots of responses the legend can also become unwieldy.
# You can turn off the part related to the model with model.legend=FALSE.

# You can also specify a ylim arguement. To control the x axis you can specify
# x.min.max, a 2-vector of minimum and maximum values (if the calculated xlim
# is within this range it will not be changed).

# NB1: Although these plots give a good global idea of fit, they can miss 
#      important finer details. The advanced lesson dmc_5_5_ggplot.R shows how
#      to use the very flexible ggplot2 package to plot observed and predicted
#      RT quantile, means and SDs as well as accuracy for follow-up checks. This
#      is particularly useful in complex patterns when you want to understand
#      patterns across design factors. This is complemented by the posterior 
#      predictive tests described below.

# The coda summary function shows that the sampling has done a good job of 
# parameter recovery 
summary.dmc(samples4)
# Iterations = 1:200
# Thinning interval = 1 
# Number of chains = 15 
# Sample size per chain = 200 
# 
# 1. Empirical mean and standard deviation for each variable,
#    plus standard error of the mean:
# 
#                   Mean       SD  Naive SE Time-series SE
# meanlog.true  -0.99106 0.008566 1.564e-04      2.372e-04
# meanlog.false -0.01837 0.012426 2.269e-04      3.210e-04
# sdlog.true     1.00558 0.008122 1.483e-04      2.187e-04
# sdlog.false    1.00433 0.010138 1.851e-04      2.509e-04
# t0             0.20029 0.001048 1.914e-05      2.642e-05
# 
# 2. Quantiles for each variable:
# 
#                   2.5%      25%      50%      75%   97.5%
# meanlog.true  -1.00878 -0.99643 -0.99100 -0.98549 -0.9743
# meanlog.false -0.04281 -0.02671 -0.01839 -0.01004  0.0059
# sdlog.true     0.99010  0.99996  1.00535  1.01089  1.0219
# sdlog.false    0.98452  0.99771  1.00428  1.01120  1.0239
# t0             0.19815  0.19958  0.20029  0.20101  0.2022

# In cases where you are doing a parameter recovery study the following function
# provides a quick check. Recall the true values are
p.vector  <- c(meanlog.true=-1,meanlog.false=0,
               sdlog.true=1,sdlog.false=1,t0=.2)

# The following function prints true values, 95% credible intervals and the bias
# in estimates (Median estimated minus true). Note estimation is said have
# achieved "coverage" if the true value falls in the 95% CI 95% of the time. 
# This indicates that the level of variability, as well as the central tendency
# has been well estimated. In this case it coverage was achieved in all
# cases (note that we had to set digits=3 as at the default 2 everything is
# 0.2 for t0)
check.recovery.dmc(samples4,p.vector,digits=3)
#                meanlog.true meanlog.false sdlog.true sdlog.false    t0
# True                 -1.000         0.000      1.000       1.000 0.200
# 2.5% Estimate        -1.009        -0.043      0.990       0.985 0.198
# 50% Estimate         -0.991        -0.018      1.005       1.004 0.200
# 97.5% Estimate       -0.974         0.006      1.022       1.024 0.202
# Median-True           0.009        -0.018      0.005       0.004 0.000

# This function can also be used to produce a nicely formatted summary in 
# fits to real data by simply omitting the true parameter vector.
check.recovery.dmc(samples4)
#                meanlog.true meanlog.false sdlog.true sdlog.false  t0
# 2.5% Estimate         -1.01         -0.04       0.99        0.98 0.2
# 50% Estimate          -0.99         -0.02       1.01        1.00 0.2
# 97.5% Estimate        -0.97          0.01       1.02        1.02 0.2

# It is often useful to view the estimates graphically and also impose a plot
# of the prior to make sure that the data "updated" the prior (i.e., to make 
# sure that you arent just getting back your prior assumptions). This is 
# particularly useful in real data where you cant check parameter recovery so 
# this is the only check you have available on whether the data are informative
# for your model.  
plot.dmc(samples4,p.prior=p.prior,rug.lty=0)
# Clearly the updating is very strong as the prior (red lines) is completely 
# dominated by the posterior. rug.lty=0 turns of CODAs default "rug"
# plot which shows posteior values as ticks on the x-axis.

# You can also view the degree of correlation among parameters using a pairs
# plot, for the LNR this correlation is fairly mild.
pairs.dmc(samples4)

# Sometimes you might want to check if the model reproduces a particular
# feature of your data. For example, we might construct a robust measure of
# skew from the data quartiles
round(quantile(samples4$data$RT,probs=c(.25,.5,.75)),2)
#  25%  50%  75% 
# 0.37 0.51 0.75 

# Positive skew is indicated by the 50% (median)m - 25% being smaller than the
# 75% to 50%. This can be calculated by the following function, where positive 
# values indicate positive skew.
qskew <- function(data) {
  out <- diff(diff(quantile(data$RT,probs=c(.25,.5,.75)))) 
  names(out) <- "qskew"
  out
}

round(qskew(samples$data),2)
# qskew 
# 0.1 

# ppp.dmc applies fun to posterior predictive data and calculates the probability
# that the observed value is greater than the posterior predictive values. 
# Extreme values of this "posterior predictive probability" (ppp) indicate misfit. 
# Optionally the posterior predictive distribution of fun can be plotted.
# By default 500 values are simulated and the ppp value returned. In this case it
# is close to 0.5 and the observed value is in the middle of the posterior 
# predictive distribution.
qskew.save <- ppp.dmc(samples4,fun=qskew)
# [1] 0.81

# NB1: The value may vary slightly with identical data due to sampling noise, as 
#      n.post gets bigger this reduces.
# NB2: This approach has been criticised as only being able to detect large problems
#      unless the data are very precise.

# By default, ppp.dmc prints the posterior predictive p value and invisibly 
# returns the vector of samples that make up the plot with the observed value 
# as an attribute.
attr(qskew.save,"observed")
#     qskew 
# 0.1041492 

# It can then be used to calculate things like credible intervals
round(quantile(qskew.save,probs=c(.025,.975)),3)
#  2.5% 97.5% 
# 0.093 0.108 

# Clearly the data value falls well within this interval.


### Automating the process: CAVEAT EMPTOR: this may not always work so well ...

# First lets set up a fresh set of samples.
samples.auto <- samples.dmc(nmc=100,p.prior,data=samples$data)

# First repeatedly get a new set of nmc samples until there are no stuck chains. 
# When verbose=TRUE the output of pick.stuck.dmc is reported. 

# In this case it takes three cycles to finish. Note this is not guaranteed to 
# work and will give up after max.try = 100 cycles by default.
samples1.auto <- run.unstuck.dmc(samples.auto,p.migrate=.05,verbose=TRUE)
# 10  20  30  40  50  60  70  80  90  100  
# Deviation of mean chain log-likelihood from median of means
#        12         9         8        14         3        13         5 
#  25946.02  19067.43  12910.50   8052.03   7779.07   7399.90   4312.62 
#        11        10         7         1         4         6         2 
#      0.00    -97.41   -855.46  -2463.80  -8493.15 -10982.28 -11520.02 
#        15 
# -11783.65 
# Bad chains: 12 9 8 14 3 13 5
# ...
# 10  20  30  40  50  60  70  80  90  100  
# Deviation of mean chain log-likelihood from median of means
#     9    13    11    10     3    15    12     1     5     4     6     8 
#  1.42  0.57  0.40  0.27  0.20  0.12  0.08  0.00 -0.31 -0.32 -0.77 -0.78 
#    14     2     7 
# -0.85 -0.91 -0.91 
# Bad chains: None

# Note: Sometimes migration can cause the posterior log-likelihoods to be higher
# than they should be, so when you start running without migration you will see
# them decrease. When the decrease is not very large the following auto 
# convergence procedure can fail to discard these samples. Also sometimes 
# there may still be a few iterations at the start of the final run with 
# migration on where a chain is still being pulled in, but this is not detected
# by pick.stuck.dmc. In such cases you may want to use run.unstuck.dmc with 
# argument end.no.migarte=TRUE, which casues a final set of samples to be taken
# with no migration.

# The following procedure adds samples to its input, and either keeps the longer series or
# or discards the initial nmc, depending on which has the best gelman.diag
# mpsrf. Best to set nmc fairly short so it can grow the series to an 
# appropriate length then remove earlier non-stationary sections. The procedure
# runs until gelman.diag < cut (1.1 by default). It is not guaranteed to work  
# and will give up after max.try = 100 cycles by default. Should be run
# after the stuck chains are dealt with as migration is assumed to be off. 

# We use samples.dmc to start afresh from the final output of the run.unstuck
# using quite a short series length. We use the default 1.1 cut (note that by
# default split=TRUE to check for stationarity). In this case no initial
# samples were discarded (index always increases below) and length goes up to 
# 350 before convergence. Note that we can specify the size of the increments
# either from the size of the input object (as here) or using the nmc= argument
# to run.converge. 
samples2.auto <- run.converge.dmc(samples.dmc(samples=samples1.auto,nmc=50),
                                  cut=1.1,verbose=TRUE,nmc=50)
# 10  20  30  40  50  
# MPSRF:  1.78453977093903 
# 60  70  80  90  100  
# [1] "N = 100 Multivariate psrf achieved = 1.39081576180761"
# 110  120  130  140  150  
# ...
# 360  370  380  390  400  
# [1] "Final multivariate psrf = 1.07814946020417"
# Effective sample size
#  meanlog.true meanlog.false    sdlog.true   sdlog.false            t0 
#           345           370           362           355           389 

# You can also ask for a minimum effective sample size. So far we have a little
# less than 400. You can request to exceed a value of either the mean or minimum
# over parameters. Here we ask for the minimum to exceed 500. We also specify
# a smaller value of nmc, otherwise increments would be in the size of the input
# object (400).
samples3.auto <- run.converge.dmc(samples2.auto,minN=500,nmc=50,verbose=TRUE)
# For the mean use e.g.,run.converge.dmc(samples2.auto,meanN=500)

# save_data (samples4,pp,samples.auto,
#            samples1.auto,samples2.auto,samples3.auto,file="dmc_3_2.RData")
