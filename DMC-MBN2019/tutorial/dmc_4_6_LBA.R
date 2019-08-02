##################  DMC Lesson 4: Multiple Subjects

### Lesson 4.6: Hierarchical LBA Model 2 x 2 with a rate effect on factor F.


rm(list=ls())

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("LBA","lba_B.R")

# load_data ("dmc_4_6_fixed.RData")
# load_data ("dmc_4_6_random.RData")

# 2x2 (Stimulus and factor F) specified in model.dmc call
model <- model.dmc(p.map     = list(A="1",B="R",t0="1",mean_v=c("F","M"),sd_v="M",st0="1"), 
                   match.map = list(M=list(s1=1,s2=2)),
                   factors=list(S=c("s1","s2"),F=c("f1","f2")),
                   constants = c(sd_v.false=1,st0=0), 
                   responses = c("r1","r2"),
                   type="norm")
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"               "B.r1"            "B.r2"            "t0"             
# [5] "mean_v.f1.true"  "mean_v.f2.true"  "mean_v.f1.false" "mean_v.f2.false"
# [9] "sd_v.true"      
# 
# Constants are (see attr(,"constants") ):
# sd_v.false        st0 
#          1          0 


# Population distribution, rate effect on F
pop.mean <- c(A=.4, B.r1=.6, B.r2=.8, 
              mean_v.f1.true=1.5, mean_v.f2.true=1, mean_v.f1.false=0, mean_v.f2.false=0, 
              sd_v.true = .25,  t0=.3)
pop.scale <-c(A=.1, B.r1=.1, B.r2=.1, 
              mean_v.f1.true=.2, mean_v.f2.true=.2, mean_v.f1.false=.2, mean_v.f2.false=.2, 
              sd_v.true = .1, t0=.05)  
pop.prior <- prior.p.dmc(
  dists = rep("tnorm",9),
  p1=pop.mean,p2=pop.scale,
  lower=c(0,0,0,NA,NA,NA,NA,0,.1),upper=c(NA,NA,NA,NA,NA,NA,NA,NA,1)
)

##  Check population distributions
par(mfcol=c(2,5)); for (i in names(pop.prior)) plot.prior(i,pop.prior)


# Simulate some data
raw.data <- h.simulate.dmc(model, p.prior = pop.prior, n = 250, ns = 40)
data.model <- data.model.dmc(raw.data, model)

# Take a look at the first 2 subjects data
par(mfcol=c(2,4)) 
for (i in 1:2) { # Upper=biased to response, Lower = biased away. First column = greater rate, second lesser
  plot.cell.density(data.cell=data.model[[i]][data.model[[i]]$F=="f1" & data.model[[i]]$S=="s1",],C="r1")
  plot.cell.density(data.cell=data.model[[i]][data.model[[i]]$F=="f1" & data.model[[i]]$S=="s2",],C="r2")
  plot.cell.density(data.cell=data.model[[i]][data.model[[i]]$F=="f2" & data.model[[i]]$S=="s1",],C="r1")
  plot.cell.density(data.cell=data.model[[i]][data.model[[i]]$F=="f2" & data.model[[i]]$S=="s2",],C="r2")
}  

# Take a look at parameters
ps <- round( attr(raw.data, "parameters"), 2)
round(apply(ps,2,mean),2)
#               A            B.r1            B.r2  mean_v.f1.true  mean_v.f2.true 
#            0.40            0.57            0.83            1.53            1.03 
# mean_v.f1.false mean_v.f2.false       sd_v.true              t0 
#           -0.04            0.01            0.27            0.30 
round(apply(ps,2,sd),2)
#               A            B.r1            B.r2  mean_v.f1.true  mean_v.f2.true 
#            0.09            0.09            0.10            0.19            0.16 
# mean_v.f1.false mean_v.f2.false       sd_v.true              t0 
#            0.21            0.20            0.11            0.06 

### FIT FIXED EFFECTS

# specify a broader prior than the true population distribution
p.prior <- prior.p.dmc(
  dists = rep("tnorm",9),
  p1=pop.mean,                           
  p2=pop.scale*5,
  lower=c(0,0,0,NA,NA,NA,NA,0,.1),upper=c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
)
##  Check population distributions
par(mfcol=c(2,5)); for (i in names(p.prior)) plot.prior(i,p.prior)

# Start with thin of 5
samples  <- h.samples.dmc(nmc = 50, p.prior, data.model, thin = 5)
# Remember there is no output as this runs, except after each subject is done
# the number of iterations for that subject is printed; it will take a while!

# We will do some timing to be used later on when considering hierarchical fitting.
system.time({samples  <- h.run.unstuck.dmc(samples, p.migrate = .05, cores = 27)})
#    user  system elapsed 
#  65.575  14.186 157.742 

# Now run without migration, and turning up thinning to 10
system.time({samples1 <- h.run.converge.dmc(h.samples.dmc(nmc=100, samples=samples), 
                                            nmc=50,cores=27)})
#     user   system  elapsed 
#  133.917   35.552 4046.275 


# All chains are converged. 
gelman.diag.dmc(samples1)
#   15   24   21   13   22   12    3   33   34    7    6   37    5    2   36 
# 1.04 1.04 1.04 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 1.05 
#   11   39   16   35   18   32   14   30   25   27   28   38   31   10   29 
# 1.05 1.05 1.05 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 
#   23   20    8   17   19   26    9   40    1    4 
# 1.06 1.06 1.06 1.06 1.06 1.07 1.07 1.07 1.07 1.07 
# Mean
# [1] 1.06

# Parameter chains look generally well converged and mixed, but with the
# occasional wayward chain.
plot.dmc(samples1[[1]],pll.chain=TRUE)
for (i in 1:1) {
  plot.dmc(samples1,density=FALSE,smooth=FALSE,subject=i,layout=c(2,5))
  #  Sys.sleep(0.1) 
}

# Effective sample size OK but clearly could have been thinned further.
es <- effectiveSize.dmc(samples1)
round(apply(data.frame(es),1,mean))        
#               A            B.r1            B.r2              t0  mean_v.f1.true 
#             982             825             824             882             875 
#  mean_v.f2.true mean_v.f1.false mean_v.f2.false       sd_v.true 
#             862             906             894             926 
round(apply(data.frame(es),1,min))
#               A            B.r1            B.r2              t0  mean_v.f1.true 
#             576             548             552             575             574 
#  mean_v.f2.true mean_v.f1.false mean_v.f2.false       sd_v.true 
#             588             597             593             597 

# We can use get.thin to look at mean and minimum effective sample sizes and 
# an estimate of what thinning might be needed to remove all autocorrelation
# as actual n divided by either minimum or mean effective n. 
get.thin(samples1)
# Minimum Effective Size
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 1301  755  609  961  894  548  729  648  739  732  710  839  574  617  616 
#   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30 
#  675 1197  719 1380  633 1581  569  788 1097  645  913  701  677 1142  580 
#   31   32   33   34   35   36   37   38   39   40 
#  639  640 1244 1055  909  857  855  618  569  796 
# Mean Effective Size
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 1425  815  665 1028  991  597  770  688  805  756  753  881  614  686  674 
#   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30 
#  754 1320  780 1497  698 1798  616  855 1242  687  955  761  696 1253  602 
#   31   32   33   34   35   36   37   38   39   40 
#  701  680 1311 1174  971  947  903  653  596  857 
# Thin
#       1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
# mean 87 50 43 56 50 45 47 39 54 43 43 43 37 47 44 50 51 52 50 43 68 37 44
# min  95 54 47 60 56 49 50 42 58 44 46 45 40 53 48 56 56 56 54 47 77 40 48
#      24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
# mean 76 43 48 67 43 50 40 56 54 41 64 58 46 43 39 32 36
# min  86 46 50 73 44 54 42 61 57 43 72 62 50 46 42 33 39

# If we were using these fits we might want more samples and greater thinning, 
# but here we mainly use as a preamble to hierarchical fitting so will stop here.

# Good fits on average
pp <- h.post.predict.dmc(samples1) 
plot.pp.dmc(pp)

# Here we use an lapply to get all fits, the assign to tmp is just to stop
# an empty list output. Most fits are pretty good, even subject 1. 
tmp <- lapply(pp, function(x){plot.pp.dmc(x, style="cdf") })

# Check parameter recovery with the following convenience function which takes
# averages over subjects.
h.check.recovery.dmc(samples1,ps)
#                    A  B.r1  B.r2 mean_v.f1.true mean_v.f2.true
# True           0.403 0.565 0.827          1.532          1.030
# 2.5% Estimate  0.240 0.414 0.630          0.194          1.273
# 50% Estimate   0.422 0.658 0.940          0.273          1.658
# 97.5% Estimate 0.577 0.990 1.336          0.341          2.040
# Median-True    0.019 0.093 0.114         -0.023          0.126
#                mean_v.f1.false mean_v.f2.false sd_v.true    t0
# True                    -0.038           0.009     0.266 0.296
# 2.5% Estimate            0.842          -0.563    -0.533 0.219
# 50% Estimate             1.130           0.159     0.168 0.283
# 97.5% Estimate           1.430           0.721     0.702 0.346
# Median-True              0.099           0.198     0.159 0.017

# If the result is assigned then a list of matrices for each subject is saved.

# save_data(p.prior,raw.data,samples,samples1,pp,est,file="dmc_4_6_fixed.RData")


### FIT RANDOM EFFECTS

# Use all truncated normal priors for locations
mu.prior <- prior.p.dmc(
  dists = rep("tnorm",9),
  p1=pop.mean,                           
  p2=c(1,1,1,2,2,2,2,1,1),
  lower=c(0,0,0,NA,NA,NA,NA,0,.1),upper=c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
)
par(mfcol=c(2,5)); for (i in names(mu.prior)) plot.prior(i,mu.prior)


sigma.prior <- prior.p.dmc(
  dists = rep("beta", length(p.prior)),
  p1=c(A=1, B.r1=1, B.r2=1, 
       mean_v.f1.true=1, mean_v.f2.true=1, mean_v.f1.false=1, mean_v.f2.false=1, 
       sd_v.true = 1, t0=1),p2=rep(1,9)
)
par(mfcol=c(2,5)); for (i in names(sigma.prior)) plot.prior(i,sigma.prior)

# Make a hyper-prior list
pp.prior <- list(mu.prior, sigma.prior)

# Make a new samples object from the 40 subject data, longer and with more
# immediate thinning that LNR given likely more autocorrelation
hsamples <- h.samples.dmc(nmc = 100, p.prior, data.model, thin = 5, pp.prior = pp.prior)

# Fit with migration at both levels, finishes in two iterations.
system.time({hsamples <- h.run.unstuck.dmc(hsamples, cores=27, report=10, 
                                           p.migrate = .05, h.p.migrate = .05)})
#     user   system  elapsed 
# 3292.390  738.729 8783.620 

# Finishes after 5 tries with 200 samples 
system.time({hsamples1 <- h.run.converge.dmc(h.samples.dmc(nmc=100, samples=hsamples), 
                                             nmc=50,cores=27)})
#      user    system   elapsed 
#  6898.445  1487.189 14638.732

# Subject 1 is no longer a problem
h.gelman.diag.dmc(hsamples1)
# hyper    39     2    40     3    16    17     8    34     4    12    25    37 
#  1.03  1.05  1.05  1.05  1.05  1.05  1.05  1.05  1.05  1.05  1.05  1.06  1.06 
#    21    18     5    22    31    36    19    27    13    29    14     9    32 
#  1.06  1.06  1.06  1.06  1.06  1.06  1.06  1.06  1.06  1.06  1.06  1.06  1.06 
#    38    10    11     1    28    26    35     6    33    20    15    23     7 
#  1.06  1.06  1.06  1.06  1.07  1.07  1.07  1.07  1.07  1.07  1.07  1.07  1.07 
#    30    24 
#  1.08  1.08 

# Now lets try using start values from the fixed effects run. For the hyper 
# make.hstart creates a prior object that can be used to sample start points
# from. Mean of mu prior is mean of subject means, mean of sigma prior is SD of
# subject means. By default sets prior sd for mu and sigma prior to 1% of mean 
# (arguments mu.sd=1, sigma.sd=1, set show.sd=TRUE to see values)
hstart <- make.hstart(samples1)
# Mu prior
# Mean
#               A            B.r1            B.r2              t0  mean_v.f1.true 
#            0.42            0.67            0.95            0.27            1.66 
#  mean_v.f2.true mean_v.f1.false mean_v.f2.false       sd_v.true 
#            1.13            0.14            0.15            0.28 
# Sigma prior
# Mean
#               A            B.r1            B.r2              t0  mean_v.f1.true 
#            0.12            0.14            0.13            0.06            0.21 
#  mean_v.f2.true mean_v.f1.false mean_v.f2.false       sd_v.true 
#            0.21            0.27            0.26            0.12

# theta1 get last sampled thetas from each subject           
theta1 <- make.theta1(samples1)

# Note that these functions also work with hierarchical models with a group 
# parameterization.

# The start points are given via the hstart.prior and theta1 arguments. 
hsamplesA <- h.samples.dmc(nmc=100,p.prior,data.model,thin=5, 
                           pp.prior = pp.prior, hstart.prior=hstart,theta1=theta1)

# No stuck chains on initial run
system.time({hsamplesA <- h.run.unstuck.dmc(hsamplesA, cores=27, report=10, 
                                            p.migrate = .05, h.p.migrate = .05)})
#     user   system  elapsed 
# 1694.934  355.858 3762.415 

# Finished after trys with 200 samples
system.time({hsamplesA1 <- h.run.converge.dmc(h.samples.dmc(nmc=100, samples=hsamplesA), 
                                              nmc=50,cores=27)})
#      user    system   elapsed 
#  5992.926  1213.604 12278.385

# Similar quality solution.
h.gelman.diag.dmc(hsamplesA1)
# hyper    29    21     6    31     9    33    35    30    40    22    16    37 
#  1.03  1.05  1.05  1.05  1.05  1.06  1.06  1.06  1.06  1.06  1.06  1.06  1.06 
#    38    34    28    18     1     4    17     3     2    23    15    19    13 
#  1.06  1.06  1.06  1.06  1.06  1.06  1.06  1.06  1.06  1.06  1.06  1.06  1.07 
#    24    14    39    20    36    12     7    11    32     5    10    26    27 
#  1.07  1.07  1.07  1.07  1.07  1.07  1.07  1.07  1.07  1.07  1.07  1.07  1.07 
#     8    25 
#  1.08  1.08 

# Total time for run from random start values: 
#  8784+14638=23422
# Total time with fixed effect starts (including time to get them):
#  158+4046+3762+12278=20244
#
# The speedup is in part due to the more efficient usage of cores in the fixed
# effect fits (as each participant takes a long time to fit) compared to the
# random effect fits (one core per chain, where each chain is quick so there
# is lots of idle time due to scatter and gather operations). 

# Sometime drawing start points from the prior can also fail. You can try a 
# tighter prior (argument h.start.prior) but if that fails using the fixed
# effect based start points is a good option. It also allows checking of the
# effects of shrinkage.

# Chains look very nice.
plot.dmc(hsamples1,hyper=TRUE,pll.chain=TRUE)
plot.dmc(hsamples1,hyper=TRUE,layout=c(2,9))

# Good yield at hyper given nominal 27*200=5400
effectiveSize.dmc(hsamples1,hyper=TRUE)
#              A.h1            B.r1.h1            B.r2.h1  mean_v.f1.true.h1 
#               2700               2368               2286               2338 
#  mean_v.f2.true.h1 mean_v.f1.false.h1 mean_v.f2.false.h1       sd_v.true.h1 
#               2414               2241               1997               2311 
#              t0.h1               A.h2            B.r1.h2            B.r2.h2 
#               2425               2156               2335               2331 
#  mean_v.f1.true.h2  mean_v.f2.true.h2 mean_v.f1.false.h2 mean_v.f2.false.h2 
#               2222               2401               2109               2198 
#       sd_v.true.h2              t0.h2 
#               2051               2315 

# Less so at individual subject level  
es <- effectiveSize.dmc(hsamples1)
round(apply(data.frame(es),1,mean))        
#               A            B.r1            B.r2              t0  mean_v.f1.true 
#             963             962             959             963             976 
#  mean_v.f2.true mean_v.f1.false mean_v.f2.false       sd_v.true 
#             973             985             971             977 
round(apply(data.frame(es),1,min))        
#               A            B.r1            B.r2              t0  mean_v.f1.true 
#             846             846             854             835             843 
#  mean_v.f2.true mean_v.f1.false mean_v.f2.false       sd_v.true 
#             837             860             839             831 

# Note that you can use get.thin to check effective sample size as well as well 
# as well as the level of autocorrelation, both at the hyper level
get.thin(hsamples1,hyper=TRUE)
#               A.h1            B.r1.h1            B.r2.h1  mean_v.f1.true.h1 
#               2700               2368               2286               2338 
#  mean_v.f2.true.h1 mean_v.f1.false.h1 mean_v.f2.false.h1       sd_v.true.h1 
#               2414               2241               1997               2311 
#              t0.h1               A.h2            B.r1.h2            B.r2.h2 
#               2425               2156               2335               2331 
#  mean_v.f1.true.h2  mean_v.f2.true.h2 mean_v.f1.false.h2 mean_v.f2.false.h2 
#               2222               2401               2109               2198 
#       sd_v.true.h2              t0.h2 
#               2051               2315 
# Thin
# mean  min 
#    2    3

# The same can be done at the individual level where the autocorrelation is
# about twice as big. Note this is MUCH less than the fixed effects fits as 
# they were still some way off being fully converged. In most cases where
# the fixed effects fits are better converged get.thin proivdes a better guide
# for what might be used in hierarchical sampling, but it is probably best to
# always use a lower value of thinning just in case.
get.thin(hsamples1)
# Minimum Effective Size
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
#  892  900 1018  911  900  846  872  937  831  872  938  973  892  843  942 
#   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30 
#  885  886  908  906  896  907  871  875  900  835  959  837  858  936  866 
#   31   32   33   34   35   36   37   38   39   40 
#  857  839  904  921  946  908  932  880  846  919 
# Mean Effective Size
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
#  954  968 1067  957  962  937  944  985  894  961 1062  995  931  928  979 
#   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30 
#  957 1029  971  979  940  990  981  965  933  925 1010  882  961 1008  968 
#   31   32   33   34   35   36   37   38   39   40 
#  950  916 1003  979 1007  959  978 1026  985  970 
# Thin
#      1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
# mean 6 6 5 6 6 6 6 5 6  6  5  5  6  6  6  6  5  6  6  6  5  6  6  6  6  5
# min  6 6 5 6 6 6 6 6 6  6  6  6  6  6  6  6  6  6  6  6  6  6  6  6  6  6
#      27 28 29 30 31 32 33 34 35 36 37 38 39 40
# mean  6  6  5  6  6  6  5  6  5  6  6  5  5  6
# min   6  6  6  6  6  6  6  6  6  6  6  6  6  6


# Fits are nice.
hpp <- h.post.predict.dmc(hsamples1)
plot.pp.dmc(hpp)
tmp <- lapply(hpp, function(x){plot.pp.dmc(x, style="cdf") })

# Good updating, as usual scale less so than location.
plot.dmc(hsamples1,hyper=TRUE,p.prior=pp.prior,layout=c(2,9))


# Parameter recovery at the hyper level is conveniently checked with the 
# following function. Note that par.type=1 will get the mean (h1) parameters
h.check.recovery.dmc(hsamples1,pop.mean,ptype=1,hyper=TRUE)
#                    A  B.r1  B.r2 mean_v.f1.true mean_v.f2.true
# True           0.400 0.600 0.800          1.500          1.000
# 2.5% Estimate  0.382 0.588 0.866          1.551          1.039
# 50% Estimate   0.419 0.627 0.902          1.613          1.092
# 97.5% Estimate 0.457 0.666 0.939          1.680          1.147
# Median-True    0.019 0.027 0.102          0.113          0.092
#                mean_v.f1.false mean_v.f2.false sd_v.true     t0
# True                     0.000           0.000     0.250  0.300
# 2.5% Estimate            0.029           0.037     0.219  0.257
# 50% Estimate             0.103           0.112     0.271  0.278
# 97.5% Estimate           0.177           0.186     0.312  0.298
# Median-True              0.103           0.112     0.021 -0.022

# Theses values are uniformly closer to true values than fixed effect estimates.

# For scale estimates we use ptype=2
h.check.recovery.dmc(hsamples1,pop.scale,ptype=2,hyper=TRUE)
#                    A  B.r1  B.r2 mean_v.f1.true mean_v.f2.true
# True           0.100 0.100 0.100          0.200          0.200
# 2.5% Estimate  0.082 0.087 0.077          0.147          0.130
# 50% Estimate   0.106 0.110 0.100          0.185          0.163
# 97.5% Estimate 0.139 0.145 0.131          0.240          0.209
# Median-True    0.006 0.010 0.000         -0.015         -0.037
#                mean_v.f1.false mean_v.f2.false sd_v.true    t0
# True                     0.200           0.200     0.100 0.050
# 2.5% Estimate            0.146           0.150     0.097 0.048
# 50% Estimate             0.191           0.194     0.125 0.060
# 97.5% Estimate           0.254           0.259     0.173 0.081
# Median-True             -0.009          -0.006     0.025 0.010

# The improvement in scale estimates is even more marked! Note that if the 
# function is assigned the printed matrix is saved.

# This function can also be used to look at the average of individual parameter 
# recovery, which is typically similar-to the hyper.
h.check.recovery.dmc(hsamples1,ps)
#                    A  B.r1  B.r2 mean_v.f1.true mean_v.f2.true
# True           0.403 0.565 0.827          1.532          1.030
# 2.5% Estimate  0.324 0.535 0.794          1.483          0.998
# 50% Estimate   0.420 0.627 0.901          1.612          1.091
# 97.5% Estimate 0.511 0.726 1.016          1.745          1.189
# Median-True    0.017 0.061 0.074          0.079          0.061
#                mean_v.f1.false mean_v.f2.false sd_v.true     t0
# True                    -0.038           0.009     0.266  0.296
# 2.5% Estimate           -0.138          -0.114     0.243  0.237
# 50% Estimate             0.105           0.112     0.276  0.279
# 97.5% Estimate           0.338           0.326     0.311  0.316
# Median-True              0.143           0.103     0.010 -0.017

# In most cases this is better than for fixed effects.

# Note that the same function can be used to check hyper and average parameters
# in real data simply by omitting the second argument.

# save_data (raw.data,pp.prior,hsamples,hsamples1,hpp,hest,file="dmc_4_6_random.RData")

