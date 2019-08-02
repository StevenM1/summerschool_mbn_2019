rm(list=ls())
source ("dmc/dmc.R")

# # Load precomputed results of time-consuming computations 
load_data("DMCpaper2.RData")

### TAKE A LOOK AT THE DATA ----

# Score accuracy
tmp <- datgf$R; levels(tmp) <- c(NA,"left","right")
crct <- datgf$S==tmp

# Plot the data aggregated over subjects. We can see that:
# 1) Easy (~87.5%) more accurate than hard (~75%)
# 2) Very few go failures (around half a percent)
# 3) Stopping more successful for hard than easy
is.go <- datgf$SS=="go"; par(mfrow=c(2,4))
is.in <- datgf$D=="easy" & datgf$S=="left" & is.go; main="GO left easy"
plot.cell.density(datgf[is.in,c("R","RT")],C=crct[is.in],digits=3,main=main)
is.in <- datgf$D=="easy" & datgf$S=="right" & is.go; main="GO right easy"
plot.cell.density(datgf[is.in,c("R","RT")],C=crct[is.in],digits=3,main=main)
is.in <- datgf$D=="hard" & datgf$S=="left" & is.go; main="GO left hard"
plot.cell.density(datgf[is.in,c("R","RT")],C=crct[is.in],digits=3,main=main)
is.in <- datgf$D=="hard" & datgf$S=="right" & is.go; main="GO right hard"
plot.cell.density(datgf[is.in,c("R","RT")],C=crct[is.in],digits=3,main=main)
is.in <- datgf$D=="easy" & datgf$S=="left" & !is.go; main="STOP left easy"
plot.cell.density(datgf[is.in,c("R","RT")],C=crct[is.in],digits=3,main=main)
is.in <- datgf$D=="easy" & datgf$S=="right" & !is.go; main="STOP right easy"
plot.cell.density(datgf[is.in,c("R","RT")],C=crct[is.in],digits=3,main=main)
is.in <- datgf$D=="hard" & datgf$S=="left" & !is.go; main="STOP left hard"
plot.cell.density(datgf[is.in,c("R","RT")],C=crct[is.in],digits=3,main=main)
is.in <- datgf$D=="hard" & datgf$S=="right" & !is.go; main="STOP right hard"
plot.cell.density(datgf[is.in,c("R","RT")],C=crct[is.in],digits=3,main=main)

# Check race model prediction that failed stop ("signal-respond") RT is faster than
# go RT, breaking down by difficulty and response accuracy
round(tapply(datgf$RT,cbind.data.frame(D=datgf$D,CORRECT=crct,GO=is.go),mean,na.rm=TRUE),3)

# , , GO = FALSE
# 
#         CORRECT
# D      FALSE  TRUE
# hard 0.632 0.685
# easy 0.641 0.664
# 
# , , GO = TRUE
# 
#         CORRECT
# D      FALSE  TRUE
# hard 0.904 0.880
# easy 0.883 0.789

# Note that errors are slower than corrects for go RTs but the opposite is true
# for signal-respond RTs.

# Check at the individual level, taking GO-STOP differences collapsed over difficulty
# and accuracy with the expectation of positive differences according to the race model. 
# All differences are positive, supporting the applicability of the race model.
hist(apply(tapply(datgf$RT,cbind.data.frame(s=datgf$s,GO=is.go),mean,na.rm=TRUE),1,diff),breaks=seq(-.05,1.05,.05),xlab="GO-STOP RT",main="")

### DEFINE THE MODEL ----

# Note that attention failure probability parameters (gf=go failure and tf =
# trigger failure) are estimated on the probit scale. 
load_model ("EXG-SS","exgSSprobit.R")

# Define the experimental factors
factors <- list(S=c("left","right"),SS=c("go","stop"),D=c("hard","easy"))

# Make a model description
model <- model.dmc(p.map=list(mu=c("D","M"),sigma=c("D","M"),tau=c("D","M"),
                              muS="1",sigmaS="1",tauS="1",tf="1",gf="1"),
                   match.map=list(M=list(left="LEFT",right="RIGHT",left="NR")), 
                   factors=factors,responses=c("NR","LEFT","RIGHT"),
                   type="exgss")

#  [1] "mu.easy.true"     "mu.hard.true"     "mu.easy.false"    "mu.hard.false"   
#  [5] "sigma.easy.true"  "sigma.hard.true"  "sigma.easy.false" "sigma.hard.false"
#  [9] "tau.easy.true"    "tau.hard.true"    "tau.easy.false"   "tau.hard.false"  
# [13] "muS"              "sigmaS"           "tauS"             "tf"              
# [17] "gf"              
# 
# Constants are (see attr(,"constants") ):
# numeric(0)
# 
# Model type = exgss

# Bind the model to the data
data.model <- data.model.dmc(datgf,model)

### FITTING TO INDIVIDUALS ("FIXED EFFECTS") ----

# Construct a vague prior for individual subject fitting.

p1 <- p2 <- attributes(model)$p.vector # easy way to get a named vector
# Set mean of truncated normal priors
p1[1:length(p1)] <- c(rep(.5,4),rep(.2,8),.5,.1,.1,-1.5,-1.5)
# Set standard deviation of truncated normal priors
p2[1:length(p2)] <- rep(1,length(p2))

# Make priors
p.prior <- prior.p.dmc(
  dists = rep("tnorm",length(p1)),p1=p1,p2=p2,
  lower=c(rep(0,15),NA,NA),upper=c(rep(4,15),NA,NA)
)
# Plot the priors
par(mfcol=c(3,6)); for (i in names(p.prior)) plot.prior(i,p.prior)

# Make an object to contain samples with thinning of 5 (a guess)
samples <- h.samples.dmc(nmc=120,p.prior,data=data.model,thin=5)

# Run the sampling, VERY SLOW, run in a batch file, not in an R Session.
samples <- h.RUN.dmc(samples,cores=4)

gelman.diag.dmc(samples) # 1.05-1.1, Mean 1.08
get.thin(samples) # Thinning required at around 6-10 on top of 5

### HIERARCHICAL FITTING  ----

# Make start points based on individual fits
# First for the population-level parameters 
hstart <- make.hstart(samples)
# Mu prior
# Mean
#     mu.easy.true     mu.hard.true    mu.easy.false    mu.hard.false  sigma.easy.true 
#             0.57             0.61             0.94             0.83             0.08 
#  sigma.hard.true sigma.easy.false sigma.hard.false    tau.easy.true    tau.hard.true 
#             0.10             0.24             0.21             0.27             0.39 
#   tau.easy.false   tau.hard.false              muS           sigmaS             tauS 
#             1.04             0.83             0.24             0.07             0.08 
#               tf               gf 
#            -2.06            -2.81 
# Sigma prior
# Mean
#     mu.easy.true     mu.hard.true    mu.easy.false    mu.hard.false  sigma.easy.true 
#             0.18             0.22             0.54             0.50             0.06 
#  sigma.hard.true sigma.easy.false sigma.hard.false    tau.easy.true    tau.hard.true 
#             0.08             0.21             0.20             0.15             0.23 
#   tau.easy.false   tau.hard.false              muS           sigmaS             tauS 
#             0.51             0.48             0.06             0.05             0.10 
#               tf               gf 
#             0.65             0.50 

#Now for the participant-level parameters
theta1 <- make.theta1(samples)

# Make hyper-prior distributions

# Prior on the population means (i.e., locations); essentially identical to the 
# individual subject priors (line 107) with a small (and inconsequential) 
# tweak to the sigma and tau parameters
p.vector <- attributes(model)$p.vector
p.mu.mu <- p.vector
p.mu.mu[1:length(p.mu.mu)] <- c(rep(.5,4),rep(.1,8),.5,.1,.1,-1.5,-1.5)
p.mu.sigma <- p.mu.mu
p.mu.sigma[1:length(p.mu.sigma)] <- rep(1,length(p.vector))
mu.prior <- prior.p.dmc(
  p1=p.mu.mu,p2=p.mu.sigma,
  lower=c(rep(0,15),NA,NA),upper=c(rep(4,15),NA,NA))
# par(mfcol=c(3,6)); for (i in names(mu.prior)) plot.prior(i,mu.prior)

# Prior on the population standard deviations (i.e., scales); 
# all identical exponentials with a mean of 1
p.sigma.shape <- rep(1,17); names(p.sigma.shape) <- names(p.vector)
p.sigma.scale <- rep(1,17);names(p.sigma.scale) <- names(p.vector)
sigma.prior <- prior.p.dmc(p1=p.sigma.shape,p2=p.sigma.scale,
                           dists=rep("gamma",length(p.vector)))
# par(mfcol=c(3,6)); for (i in names(sigma.prior)) plot.prior(i,sigma.prior)

# This makes a figure 12a, which combines the priors for the population means (locations) 
# with one panel for the population standard deviation (scale) prior.
plotPrior <- prior.p.dmc(
  dists=c(rep("tnorm",length(p.vector)),"gamma"),
  p1=c(p.mu.mu,SCALE=1),p2=c(p.mu.sigma,1),
  lower=c(rep(0,15),NA,NA,NA),upper=c(rep(4,15),NA,NA,NA))
par(mfcol=c(3,6)); for (i in names(plotPrior)) plot.prior(i,plotPrior)

# Make a hyper-prior list
pp.prior=list(mu.prior,sigma.prior)

# Initialize samples
hsamples <- h.samples.dmc(nmc=100,p.prior,data.model,pp.prior,
                          hstart.prior=hstart,theta1=theta1)

# VERY SLOW, run in a batch file, not in an R Session
hsamples <- h.run.unstuck.dmc(hsamples,cores=4,report=1,p.migrate=0.05,h.p.migrate=0.05)
hsamples <- h.run.converge.dmc(samples=h.samples.dmc(samples=hsamples,nmc=120,thin=25),
                               nmc=40,cores=4,report=1,verbose=TRUE)
hsamples <- h.run.dmc(h.samples.dmc(samples=hsamples,nmc=500),report=1,cores=4)
# Chains don't look OK yet; get 300 new posterior samples
hsamples <- h.run.dmc(h.samples.dmc(samples=hsamples,nmc=300),report=1,cores=4)

# Mixing is good at population (hyper) level
gelman.diag.dmc(hsamples,hyper=TRUE) # Mean 1
# and individual subject level
gelman.diag.dmc(hsamples)
#  s05  c37  c09  s06  c15  c62  c49  c50  c19  c64  c31  c46  c39  c59 
# 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 
#  c48  c21  c44  c33  c18  c12  c07  c45  c41  s07  c66  c51  s08  c36 
# 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 
#  c65  c27  c11  c29  s11  c23  c03  c01  c60  c30  s01  c35  c22  c17 
# 1.01 1.01 1.01 1.01 1.02 1.02 1.02 1.02 1.02 1.02 1.03 1.03 1.04 1.05 
#  c52  c26  c57  c08  c63 
# 1.05 1.05 1.06 1.06 1.10 
# Mean
# [1] 1.02

# Posterior likelihoods flat with some slow waves despite heavy thinning
plot.dmc(hsamples,hyper=TRUE,pll.chain=TRUE)

# Chains for the population-level (hyper) parameters look nice
plot.dmc(hsamples,hyper=TRUE,layout=c(3,6))

# This makes the prior and posterior plots in Figure 12b and c
plot.dmc(hsamples,layout=c(3,6),p.prior=pp.prior,location=TRUE)
plot.dmc(hsamples,layout=c(3,6),p.prior=pp.prior,scale=TRUE)
# Group mean (location) and group standard deviation (scale) 
# for tau_false poorly updated

### STANDARD GOODNESS OF FIT CHECKS ----

# Create posterior predictives to check goodness of fit;
# This can be SLOW unless you have lots of cores available
pp <- h.post.predict.dmc(hsamples,cores=4, gglist=TRUE)

# We can now plot observed response proportions and the corresponding 
# 95% credible intervals of the posterior predictions 
# (i.e., the 2.5% and 97.5% quantiles of the predicted data).
# See lesson 5.5 for more detail on these functions. Because the proportion of
# non-responses on stop trials is of interest we must set the include.NR 
# argument to TRUE (by default proportions of trials with an RT response are
# plotted). The graph shows results at the group level, obtained by first
# generating predicted data for each subject and for each posterior sample, 
# then taking the average over subjects for each posterior sample. 

#set ggplot theme 
theme_set(theme_simple())
ggplot.RP.dmc(pp,include.NR=TRUE)  # Figure 13a

# Next plot the RT distributions, by default as the 10th, 50th and 90th 
# percentiles. Note that the NR response level is dropped because by definition
# it has no RT data.
ggplot.RT.dmc(pp) # Figure 13b

# Plot average cdf (Figure 13c)
plot.pp.dmc(pp,layout=c(2,4),ylim=c(0,1),fits.pcol="grey",
            model.legend=FALSE,dname="",aname="",x.min.max=c(.25,2))

# The following plot functions require saving posterior predictve simulations;
# Again this can be SLOW
pp.sim <- h.post.predict.dmc(hsamples,save.simulation = TRUE,cores=4)

### INDIVIDUAL INHIBITION AND SRRT FUNCTIONS ----

# DMC plots two functions commonly used in the stop-signal paradigm, the 
# inhibition function (IF: the probability of responding on stop-signal trials
# as a function of SSD) and the median signal-respond RT (SRRT; i.e., RTs on failed stop-signal trials)
# as a function of SSD. These can be plot with and without fits.

# As shown here, different subjects can have very different ranges of SSDs
par(mfcol=c(2,2))
plot_SS_if.dmc(data=hsamples[[1]]$data,main="IF: Subject 1")
plot_SS_if.dmc(data=hsamples[[3]]$data,main="IF: Subject 3")

# If posterior predictions are included (using the sim argument), DMC also reports
# (1) a table with the number of data points per SSD and
# the proportion of simulated data points (i.e., response rates) 
# higher than the observed data (p; i.e., posterior predictive p value);
# (2) "violin" plots for each SSD, these are density plots of the 
# predicted response rates reflected around the vertical axis.
plot_SS_if.dmc(data=hsamples[[1]]$data,sim=pp.sim[[1]],main="Subject 1")
# For Subject 3, the table is missing p values where the posterior predictions
# never produced a failed inhibition
plot_SS_if.dmc(data=hsamples[[3]]$data,sim=pp.sim[[3]],main="Subject 3")
# Subject 3 
#     SSD  n    p
# 1 0.100  1   NA
# 2 0.133  9 0.90
# 3 0.167 27 0.86
# 4 0.200 42 0.47
# 5 0.233 33 0.34
# 6 0.267 18 0.75
# 7 0.300  9 0.66
# 8 0.333  3 0.97
# 9 0.367  1   NA

# For SRRTs, only SSDs with at least one signal-respond RT are plotted.
par(mfcol=c(2,3))
plot_SS_srrt.dmc(data=hsamples[[1]]$data,main="Subject 1")
plot_SS_srrt.dmc(data=hsamples[[3]]$data,main="Subject 3")
# If posterior predictions are included (using the sim argument), DMC also reports
# (1) a table with the number of observed signal-respond RTs per SSD (nrt), 
# the posterior predictive p values, and the number of predicted data sets 
# with at least one signal-respond RT (n.sim);
# sometimes n.sim can be very low and so the p value is unreliable. 
# (2) "violin" plots for each SSD, these are density plots of the 
# predicted SRRTs reflected around the vertical axis.
# Below we show only the table for Subject 3:
plot_SS_srrt.dmc(data=hsamples[[1]]$data,sim=pp.sim[[1]],main="Subject 1")
plot_SS_srrt.dmc(data=hsamples[[3]]$data,sim=pp.sim[[3]],main="Subject 3")
# Subject 3 
#     SSD nrt     p n.sim
# 1 0.133   1 1.000    90
# 2 0.167   7 0.230   100
# 3 0.200  20 0.590   100
# 4 0.233  21 0.370   100
# 5 0.267  11 0.990   100
# 6 0.300   7 0.510   100
# 7 0.333   2 0.610   100
# 8 0.367   1 0.052    97

# Plots can also be made using only correct RTs. For subject 3, this causes the
# SSDs beyond .3 to be omitted as they contained only errors RTs. 
plot_SS_srrt.dmc(hsamples[[1]]$data,sim=pp.sim[[1]],do.correct=TRUE,main="Subject 1")
plot_SS_srrt.dmc(hsamples[[3]]$data,sim=pp.sim[[3]],do.correct=TRUE,main="Subject 3")
#     SSD nrt    p n.sim
# 1 0.133   1 1.00    90
# 2 0.167   7 0.23   100
# 3 0.200  15 0.44   100
# 4 0.233  19 0.37   100
# 5 0.267   6 0.43   100
# 6 0.300   6 0.39   100

### AVERAGE INHIBITION AND SRRT FUNCTIONS ----

# To make these plots, we have to average over ranges of SSDs or else the result 
# is very noisy. This causes an issue when deriving averages because, as 
# illustrated above, different subjects can have very different SSD ranges.
# DMC averages in two ways, either in terms of absolute time or by dividing each
# subject's SSDs up according to percentile ranges. The following examples show
# that for the inhibition function the percentile range method is best (it is
# the default) as the absolute method flattens the function in a way that is 
# not representative of individual subjects. For the SRRT function, 
# the opposite it true, so the absolute method is the default.
par(mfrow=c(2,2))
plot_SS_if.dmc(data=hsamples,n.intervals=6,violin=.25,percentile.av=FALSE)
plot_SS_if.dmc(data=hsamples,n.intervals=6,violin=.25)
plot_SS_srrt.dmc(data=hsamples,n.intervals=6,violin=.25)
plot_SS_srrt.dmc(data=hsamples,n.intervals=6,violin=.25,percentile.av=TRUE)

# We can explicitly specify ranges with the probs argument (NB. you don’t need to specify 
# 0 or 1) or just provide n.intervals so you get evenly spaced ranges (here 
# 20%). Note that the posterior predictive p values can be sensitive to how intervals are defined.
plot_SS_if.dmc(data=hsamples,sim=pp.sim,n.intervals=5,violin=.25)
# Average absolute SSD for each percentile range
#   (0,20%]  (20,40%]  (40,60%]  (60,80%] (80,100%] 
#     0.295     0.399     0.450     0.505     0.579 
# 
# Inhibition function 
#         SSD    n    p
# 1   (0,20%] 1363 0.22
# 2  (20,40%] 1356 0.05
# 3  (40,60%] 1322 0.13
# 4  (60,80%] 1356 0.96
# 5 (80,100%] 1363 0.96

# Note that you can use your own x-axis labels. For example say you wanted to use the
# average times above in ms (x-axis labels were added after an initial 
# run without them to find values)
plot_SS_if.dmc(data=hsamples,sim=pp.sim,n.intervals=5,violin=.25,
               xlabels=c(".295",".399",".450",".505",".579"),xlab="SSD (ms)")

# Here we use 12 percentile intervals. This is used to create Figure 14a. Note the 
# under-estimation at SSD = ~0.4s and over-estimation at SSD = ~.55s, confirming 
# the findings with the first set of intervals.
plot_SS_if.dmc(data=hsamples,sim=pp.sim,n.intervals=12,violin=.25,xlab="SSD")
# Average absolute SSD for each percentile range
#    (0,8%]   (8,17%]  (17,25%]  (25,33%]  (33,42%]  (42,50%]  (50,58%] 
#     0.237     0.328     0.367     0.397     0.419     0.441     0.463 
#  (58,67%]  (67,75%]  (75,83%]  (83,92%] (92,100%] 
#     0.483     0.508     0.534     0.564     0.610 
# 
# Inhibition function 
#     SSD   n    p
# 1  .237 564 0.64
# 2  .328 564 0.24
# 3  .367 564 0.37
# 4  .397 564 0.00
# 5  .419 563 0.06
# 6  .441 564 0.45
# 7  .563 558 0.78
# 8  .483 564 0.20
# 9  .508 563 0.42
# 10 .534 564 1.00
# 11 .564 564 0.99
# 12 .610 564 0.25

# There almost exactly equal numbers of observations in each interval because a 
# small random perturbation was added to each SSD to make it unique before 
# cutting them up, so that where the discrete set of SSDs breaks into categories 
# unevenly they are proportionally and randomly allocated. 
# NOTE:  this randomness can cause some small fluctuation in p values
# over repeated runs.

# For SRRT, the functionality is much as for the IF function, except that everything 
# is done on signal-respond trials only.

# Percentile averaging produces somewhat different results to  
# absolute averaging, a difference that is even more marked with many intervals, 
# where absolute averaging is much smoother.
par(mfrow=c(2,2))
plot_SS_srrt.dmc(data=hsamples,sim=pp.sim,n.intervals=6,violin=.25,percentile.av=TRUE)
plot_SS_srrt.dmc(data=hsamples,sim=pp.sim,n.intervals=6,violin=.25)
plot_SS_srrt.dmc(data=hsamples,sim=pp.sim,n.intervals=20,violin=.25,percentile.av=TRUE)
plot_SS_srrt.dmc(data=hsamples,sim=pp.sim,n.intervals=20,violin=.25)  

# Note the messy intervals due to random 
# perturbation of SSDs; for nice presentation, we can use xlabels to replace this.
# This is used to create Figure 14b
par(mfrow=c(1,1))
tmp <- plot_SS_srrt.dmc(data=hsamples,sim=pp.sim,n.intervals=20,
                        violin=.5,percentile.av=FALSE)
tmp <- round(attr(tmp,"cuts"),2)
xlabels <- paste("(",tmp[-length(tmp)],",",tmp[-1],"]",sep="")

plot_SS_srrt.dmc(data=hsamples,sim=pp.sim,n.intervals=20,
                 violin=.5,percentile.av=FALSE,xlabels=xlabels,xlab="SSD (s)")
# Median signal-respond RTs 
# SSD nrt    p n.sim
#1    (0,0.17]  160 0.13   100
#2   (0.17,0.2] 160 0.95   100
#3   (0.2,0.23] 159 0.50   100
#4  (0.23,0.27] 160 0.51   100
#5   (0.27,0.3] 159 0.77   100
#6   (0.3,0.33] 160 0.81   100
#7  (0.33,0.33] 159 0.66   100
#8  (0.33,0.37] 160 0.68   100
#9   (0.37,0.4] 159 0.71   100
#10  (0.4,0.43] 160 0.41   100
#11 (0.43,0.47] 160 0.24   100
#12  (0.47,0.5] 159 0.05   100
#13  (0.5,0.53] 160 0.34   100
#14 (0.53,0.57] 159 0.34   100
#15  (0.57,0.6] 160 0.29   100
#16  (0.6,0.63] 159 0.43   100
#17 (0.63,0.67] 160 0.00   100
#18 (0.67,0.73] 159 0.25   100
#19  (0.73,0.9] 160 0.15   100
#20  (0.9,1.33] 160 0.12   100

# Fit is generally good with only one interval (interval 17) with an extreme p-value 

# We can also plot correct RTs only:
par(mfrow=c(1,2))
plot_SS_srrt.dmc(data=hsamples,sim=pp.sim,n.intervals=6,violin=.25,
                 percentile.av=TRUE,do.correct=TRUE)
plot_SS_srrt.dmc(data=hsamples,sim=pp.sim,n.intervals=6,violin=.25,
                 do.correct=TRUE)

### POSTERIOR CORRELATIONS ----

# There is little correlation among the population-level mean (location) parameters
location.r <- pairs.dmc(hsamples,location=TRUE,do.plot=FALSE)
summary(location.r)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -0.0217354 -0.0057561  0.0009818  0.0006997  0.0062785  0.0297411

# There is little correlation among the population-level standard deviation (scale) parameters
scale.r <- pairs.dmc(hsamples,scale=TRUE,do.plot=FALSE)
summary(scale.r)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -0.0214038 -0.0058707 -0.0001101  0.0001704  0.0049704  0.0250429 

# The correlations among population-level (hyper) parameters are typically small,
# with a few exceptions:
hyper.r <- pairs.dmc(hsamples,hyper=TRUE,do.plot=FALSE)
hist(hyper.r,breaks="fd")

# The exceptions are the correlations between the 
# 12 ex-Gaussian location and scale parameters of the same type:
head(round(sort(hyper.r),2),12)

# To plot correlations at the subject level, we need to standardize each subject's
# parameters otherwise individual differences wash out the underlying 
# correlations. This is the default action of pairs.dmc when given a list of 
# subjects and when not focusing on the population-level parameters.
subject.r <- pairs.dmc(hsamples,do.plot=FALSE)

# Moderate negative (~ -.45 to -.75) tau-mu and tau-sigma correlations and 
# moderate positive mu-sigma correlations for go runners within each 
# of easy vs. hard x true vs. false.
# Weaker negative (~ -.65 and -.2) tauS~mus and tauS~sigmaS correlations for 
# stop runner, and ~ -.25 tau~tf correlation: 
head(round(sort(subject.r),2),11)
tail(round(sort(subject.r),2),4)

# These plots focus on the higher correlations
pairs.dmc(hsamples,thin=100,
          use=c("mu.hard.true","sigma.hard.true","tau.hard.true"))
pairs.dmc(hsamples,thin=100,
          use=c("mu.hard.false","sigma.hard.false","tau.hard.false"))
pairs.dmc(hsamples,thin=100,
          use=c("mu.easy.true","sigma.easy.true","tau.easy.true"))
pairs.dmc(hsamples,thin=100,
          use=c("mu.easy.false","sigma.easy.false","tau.easy.false"))
pairs.dmc(hsamples,thin=100,
          use=c("muS","sigmaS","tauS","tf"))

### POSTERIOR INFERENCE ----

# Comparing hard vs. easy ex-Gaussian parameters for the matching runner we 
# see that hard is always higher:
compare.p(hsamples,pnames=c("mu.hard.true","mu.easy.true"),pretty=c("hard","easy"))
#        hard  easy contrast
# 2.5%  0.593 0.562    0.024
# 50%   0.600 0.567    0.033
# 97.5% 0.608 0.572    0.042
# p.gt.0 
#      1

compare.p(hsamples,pnames=c("sigma.hard.true","sigma.easy.true"),pretty=c("hard","easy"))
#        hard  easy contrast
# 2.5%  0.086 0.069    0.012
# 50%   0.092 0.073    0.019
# 97.5% 0.097 0.076    0.026
# p.gt.0 
#      1 

compare.p(hsamples,pnames=c("tau.hard.true","tau.easy.true"),pretty=c("hard","easy"))
#        hard  easy contrast
# 2.5%  0.379 0.265    0.103
# 50%   0.393 0.273    0.120
# 97.5% 0.408 0.282    0.136
# p.gt.0 
#      1

# As discussed in tutorial 5.6, these tests account for the within-subject 
# correlations by testing the posterior distribution of the differences averaged 
# across the participants. Tests of the population-level (hyper) parameters 
# have great uncertainty because they do not take account of these correlations as shown
# in the following example where we test the population mean of the mu.true parameters
compare.p(hsamples,pnames=c("mu.hard.true","mu.easy.true"),
          pretty=c("hard","easy"),hyper=TRUE)
#        hard  easy contrast
# 2.5%  0.538 0.521   -0.045
# 50%   0.599 0.567    0.032
# 97.5% 0.658 0.614    0.107
# p.gt.0 
#  0.802 

# Proper tests of within-subjects effects at the population level would require a 
# differently structured hierarchical model (something we are working on). At
# present, tests of the population-level parameters should only be used for between-subject 
# comparisons, and inferences about the recommended within-subject tests 
# should be understood as applying to the sample of participants rather than
# to the population.

# Greater accuracy in easy seems to be mediated by much greater tau difference
# between match and mismatch. 
fun <- function(x){
  hard <- diff(x[c("tau.hard.true","tau.hard.false")])             
  easy <- diff(x[c("tau.easy.true","tau.easy.false")])             
  c(easy,hard,easy-hard)
}
compare.p(hsamples,show.plot=TRUE,pretty=c("easy","hard"),fun=fun)
#        easy  hard contrast
# 2.5%  0.892 0.610    0.223
# 50%   0.976 0.656    0.320
# 97.5% 1.061 0.704    0.414
# p.gt.0 
#      1

### PRIORS FROM POSTERIORS ----

# The following function takes the posterior samples for one participant and
# uses them to update the p.prior object with parameters corresponding to the
# best fitting truncated normal distributions
post.prior1 <- make.tnorm.prior(p.prior,hsamples[[1]])

# Fitting is done by optimization and is not always guaranteed to work, so it
# is best to plot the results to check. In this case the fits are quite good
plot.dmc(hsamples[[1]],p.prior=post.prior1,layout=c(3,6))

# The following convenience function allows you to extract the numeric values 
# from the p.prior object as a named matrix

ppars <- get.ppars(post.prior1)
round(ppars,3)
#                    mean    sd
# mu.hard.true      0.565 0.038
# mu.easy.true      0.456 0.013
# mu.hard.false     0.746 0.165
# mu.easy.false     0.989 0.168
# sigma.hard.true   0.128 0.028
# sigma.easy.true   0.049 0.012
# sigma.hard.false  0.234 0.082
# sigma.easy.false  0.245 0.084
# tau.hard.true     0.312 0.051
# tau.easy.true     0.257 0.022
# tau.hard.false    1.117 0.553
# tau.easy.false    2.699 1.108
# muS               0.217 0.018
# sigmaS            0.029 0.020
# tauS              0.026 0.017
# tf               -2.235 0.385
# gf               -3.583 0.563

# Suppose we wished to make the priors vaguer by doubling their standard deviation.
# This can be done by altering the ppar output then using it to change the post.prior object
ppars$sd <- ppars$sd*2
post.prior1 <- assign.ppars(post.prior1,ppars)
plot.dmc(hsamples[[1]],p.prior=post.prior1,layout=c(3,6))

# Suppose we wanted to use a neater set of values
ppars$mean <- c(.55,.45,.75,1,.1,.05,.2,.2,.3,.25,1,2,.2,.03,.03,-2,-3)
ppars$sd   <- c(.1,.04,.3,.3,.05,.05,.2,.2,.1,.1,1,2,.04,.04,.04,1,1)
post.prior1 <- assign.ppars(post.prior1,ppars)

# A plot shows these look fairly similar
plot.dmc(hsamples[[1]],p.prior=post.prior1,layout=c(3,6))

# We can do the same thing at the population level. Here we use the scale.sd argument
# to make.tnorm.prior to double the standard deviation of the priors in one step.
# This gives the new hyper priors for the population means (locations):
post.prior.location <- make.tnorm.prior(pp.prior[[1]],hsamples,scale.sd=2)
# This gives the new hyper priors for the population standard deviations (scales);
# Here you have to specify hpar=2
post.prior.scale <- make.tnorm.prior(pp.prior[[2]],hsamples,scale.sd=2,hpar=2)

# The two new priors are then combined to make a hyper prior
post.pp.prior <- list(post.prior.location,post.prior.scale)

# The makes Figure 15a and 15b  
plot.dmc(hsamples,p.prior=post.pp.prior,layout=c(3,6),location=TRUE)
plot.dmc(hsamples,p.prior=post.pp.prior,layout=c(3,6),scale=TRUE)

# We might also want to obtain priors at the participant level.
# To do so, wefirst collapse all participant-level posterior samples
theta <- collapse.subjects(hsamples)$theta

# Fitting and plotting can take a little while as there are a large number of 
# samples (~700K)
post.prior.individual <- make.tnorm.prior(p.prior,theta,scale=2)

# To obtain a plot, we also need to let plot.dmc know to collapse, producing
# Figure 15c
plot.dmc(hsamples,layout=c(3,6),p.prior=post.prior.individual,collapse.subjects=TRUE)

# Save the values to make Table 1
ppars.location <- get.ppars(make.tnorm.prior(pp.prior[[1]],hsamples))
ppars.scale <- get.ppars(make.tnorm.prior(pp.prior[[2]],hsamples,hpar=2))
ppars.individual <- get.ppars(make.tnorm.prior(p.prior,theta))
write.csv(round(ppars.location,3),"tmp.csv")
write.csv(round(ppars.scale,3),"tmp1.csv")
write.csv(round(ppars.individual,3),"tmp2.csv")

# save_data(datgf,samples,hsamples,pp,pp.sim,file="DMCpaper2.RData")

