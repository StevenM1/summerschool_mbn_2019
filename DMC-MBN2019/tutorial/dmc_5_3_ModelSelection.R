##################  DMC Lesson 5: Advanced Multiple Subjects

#### Lesson 5_3: Hierarchical Model selection with WAIC, DIC and BPIC, LNR model


# Here we compare two LNR models for a design with a stimulus manipulation (S)
# and a second manipulation (F). In model 0 neither factor has an effect,
# in model 1 F has an effect on meanlog. Unlike lesson 3.5, we only use the 
# WAIC criteria to compare two random effects (hierarchical) LNR models.

#  N.B: The two models are estimated as random effects models. However, we 
#       conduct model selection using a fixed effects version of WAIC. WAIC 
#       measures how well the model predicts new data. Fixed effects WAIC
#       is how well the model predicts new data from the SAME participants.

# All of the fitting was done beforehand (see commented out code at the end of 
# this file). Modify that code if you want to redo it yourself.

rm(list=ls())
# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("LNR","lnrPP.R")

load_data ("dmc_5_3.RData")


# #### Plot Data ####

# posterior predictives for each subject
pp00 <- h.post.predict.dmc(hsamples0.0)
pp11 <- h.post.predict.dmc(hsamples1.1)
pp01 <- h.post.predict.dmc(hsamples0.1)
pp10 <- h.post.predict.dmc(hsamples1.0)

# Goodness of fit can be assessed by plotting the cdf for the empirical data
# and posterior predictive averaged over subjects. 

plot.pp.dmc(pp00) # Model 0 Predictives and Data 0
plot.pp.dmc(pp11) # Model 1 Predictives and Data 1
plot.pp.dmc(pp01) # Model 0 Predictives and Data 1
plot.pp.dmc(pp10) # Model 1 Predictives and Data 0
# And for subjects e.g., 
s=1
plot.pp.dmc(pp00[[s]],"cdf") # Model 0 Predictives and Data 0
plot.pp.dmc(pp11[[s]],"cdf") # Model 1 Predictives and Data 1
plot.pp.dmc(pp01[[s]],"cdf") # Model 0 Predictives and Data 1
plot.pp.dmc(pp10[[s]],"cdf") # Model 1 Predictives and Data 0

#### Model Selection M0D0 ####

# To calculate WAIC with the first method (see Lesson 4_4) we need to create 
# iterations by trials pointwise log-likelihood matrix (i.e., points
# are trials). Recall these can be big objects, here thinning is set high
# for speed, reducing the size of the objects to 37MB (without thinning they 
# are 1.83 GB!). The effect on estimates is only small.
group_trial_0.0 <- group_trial_log_likes(hsamples0.0,thin=50)

# Note that you can run e.g.,
#  group_trial_log_likes(hsamples0.0,thin=50,get.size=TRUE)
# to get a quick estimate of the size of the resulting object.

# This can be a bit slow. To speed up you can use multicore. In this case it 
# doesn’t save much time, but in slower cases it may help, although be careful 
# about the amount of RAM that is used. The first subject is done with a single
# core and the size of the output checked, then the multicore run is launched
# for the remaining subjects.
group_trial_0.0 <- group_trial_log_likes(hsamples0.0,thin=50,cores=4)
group_trial_1.1 <- group_trial_log_likes(hsamples1.1,thin=50,cores=4)
group_trial_0.1 <- group_trial_log_likes(hsamples0.1,thin=50,cores=4)
group_trial_1.0 <- group_trial_log_likes(hsamples1.0,thin=50,cores=4)

# For the second method we create an iterations by subjects pointwise 
# log-likelihoods matrix (i.e., points are subjects). Remember to only use
# this method with a large subject sample size
group_subject_0.0 <- group_subject_log_likes(hsamples0.0)
group_subject_0.1 <- group_subject_log_likes(hsamples0.1)
group_subject_1.0 <- group_subject_log_likes(hsamples1.0)
group_subject_1.1 <- group_subject_log_likes(hsamples1.1)


# The usual waic function can then be used to calculate the group waic
M0D0_G <- waic.dmc(group_trial_0.0,save=T)
# Again the usual waic function can then be used to calculate the group waic
M0D0_S <- waic.dmc(group_subject_0.0,save=T)

# In this case with a fair number of subjects the two methods are
# in reasonable agreement with respect to WAIC but not in terms of number of
# parameters. The trial version is close to the nominal value (40 x 5 =200).
cbind(M0D0_G,M0D0_S)
#              M0D0_G         M0D0_S     
# elpd_waic    -81854.62      -81822.52  
# p_waic       190.4481       98.4521    
# waic         163709.2       163645     
# se_elpd_waic 351.0542       3307.189   
# se_p_waic    3.19809        1.253472   
# se_waic      702.1084       6614.379   
# ...

# Similar results are obtained for the other models

M1D1_G <- waic.dmc(group_trial_1.1,save=T)
M1D1_S <- waic.dmc(group_subject_1.1,save=T)
cbind(M1D1_G,M1D1_S)
#              M1D1_G         M1D1_S     
# elpd_waic    -66696.78      -66655.92  
# p_waic       264.5097       138.3961   
# waic         133393.6       133311.8   
# se_elpd_waic 365.6809       2277.159   
# se_p_waic    3.135385       1.343365   
# se_waic      731.3618       4554.319   
# ...

M0D1_G <- waic.dmc(group_trial_0.1,save=T)
M0D1_S <- waic.dmc(group_subject_0.1,save=T)
cbind(M0D1_G,M0D1_S)
#              M0D1_G         M0D1_S     
# elpd_waic    -71412.12      -71383.31  
# p_waic       185.8336       98.46823   
# waic         142824.2       142766.6   
# se_elpd_waic 365.5664       2130.272   
# se_p_waic    2.772893       1.175869   
# se_waic      731.1327       4260.544   
# ...

M1D0_G <- waic.dmc(group_trial_1.0,save=T)
M1D0_S <- waic.dmc(group_subject_1.0,save=T)
cbind(M1D0_G,M1D0_S)
#              M1D0_G         M1D0_S     
# elpd_waic    -81892.09      -81856.98  
# p_waic       255.2377       137.1047   
# waic         163784.2       163714     
# se_elpd_waic 351.0943       3305.995   
# se_p_waic    3.23938        1.863431   
# se_waic      702.1885       6611.99    
# ...

# We have calculated a pointwise matrix at the trial and subject level for 
# each of the four models. So now we can use the "loocompare.dmc" function to
# compare the outputs of waic.dmc for two models.

# First lets use the trial level WAIC calculations. We will compare the Model 0 
# and Model 1, which have both been fit to Data 0. Recall that the calculation 
# is (second - first), so larger values 
# mean that the first model is more adequate.

# Here we have a strong indication that the true model (Model 0), is more 
# adequate than the false model (Model 1).
loocompare.dmc(M0D0_G,M1D0_G,digits=3)
#         waic_diff        se 
#              74.9      16.5 

# Now lets compare the Model 1 and Model 0, which have both been fit to
# Data 1. As expected we have very strong evidence that the true model (Model 1) 
# is more adequate than the false model (Model 0)
loocompare.dmc(M1D1_G,M0D1_G,digits=3)
#         waic_diff        se 
#              9431       179

# Now lets look at the WAIC results for the subject level. We have a 40 subjects
# and should see the stability of the estimates decrease (higher standard error) 

# Again the true model (Model 0) is more adequate than the false model (Model 1).
# However, the standard error estimates suggest we have not lost much stability 
# by using the subject level WAIC
loocompare.dmc(M0D0_S,M1D0_S,digits=3)
#         waic_diff        se 
#              68.9      12.1

# Model weights confirm a strong preference for the simpler model. 
loocompare.dmc(list(m0d0=M0D0_S,m1d0=M1D0_S),digits=3)
#        m0d0     m1d0
# IC-min    0 6.89e+01
# w         1 1.08e-15

# For the second comparison the true model (Model 1) is more adequate than the 
# false model (Model 0). The result is like the trial level comparison, but
# we can see that the standard error of the estimate has increased by a factor 
# of about 6. This demonstrates the stability we lose when using subjects in 
# the pointwise calculation.
loocompare.dmc(M1D1_S,M0D1_S,digits=3)
#         waic_diff        se 
#              9455      1095 


# save_data(pp00,hsamples0.0,pp11,hsamples1.1,pp01,hsamples0.1,pp10,hsamples1.0,
#      file="dmc_5_3.RData")

#### DIC model selection.

# The following function applies IC.dmc to each subject and prints the sum,
# and invisibly returns the results for each subject. This analysis is fast
# unless fast=FALSE, in which case the posterior likelihoods are re-calculated.
# It returns the minimum summed deviance (a measure of best fit) and one of 
# two ICs

# By default, BPIC is calculated:

h.IC.dmc(hsamples0.0)
# Summed Minimum Deviance and BPIC
# [1] 163320.0 163895.5
h.IC.dmc(hsamples1.0)
# Summed Minimum Deviance and BPIC
# [1] 163262.6 164057.7

# The true model is chosen by a wide margin 164057.7-163895.5 = 162.2. This is
# due to the complexity penalty as in terms of pure fit the more complex model
# is better: 163320-163262.6 = 57.4

h.IC.dmc(hsamples1.1)
# Summed Minimum Deviance and BPIC
# [1] 132853.3 133656.7
h.IC.dmc(hsamples0.1)
# Summed Minimum Deviance and BPIC
# [1] 142441.3 143018.9

# The true model is chosen by an even wider margin in this easy case
# 143018.9-133656.7 = 9362.2, due to much better fit 142441.3-132853.3 = 9588,
# which the complexity penalty does not approach.

# Now DIC

h.IC.dmc(hsamples0.0,DIC=TRUE)
# Summed Minimum Deviance and DIC
# [1] 163320.0 163702.1
h.IC.dmc(hsamples1.0,DIC=TRUE)
# Summed Minimum Deviance and DIC
# [1] 163262.6 163791.0

# As expeced the margin is narrower given the weaker simplicity penalty
# 163791-163702.1 = 88.9

h.IC.dmc(hsamples1.1,DIC=TRUE)
# Summed Minimum Deviance and DIC
# [1] 132853.3 133388.4
h.IC.dmc(hsamples0.1,DIC=TRUE)
# Summed Minimum Deviance and DIC
# [1] 142441.3 142824.7
# and is wider when the true model is more complex
# 142824.7-133388.4 = 9436.3


# #### DATA AND FIT GENERATION
# 
# n.trials.per.cell = 5e2
# 
# #### Model 0 ####
# model0 <- model.dmc(p.map=list(meanlog= "M",sdlog="M",t0="1",st0="1"),
#                     match.map=list(M=list(s1="r1",s2="r2")),
#                     factors=list(S=c("s1","s2"),F=c("f1","f2")),
#                     constants=c(st0=0),
#                     responses=c("r1","r2"),
#                     type="lnr")
# 
# # Population distribution (used to simulate data not as a prior)
# p.mu0  <- c(meanlog.true=-.5,meanlog.false=.5,# natural scale
#             sdlog.true=log(1),sdlog.false=log(1),t0=log(.2)) # log scale
# 
# # Fairly tight distributions
# p.sigma0 <- c(.2,.2,.2,.2,.1); names(p.sigma0) <- names(p.mu0)
# p.prior0 <- prior.p.dmc(p1=p.mu0,p2=p.sigma0,
#                         untrans=c(sdlog.true="exp",sdlog.false="exp",t0="exp"))
# 
# # plot population distributions
# par(mfcol=c(2,3)); for (i in names(p.prior0)) plot.prior(i,p.prior0)
# 
# # Data
# data0 <- h.simulate.dmc(model0,p.prior=p.prior0,n=n.trials.per.cell,ns=40)
# 
# # Location prior
# p.mu.mu0 <- p.mu0
# p.mu.sigma0 <- p.sigma0*10
# mu.prior0 <- prior.p.dmc(p1=p.mu.mu0,p2=p.mu.sigma0,
#                          untrans=c(sdlog.true="exp",sdlog.false="exp",t0="exp"))
# # Scale prior
# p.sigma.shape0 <- rep(1,length(p.sigma0)); names(p.sigma.shape0) <- names(p.sigma0)
# p.sigma.scale0 <- rep(1,length(p.sigma0))
# sigma.prior0 <- prior.p.dmc(p1=p.sigma.shape0,p2=p.sigma.scale0,
#                             dists=rep("gamma",length(p.sigma0)))
# 
# # Make a hyper-prior list
# pp.prior0=list(mu.prior0,sigma.prior0)
# 
# #### Model 1 ####
# model1 <- model.dmc(p.map=list(meanlog= c("F","M"),sdlog="M",t0="1",st0="1"),
#                     match.map=list(M=list(s1="r1",s2="r2")),
#                     factors=list(S=c("s1","s2"),F=c("f1","f2")),
#                     constants=c(st0=0),
#                     responses=c("r1","r2"),
#                     type="lnr")
# 
# # Population distribution (used to simulate data not as a prior)
# p.mu1  <- c(meanlog.f1.true=-.5,meanlog.f2.true=-1,# natural scale
#             meanlog.f1.false=.5,meanlog.f2.false=0,# natural scale
#             sdlog.true=log(1),sdlog.false=log(1),t0=log(.2)) # log scale
# 
# # Fairly tight distributions
# p.sigma1 <- c(.2,.2,.2,.2,.2,.2,.1); names(p.sigma1) <- names(p.mu1)
# p.prior1 <- prior.p.dmc(p1=p.mu1,p2=p.sigma1,
#                         untrans=c(sdlog.true="exp",sdlog.false="exp",t0="exp"))
# 
# # plot population distributions
# par(mfcol=c(3,3)); for (i in names(p.prior1)) plot.prior(i,p.prior1)
# 
# # Data
# data1 <- h.simulate.dmc(model1,p.prior=p.prior1,n=n.trials.per.cell,ns=40)
# 
# # Location prior
# p.mu.mu1 <- p.mu1
# p.mu.sigma1 <- p.sigma1*10
# mu.prior1 <- prior.p.dmc(p1=p.mu.mu1,p2=p.mu.sigma1,
#                          untrans=c(sdlog.true="exp",sdlog.false="exp",t0="exp"))
# # Scale prior
# p.sigma.shape1 <- rep(1,length(p.sigma1)); names(p.sigma.shape1) <- names(p.sigma1)
# p.sigma.scale1 <- rep(1,length(p.sigma1))
# sigma.prior1 <- prior.p.dmc(p1=p.sigma.shape1,p2=p.sigma.scale1,
#                             dists=rep("gamma",length(p.sigma1)))
# 
# # Make a hyper-prior list
# pp.prior1=list(mu.prior1,sigma.prior1)
# 
# #### TRUE MODEL FITS 
# 
# #### FITTING MODEL0 TO MODEL0 DATA ####
# 
# data00 <- data.model.dmc(data0,model0)
# hsamples0.0 <- h.samples.dmc(nmc=100,p.prior0,data00,thin=10,pp.prior=pp.prior0)
# hsamples0.0 <- h.run.dmc(hsamples0.0,cores=4,report=10,blocks=NA,
#                          p.migrate=0.05,h.p.migrate=0.05)
# hsamples0.0 <- h.run.dmc(h.samples.dmc(nmc=200, samples=hsamples0.0,
#                                        thin=10), cores=4,report=10)
# # Convergence is rapid
# plot.dmc(hsamples0.0,hyper=TRUE,layout=c(2,4),smooth=FALSE)
# # R hat 
# gelman.diag.dmc(hsamples0.0,hyper=TRUE)
# # Things also look good at the individual participant level
# gelman.diag.dmc(hsamples0.0)
# # Generate summary of hyper parameters
# summary.dmc(hsamples0.0,hyper=TRUE)
# # Reasonably large effective sample size
# effectiveSize.dmc(hsamples0.0,hyper=TRUE)
# 
# #### FITTING MODEL1 TO MODEL1 DATA ####
# 
# data11 <- data.model.dmc(data1,model1)
# hsamples1.1 <- h.samples.dmc(nmc=100,p.prior1,data11,thin=10,pp.prior=pp.prior1)
# hsamples1.1 <- h.run.dmc(hsamples1.1,cores=4,report=10,blocks=NA,
#                          p.migrate=0.05,h.p.migrate=0.05)
# hsamples1.1 <- h.run.dmc(h.samples.dmc(nmc=200, samples=hsamples1.1,
#                                        thin=10), cores=4,report=10)
# # Convergence is rapid
# plot.dmc(hsamples1.1,hyper=TRUE,layout=c(2,4),smooth=FALSE)
# # R hat 
# gelman.diag.dmc(hsamples1.1,hyper=TRUE)
# # Things also look good at the individual participant level
# gelman.diag.dmc(hsamples1.1)
# # Generate summary of hyper parameters
# summary.dmc(hsamples1.1,hyper=TRUE)
# # Reasonaby large effective sample size
# effectiveSize.dmc(hsamples1.1,hyper=TRUE)
# 
# #### FITTING MODEL0 TO MODEL1 DATA ####
# 
# data01 <- data.model.dmc(data1,model0)
# hsamples0.1 <- h.samples.dmc(nmc=100,p.prior0,data01,thin=10,pp.prior=pp.prior0)
# hsamples0.1 <- h.run.dmc(hsamples0.1,cores=4,report=10,blocks=NA,
#                          p.migrate=0.05,h.p.migrate=0.05)
# hsamples0.1 <- h.run.dmc(h.samples.dmc(nmc=200, samples=hsamples0.1,
#                                        thin=10), cores=4,report=10)
# # Convergence is rapid
# plot.dmc(hsamples0.1,hyper=TRUE,layout=c(2,4),smooth=FALSE)
# # R hat 
# gelman.diag.dmc(hsamples0.1,hyper=TRUE)
# # Things also look good at the individual participant level
# gelman.diag.dmc(hsamples0.1)
# # Generate summary of hyper parameters
# summary.dmc(hsamples0.1,hyper=TRUE)
# # Reasonaby large effective sample size
# effectiveSize.dmc(hsamples0.1,hyper=TRUE)
# 
# #### FITTING MODEL1 TO MODEL0 DATA ####
# 
# data10 <- data.model.dmc(data0,model1)
# hsamples1.0 <- h.samples.dmc(nmc=100,p.prior1,data10,thin=10,pp.prior=pp.prior1)
# hsamples1.0 <- h.run.dmc(hsamples1.0,cores=4,report=10,blocks=NA,
#                          p.migrate=0.05,h.p.migrate=0.05)
# hsamples1.0 <- h.run.dmc(h.samples.dmc(nmc=200, samples=hsamples1.0,
#                                        thin=10), cores=4,report=10)
# # Convergence is rapid
# plot.dmc(hsamples1.0,hyper=TRUE,layout=c(2,4),smooth=FALSE)
# # R hat 
# gelman.diag.dmc(hsamples1.0,hyper=TRUE)
# # Things also look good at the individual participant level
# gelman.diag.dmc(hsamples1.0)
# # Generate summary of hyper parameters
# summary.dmc(hsamples1.0,hyper=TRUE)
# # Reasonably large effective sample size
# effectiveSize.dmc(hsamples1.0,hyper=TRUE)


