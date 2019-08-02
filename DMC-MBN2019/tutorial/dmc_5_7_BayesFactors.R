##################  DMC Lesson 5: Advanced Model Selection

### Lesson 5.7:  Model Selection Using Bayes Factors 

# This lesson is based on the article Gronau, Q. F., Heathcote, A., & Matzke, D.
# (2018). Computing Bayes factors for evidence-accumulation models using
# Warp-III bridge sampling. Manuscript submitted for publication. Available from
# https://psyarxiv.com/9g4et

rm(list=ls()) 
# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")

# load_data ("dmc_5_7.RData")

# This lesson gives an example comparing nested single-participant LBA models 
# and an example comparing non-nested hierarchical LBA models. Non-nested 
# single-participant and nested hierarchical comparisons are also possible. 
# Furthermore, DMC allows users to compute Bayes factors for any implemented 
# model, not just the LBA.


###################  SINGLE-PARTICIPANT EXAMPLE ################################

# We will first use the single-participant example from Gronau et al. (2018).
# In this example, we compute a Bayes factor to compare a simple LBA model
# (referred to as "full" model) to a restricted version of the model that fixes
# the parameter mean_v.true to the constant 3.55.
# The Bayes factor is given by the ratio of the models' marginal likelihoods.
# In DMC, we can compute the (log) marginal likelihood of a model using Warp-III
# bridge sampling. Once the marginal likelihood has been obtained for both
# models, we can compute the Bayes factor of interest.
load_model ("LBA","lba_B.R")

# First, let us specify the full model
LBAfull <- model.dmc(
  p.map = list(A = "1", B = "1", mean_v = "M", sd_v = "1",
               t0 = "1", st0 = "1"),
  constants = c(st0 = 0, sd_v = 1),
  match.map = list(M = list(s1 = 1, s2 = 2)),
  factors = list(S = c("s1", "s2")),
  responses = c("r1", "r2"),
  type = "norm")

# Select parameter values for data generation
pLBAfull  <- c(A = .5, B = 1, mean_v.true = 4, mean_v.false = 3, t0 = .2)

# Generate a data set with 250 trials per stimulus (i.e., a total of 500 trials)
# from the full model.
set.seed(1) # So results reproducible
data <- simulate.dmc(pLBAfull, LBAfull, n = 250)

# Bind the data and model together ready for sampling
dmLBAfull <- data.model.dmc(data, LBAfull)

# Specify priors
p.priorLBAfull <- prior.p.dmc(
  p1 = c(A = 1, B = 1, mean_v.true = 2, mean_v.false = 1, t0 = .3),                       
  p2 = c(1, 1, 3, 3, .25), lower = c(0, 0, NA, NA, .1), upper = rep(NA, 5)
)

# Create object to fill with posterior samples for the full model
sLBAfull <- samples.dmc(nmc = 300, p.priorLBAfull, dmLBAfull, thin = 10)

# Obtain posterior samples for the full model
sLBAfull <- RUN.dmc(sLBAfull, cores = 4)

# Check convergence
plot.dmc(sLBAfull)
gelman.diag.dmc(sLBAfull)
# Potential scale reduction factors:
# 
#              Point est. Upper C.I.
# A                 1.003       1.01
# B                 0.999       1.00
# mean_v.true       1.003       1.01
# mean_v.false      1.003       1.01
# t0                0.999       1.00
# 
# Multivariate psrf
# 
# 1.01

# The chains look good and the gelman.diag values indicate convergence

# We can now compute the (log) marginal likelihood for the full model using
# Warp-III bridge sampling. Note that iterations are reported on completion.
bridgeLBAfull <- bridge.sampler.dmc(samples = sLBAfull, cores = 4)

# bridge.sampler.dmc has the following arguments:
#
# samples: An object with posterior samples for a single participant.
# cores: The number of cores to use for the computation.
# repetitions: Allows the user to repeat the (log) marginal likelihood
#              estimation several times, each time based on new proposal
#              samples. 
#              N.B.: Each repetition will use the same posterior samples and,
#              ideally, we recommend to instead assess the estimation
#              uncertainty by repeatedly drawing posterior samples and then
#              computing the marginal likelihood for each set of posterior
#              samples separately.
# use_neff: Logical that specifies whether or not to use the effective sample
#           size in the optimal bridge function. Default TRUE
# silent: Logical that specifies whether or not to print the number of
#         iterations to the console. Default FALSE
# maxiter: Maximum number of iterations for the iterative scheme. Default 500
# r0: Specifies the starting value for the iterative scheme. N.B.: Due to the
#     more numerically stable implementation, this does not correspond exactly 
#     to the initial log marginal likelihood value. In general, we recommend to
#     not change this value.
# tol1: Specifies the criterion when to stop the iterative scheme. We recommend
#       to not change this value.
# tol2: Specifies the criterion when to stop the iterative scheme in case it did
#       not converge within maxiter (i.e., when the geometric mean trick is
#       used). We recommend to not change this value.
#
# See the following for a tutorial on bridge sampling.
# Gronau, Q. F., Sarafoglou, A., Matzke, D., Ly, A., Boehm, U., Marsman, M., 
# et al. (2017). A tutorial on bridge sampling, Journal of Mathematical 
# Psychology, 81, 80–97. http://doi.org/10.1016/j.jmp.2017.09.005

# Next, we specify the restricted model that fixes mean_v.true to 3.55
LBAres <- model.dmc(
  p.map = list(A = "1", B = "1", mean_v = "M", sd_v = "1",
               t0 = "1", st0 = "1"),
  constants = c(mean_v.true = 3.55, st0 = 0, sd_v = 1),
  match.map = list(M = list(s1 = 1, s2 = 2)),
  factors = list(S = c("s1", "s2")),
  responses = c("r1", "r2"),
  type = "norm")

# Bind the data and model together ready for sampling
dmLBAres <- data.model.dmc(data, LBAres)

# Specify priors for the restricted model
p.priorLBAres <- prior.p.dmc(
  p1 = c(A = 1, B = 1, mean_v.false = 1, t0 = .3),                           
  p2 = c(1, 1, 3, .25), lower = c(0, 0, NA, .1), upper = rep(NA, 4)
)

# Create object to fill with posterior samples for the restricted model
sLBAres <- samples.dmc(nmc = 300, p.priorLBAres, dmLBAres, thin = 10)

# Obtain posterior samples for the restricted model
sLBAres <- RUN.dmc(sLBAres, cores = 4)

# Check convergence
plot.dmc(sLBAres)
gelman.diag.dmc(sLBAres)
# Potential scale reduction factors:
# 
#              Point est. Upper C.I.
# A                     1          1
# B                     1          1
# mean_v.false          1          1
# t0                    1          1
# 
# Multivariate psrf
# 
# 1.01

# Again, the chains look good and the gelman.diag values indicate convergence

# We can now compute the (log) marginal likelihood for the restricted model
# using Warp-III bridge sampling
bridgeLBAres <- bridge.sampler.dmc(samples = sLBAres, cores = 4)


# Finally, we can compute the Bayes factor in favor of the full model
bf_full_res <- bf.dmc(bridgeLBAfull, bridgeLBAres)
print(bf_full_res)
# [1] 0.7608655

# The Bayes factor indicates that the data are ~ 0.8 times more likely under the
# full model than under the restricted model. This Bayes factor value is very
# close to 1, indicating that the data are almost as likely under the full model
# as under the restricted model (i.e., the evidence is ambiguous).

# In general, it may be easier to interpret Bayes factor values larger than 1
# which can be achieved by re-expressing the Bayes factor in favor of the
# restricted model. Note that the Bayes factor in favor of the full model and
# the Bayes factor in favor of the restricted model convey the same information
# (bf_full_res = 1 / bf_res_full); re-expressing the Bayes factor so that it is
# larger than 1 is simply a matter of convenience.
bf_res_full <- bf.dmc(bridgeLBAres, bridgeLBAfull)
print(bf_res_full)
 bridgeLBAfull$logml# [1] 1.314293

# Hence, expressed differently, the Bayes factor indicates that the data are
# ~ 1.2 times more likely under the restricted model than under the full model


###################  HIERARCHICAL EXAMPLE ######################################

# We will use one of the hierarchical examples from Gronau et al. (2018).
# Specifically, we will use data simulated from the LBA for a design with four
# cells, two conditions that differed in the threshold B parameter crossed with
# two stimuli, and two possible responses. This data set features
# 20 participants, each providing 200 trials per cell of the design. We will
# compare the data-generating LBA model that allows only threshold B to be
# different across conditions (i.e., B-model) to a version of the LBA that
# allows only mean drift rate mean_v.true to be different across conditions
# (i.e., V-model).

# The Bayes factor is given by the ratio of the models' marginal likelihoods.
# In DMC, we can compute the (log) marginal likelihood of a model using Warp-III
# bridge sampling. Once the marginal likelihood has been obtained for both
# models, we can compute the Bayes factor of interest.

load_model ("LBA","lba_B2v.R")

# We assume that the user has already obtained converged posterior samples for
# both models. Note that, in general, for the computation of Bayes factors using
# Warp-III bridge sampling, it is advisable to obtain more posterior samples
# than for parameter estimation. As a rough guideline, for the hierarchical
# case, a good starting point might be to aim for at least 20,000-30,000 samples
# (collapsed across all chains), with not too high autocorrelation (e.g., an
# effective size of 15,000 or greater).

# We can now compute the (log) marginal likelihood for the data-generating
# B-model using Warp-III bridge sampling. The (log) marginal likelihood for 
# the misspecified V-model can be computed in analogous fashion

# Load the posterior samples for the two models and perform bridge sampling
# WARNING: These objects are LARGE, ~ 142MB, and the computations will take a 
# few hours each! To continue the lesson quicker, users can skip the next four 
# lines and simply use the bridgeB and bridgeV results stored in dmc_5_7.RData. 
load_data("BBl_1.RData")
load_data("BVl_1.RData")
bridgeB <- h.bridge.sampler.dmc(samples = hVVl_1, cores = 4)
bridgeV <- h.bridge.sampler.dmc(samples = hBVl_1, cores = 4)

# Now we can compute the Bayes factor in favor of the data-generating B-model
bf_B_V <- bf.dmc(bridgeB, bridgeV)
print(bf_B_V)
# [1] Inf

# The Bayes factor indicates overwhelming evidence in favor of the
# data-generating B-model over the misspecified V-model

# In case the Bayes factor is very large or even (numerically) Inf, it may be
# useful to consider the log Bayes factor that can also be obtained using bf.dmc
logbf_B_V <- bf.dmc(bridgeB, bridgeV, log = TRUE)
print(logbf_B_V)
# [1] 903.9653

# save_data(sLBAfull, bridgeLBAfull,sLBAres,bridgeLBAres,
#           bridgeB,bridgeV,file="dmc_5_7.RData")

