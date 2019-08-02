##################  DMC Lesson 4: Multiple Subjects


### Lesson 4.1 Simulating multiple subjects

rm(list=ls())

# Current working directory must be set to the top-level folder  
source ("dmc/dmc.R")
load_model ("LBA","lba_B.R")

# Same as the single-subject LBA example in the early single-subject lessons, 
# but now with multiple subjects, where the subjects factor is a fixed effect 
# (i.e., each subject is dealt with separately).

# Model setup and parameter vector are as before.
model <- model.dmc(p.map=list(A="1",B="R",t0="1",mean_v="M",sd_v="M",st0="1"),
                   match.map=list(M=list(s1=1,s2=2)),
                   constants=c(sd_v.false=1,st0=0),
                   factors=list(S=c("s1","s2")),
                   responses=c("r1","r2"),
                   type="norm")

p.vector  <- c(A=3,B.r1=4,B.r2=5,
               mean_v.true=2,mean_v.false=-1, sd_v.true=0.5, t0=.2)

# First lets simulate ns=2 subjects data, with 1 data point per design cell
h.simulate.dmc(model,ps=p.vector,ns=2,n=1)
#   s  S  R       RT
# 1 1 s1 r1 2.560937
# 2 1 s2 r2 3.451212
# 3 2 s1 r1 3.175386
# 4 2 s2 r2 3.197400

# The "s" column is a factor ("s" is reserved for subject IDs and cannot be 
# used as a factor name)

# In general, n can be a single number for a balanced design or matrix for an 
# unbalanced design, where rows are subjects and columns are design cells. If 
# the matrix has one row then all subjects have the same n in each cell. If it 
# has one column then all cells have the same n, otherwise each entry specifies 
# the n for a particular design subject x design cell combination. Following
# are some examples.

# To get more data for one subject than another use an ns row, one column, 
# matrix for n, here n=1 per cell for s=1 and n=2 per cell for s=2
h.simulate.dmc(model,ps=p.vector,ns=2,n=matrix(c(1:2),ncol=1))

# As in simulate.dmc there can be a different n per cell. Here with equal
# n for each participant this is achieved with a 1 row, n cell column matrix.
# For example this design has 2 cells, n=1 for s1 and n=2 for s2 is given by:
h.simulate.dmc(model,ps=p.vector,ns=2,n=matrix(c(1:2),nrow=1))

# Uneven n per cells and subjects requires a matrix, e.g., 
# Specify ps as a ns row matrix to get different parameters for each subject
# s=1: s1=1, s2=2, s=2: s1=3, s2=4 
h.simulate.dmc(model,ps=p.vector,ns=2,n=matrix(c(1:4),nrow=2,byrow=TRUE))

# Can also have different parameters for different participants by making ps 
# a matrix with ns rows. 
data <- h.simulate.dmc(model,ps=rbind(p.vector,p.vector*.9),ns=2,n=1)

# Note that the parameters for each subject are saved as an attribute matrix
# with one row per subject (whether ps is a vector or matrix)
attr(data,"parameters")
#     A B.r1 B.r2 mean_v.true mean_v.false sd_v.true   t0
# 1 3.0  4.0  5.0         2.0         -1.0      0.50 0.20
# 2 2.7  3.6  4.5         1.8         -0.9      0.45 0.18


# data.model.dmc recognizes the presence of the "s" column and outputs a list
# where each element is a single subject. 
data.model <- data.model.dmc(data,model)
data.model
# $`1`
#    S  R       RT
# 1 s1 r1 2.553215
# 2 s2 r2 3.377843
# 
# $`2`
#    S  R       RT
# 3 s1 r1 3.398975
# 4 s2 r2 3.518320

# This allows discrete processing of subjects, for example: 
for (s in names(data.model))
  print(likelihood.dmc(attr(data.model[[s]],"parameters"),data.model[[s]]))

# SUBJECTS AS RANDOM EFFECTS
#
# Each subject can be thought of as a random effect. That is, rather than 
# thinking of each subject as being completely unrelated to all other subjects
# we treat them as coming from the same population. In practice this means that
# each subject's parameters are from a population distribution.

# This can be done with the "p.prior" (parameter prior) argument to 
# h.simulate.dmc, with exactly the same format as the prior objects introduced
# in lesson 2 (i.e., lists constructed by the prior.p.dmc function that defines 
# distributions from which parameters can be sampled. We will use the default 
# distribution (normal), requiring a mean (p1) and standard deviation (p2) 
# parameter. Here we use p.vector to give the means and the same sd=0.025 for 
# all parameters (a vector for p2 allows different sds for each parameter).
p.prior <- prior.p.dmc(p1=p.vector,p2=.025)

# You can plot the population distributions, for example the t0 distribution is:  
par(mfrow=c(2,4))
plot.prior("A",p.prior)
plot.prior("B.r1",p.prior)
plot.prior("B.r2",p.prior)
plot.prior("mean_v.true",p.prior)
plot.prior("mean_v.false",p.prior)
plot.prior("sd_v.true",p.prior)
plot.prior("t0",p.prior)

# Save off data from two subjects 
data <- h.simulate.dmc(model,p.prior=p.prior,ns=2,n=1)
# The sampled parameters are
attr(data,"parameters")


# Now lets try this with a larger set of subjects
data <- h.simulate.dmc(model,p.prior=p.prior,ns=100,n=100)

# The parameter distributions look like the priors
par(mfrow=c(2,4))
for (i in dimnames(attr(data,"parameters"))[[2]])
  plot(density(attr(data,"parameters")[,i]),xlab=i,main="")

# As do their means
round(apply(attr(data,"parameters"),2,mean),2)

# and standard deviations
round(apply(attr(data,"parameters"),2,sd),3)

###################### GETTING REAL DATA INTO DMC ##############################
#
# Setup your data frame like the output of h.simulate.dmc.
# 1) Make sure you have an "s" column identifying the data for each subject that
#    is a factor
# 2) Follow that by columns specifying the factors in your design. Again each 
#    must be a factor
# 3) The second last column give the response (again a factor)
# 4) The last columns is called RT, make sure the unit is seconds.
#
# Column names and factor levels for (2) and (3) depend on your model 
# specification so make sure they are compatible. data.model.dmc attempts to 
# check this.
#
# Note that model.dmc is very picky about overlap in factor levels, particularly
# where one level is a substring of another. This means you should avoid short
# and cryptic level names (this also helps with reading output!).
#
# Avoid superfluous columns, they are usually ignored but can sometimes cause
# problems.
#
# See advanced lesson dmc_5_1_Factorial.R for setting up models for the complex
# designs that often occur in real experiments.
