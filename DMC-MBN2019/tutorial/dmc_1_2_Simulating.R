##################  DMC Lesson 1 Models and Data
#
# Note: Before running these scripts, the current working directory must be set 
# to the top-level folder containing the dmc and tutorial subfolders 
#
# !!! RSTUDIO USERS: Plots can cause the error message: "Error in plot.new() : 
# figure margins too large". This is not a DMC bug, it just means your screen
# size is limited. You can fix this by maximising the plot pane before 
# plotting, or reducing the number of panels and plotting multiple windows 
# using par(mfrow=c(x,y)) with smaller x and y values (see below for details).

### Lesson 1.2 Simulating one subject's data, and likelihoods

# It is essential to be able to simulate data from any model you want to use for
# several reasons:
# 1) So you can play with the model, changing its parameters and seeing the 
#    effects, just as you would run experiments and analyse data to understand
#    the effects of experimental manipulations.
# 2) So you can check your estimation procedure is working. If you can’t recover
#    the parameters that you used to simulate a very large data set data by 
#    fitting to that simulated data then something is wrong.
# 3) To check if your estimation procedure will work in the experimental design
#    you intend to apply it to. Sometimes the estimation procedure works in
#    very large samples but fails in the design you are using in real a real
#    experiment. Simulating and fitting helps guide you to the design you need.

rm(list=ls()); par(mfrow = c(1,1))
source ("dmc/dmc.R")

load_model ("LBA","lba_B.R")

# Load in the model created for the previous lesson
load_data ("dmc_1_1.RData")
# load_data ("dmc_1_2.RData")

# SIMULATING DATA.
#
# The function "simulate.dmc" (in "dmc_model.R") makes a data frame based on the 
# parameter vector and the model, with n observations for each row in model. 
# n can be a scalar or a vector of length equal to the number of rows (i.e., 
# design cells) in model.

data <- simulate.dmc(p.vector,model,n=1)

# Indicative results (subject to randomness):
#    S  R       RT
# 1 s1 r1 3.308407
# 2 s2 r2 2.336722

# The function "data.model.dmc" combines the data frame (which would normally
# come from an experiment rather than being simulated) and the model so it can
# be passed as a single object to "likelihood.dmc" (in "lba_B.R" in dmc/models) 
# in order to calculate the likelihood of a parameter vector.
#
# The likelihood summarizes all of the information in the data relevant to 
# model fitting. As you will see in later lessons having a likelihood is the
# key to being able to perform Bayesian estimation.
#
# By convention, we call a specific combination of experimental data and a 
# particular model a data-model instance (data.model)

data.model <- data.model.dmc(data,model)

# This also adds "cell.index" and "cell.empty" attributes used to speed the 
# computation of the likelihood.

# "likelihood.dmc" computes a vector of likelihoods, one for each RT (in the)
# same order as in the data frame. Note that the output must always be 
# greater than or equal to the function's parameter min.like argument 
# (1e-10 by default).
likelihood.dmc(p.vector,data.model)

# At this point it is worth taking a look at lba_B.R:
#  transform.dmc allows a flexible relationship between the standard model 
#    and the model parameterization (e.g., b vs. B)
#  likelihood.dmc interfaces with the rtdists function

# SIMULATE SOME MORE REALISTIC DATA (with n=1e4 observations per row)
p.vector <- c(A=1,B.r1=1,B.r2=1,mean_v.true=1,mean_v.false=0,sd_v.true=1,t0=.2)
data1 <- simulate.dmc(p.vector,model,1e4)
data.model1 <- data.model.dmc(data1,model)

# Look at the accuracy and density for each stimulus level (correct responses 
# are in black and error responses are in red)
plot.cell.density(data.model1[data.model1$S=="s1",],C="r1",xlim=c(0,4))
plot.cell.density(data.model1[data.model1$S=="s2",],C="r2",xlim=c(0,4))

# Look at how summed log-likelihood varies around the true value using
# function "profile.dmc" in "dmc_plotting.R". This function takes a parameter  
# name and min/max values to graph, and also prints the maximum likelihood  
# estimate for the parameter, which should tend to the true value for large 
# samples
ylim=c(-35000,-32000)
par(mfrow=c(2,3))
profile.dmc("A",           .1,  2,p.vector,data.model1,ylim=ylim)
profile.dmc("B.r1",        .1,  2,p.vector,data.model1,ylim=ylim)
profile.dmc("mean_v.true", .1,  2,p.vector,data.model1,ylim=ylim)
profile.dmc("sd_v.true",   .1,  2,p.vector,data.model1,ylim=ylim)
profile.dmc("t0",          .01,.4,p.vector,data.model1,ylim=ylim)

# NB1: Rstudio sometimes give an error if you try to plot too many panels in 
#      one plot. To avoid this you can sometimes make the plotting area big 
#      before you run the command. You could also ask for fewer panels e.g.,
#      par(mfrow=c(1,3)) which will do one row of 3 panels over two pages, or 
#      you could shrink the margins around the plots, e.g., put mar=rep(2,4)
#      inside par() (this shrinks each of the four margins around a panel to 
#      a value of 2, you can use smaller or larger values).
# NB2: If you dont specify the ylim argument or set its value to NA profile.dmc
#      will choose its own scale, but this will likely differ across 
#      parameters so be careful when you compare graphs for different
#      parameters.

# SIMULATING DATA FOR UNBALANCED DESIGNS.
#
# An unbalanced design can be simulated by making n a vector of the same length
# as the number of cells in the design. Here we have 2 observations for s1 and
# 4 for s2. 
simulate.dmc(p.vector,model,n=c(2,4))

# Indicative results (subject to randomness):
#    S  R        RT
# 1 s1 r2 2.9056123
# 2 s1 r1 3.0826531
# 3 s2 r2 0.8964241
# 4 s2 r2 0.6017682
# 5 s2 r2 1.1289151
# 6 s2 r1 1.0046132

# save_data(data.model,data.model1,file="dmc_1_2.RData")
