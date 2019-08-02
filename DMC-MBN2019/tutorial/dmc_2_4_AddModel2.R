##################  DMC Lesson 2: Priors, Posteriors and Adding New Models
#
# THIS IS AN ADVANCED LESSON, BEST DONE AFTER WORKING THROUGH TO THE END OF LESSON 4
#
# Note: Before running these scripts, the current working directory must be set 
# to the top-level folder containing the dmc and tutorial subfolders 
#
#
### Lesson 2.4:  Advanced model building

rm(list=ls())
source ("dmc/dmc.R")

### EXAMPLE 1: Random threshold levels

# First we will look at an LBA model type where there are two different threshold 
# levels with either one or the other being used randomly from trial to trial.
# A lesson about this model, lbaBp2_lesson.R is available in the LBA subdirectory.

load_model("LBA","lba_Bp2.R")

# This model has two b parameters, b1 and b2 (in order of increasing size). We 
# first have to modify the random function for the standard LBA to sample n2 out
# of n trials using b2 and the remainder using b1, where n2 is drawn from a 
# binomial with probability pb2.

random.dmc <- function(n,p.df,model)
{
  n2 <- rbinom(1,n,p.df$pb2[1])
  out <- matrix(nrow=n,ncol=2)
  dimnames(out) <- list(NULL,c("rt","response"))
  if (n2<n) out[1:(n-n2),] <- rlba.norm(n-n2,A=p.df$A,b=p.df$b1,t0=p.df$t0[1],
                                        mean_v=p.df$mean_v,sd_v=p.df$sd_v,st0=p.df$st0[1],
                                        posdrift=attr(model,"posdrift"))
  if (n2>0) out[(n-n2+1):n,] <- rlba.norm(n2,A=p.df$A,b=p.df$b2,t0=p.df$t0[1],
                                          mean_v=p.df$mean_v,sd_v=p.df$sd_v,st0=p.df$st0[1],
                                          posdrift = attr(model,"posdrift"))
  data.frame(out)
}


# In likelihood.dmc we calculate the overall likelihood = probability of using 
# b1 x likelihood when using b1 plus probability of using b2 x likelihood when 
# using b2. 

likelihood.dmc <- function(p.vector,data,min.like=1e-10) 
  # Returns vector of likelihoods for each RT in data (in same order)
{
  
  if (attr(attr(data,"model"),"type")=="normBp2") dist <- "norm" else dist <- "lnorm"
  likelihood <- numeric(dim(data)[1])
  for ( i in row.names(attr(data,"model")) ) if ( !attr(data,"cell.empty")[i] )
  {
    p.df <- p.df.dmc(p.vector,i,attributes(data)$model,n1order=TRUE)
    likelihood[ attr(data,"cell.index")[[i]] ] <-
      (1-p.df$pb2[1])*n1PDF(data$RT[attr(data,"cell.index")[[i]]],
                            A=as.list(p.df$A),b=as.list(p.df$b1),
                            mean_v=p.df$mean_v,sd_v=p.df$sd_v,
                            t0=p.df$t0[1],st0=p.df$st0[1], 
                            distribution=dist,
                            args.list = list(posdrift=attr(attr(data,"model"),"posdrift")),
                            silent=TRUE) + 
      p.df$pb2[1]*n1PDF(data$RT[attr(data,"cell.index")[[i]]],
                        A=as.list(p.df$A),b=as.list(p.df$b2),
                        mean_v=p.df$mean_v,sd_v=p.df$sd_v,
                        t0=p.df$t0[1],st0=p.df$st0[1], 
                        distribution=dist,
                        args.list = list(posdrift=attr(attr(data,"model"),"posdrift")),
                        silent=TRUE) 
  }
  
  pmax(likelihood,min.like)
}


# Also transform.R had to be modified to map from the external 
# parameters A,B1,B2,pb2,mean_v,sd_v,t0, and st0 to the internal parameters A,
# b1,b2,pb2,mean_v,sd_v,t0,st0, and to make b1<b2 out of B1 (the gap from the
# top of A to b1) and B2 (the gap from b1 to b2). 

transform.dmc <- function(par.df) 
{
  # User supplied tranforms go here
  par.df$b1 <- par.df$B1+par.df$A
  par.df$b2 <- par.df$B2+par.df$b1
  
  par.df[,c("A","b1","b2","pb2","t0","mean_v","sd_v","st0")]
}

### EXAMPLE 2:  Different number of accumulators

# Next, we will look at an LBA model where there can be different numbers of 
# accumulators in different conditions. We have a factor F with level f1 where 
# participants perform a two-choice task and another f2 where they perform a 
# three choice task. In the former task they must choose among two stimuli 
# (s1 and s2) and in the latter three (also s3).

load_model ("LBA","lbaN_B.R")

# We define a "parameter" N whose sole purpose is to carry information about the
# number of accumulators into likelihood.dmc, with its values being
# set as constants of either 2 or 3 as appropriate. 

model <- model.dmc(
  p.map=list(A="1",B="1",mean_v=c("F","M"),sd_v=c("1"),t0="1",st0="1",N=c("F")),
  constants=c(st0=0,sd_v=1,N.f1=2,N.f2=3),
  match.map=list(M=list(s1="r1",s2="r2",s3="r3")),
  factors=list(S=c("s1","s2","s3"),F=c("f1","f2")),
  responses=c("r1","r2","r3"),type="norm")

# In lbaN_B.R we use N to modify how long the parameter vectors are, which in 
# turn modifies the number of parameters in the race, in both random.dmc

p.vector=c(A=0,B=1,mean_v.f1.true=2,mean_v.f2.true=3,
           mean_v.f1.false=4,mean_v.f2.false=5,t0=6)

# Let's take look at the p.df data frame that random.dmc acts upon that is
# created by a call in to p.df.dmc in simulate.dmc for a particular cell. 
# Here we show the first cell (where N=2) and the fourth (where N=3)

print.cell.p(p.vector, model)
# 
# [1] "s1.f1.r1"
#    A b t0 mean_v sd_v st0 N
# r1 0 1  6      2    1   0 2
# r2 0 1  6      4    1   0 2
# r3 0 1  6      4    1   0 2
# 
# ...
# 
# [1] "s1.f2.r1"
#    A b t0 mean_v sd_v st0 N
# r1 0 1  6      3    1   0 3
# r2 0 1  6      5    1   0 3
# r3 0 1  6      5    1   0 3
# 
# ...


# Lets focus on the first cell where N=2
p.df <- p.df.dmc(p.vector,row.names(model)[1],model,n1order=FALSE)

# In random.dmc we need to remove the final row, which makes
# the race between two accumulators rather than three:
p.df[1:p.df$N[1],]
#    A b t0 mean_v sd_v st0 N
# r1 0 1  6      2    1   0 2
# r2 0 1  6      4    1   0 2

random.dmc<- function(n,p.df,model)
{
  rlba.norm(n,
            A=p.df$A[1:p.df$N[1]],
            b=p.df$b[1:p.df$N[1]],
            mean_v=p.df$mean_v[1:p.df$N[1]],
            sd_v=p.df$sd_v[1:p.df$N[1]],
            t0=p.df$t0[1], 
            st0=p.df$st0[1],
            posdrift = attr(model,"posdrift"))
}

# and likelihood.dmc 

likelihood.dmc <- function(p.vector,data,min.like=1e-10) 
  # Returns vector of likelihoods for each RT in data (in same order)
{
  
  likelihood <- numeric(dim(data)[1])
  for ( i in row.names(attr(data,"model")) ) if ( !attr(data,"cell.empty")[i] )
  {
    p.df <- p.df.dmc(p.vector,i,attributes(data)$model,n1order=TRUE)
    likelihood[ attr(data,"cell.index")[[i]] ] <-
      n1PDF(data$RT[attr(data,"cell.index")[[i]]],
            A=as.list(p.df$A[1:p.df$N[1]]),
            b=as.list(p.df$b[1:p.df$N[1]]),
            mean_v=p.df$mean_v[1:p.df$N[1]],
            sd_v=p.df$sd_v[1:p.df$N[1]],
            t0=p.df$t0[1], 
            st0=p.df$st0[1],
            distribution="norm",
            args.list = list(posdrift=attr(attr(data,"model"),"posdrift")),
            silent=TRUE)
    
  }
  pmax(likelihood,min.like)
}

# Note that only the accumulator parameters are set with N (ie., not t0 and st0)
# and that only the first value of p.df$N.

# We can see this has the desired effect by simulating some data.
p.vector  <- c(A=.25,B=1,mean_v.f1.true=1,mean_v.f2.true=2,
               mean_v.f1.false=.25,mean_v.f2.false=.5,t0=.2)
data <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)

# Here we can see that r3 is never made for F=f1. 
table(data[,c("S","R","F")])
# , , F = f1
# 
#     R
# S      r1   r2   r3
#   s1 6521 3479    0
#   s2 3514 6486    0
#   s3 5016 4984    0
# 
# , , F = f2
# 
#     R
# S      r1   r2   r3
#   s1 7161 1415 1424
#   s2 1432 7139 1429
#   s3 1553 1517 6930

# Complications can arise when the S factor comes into the model and so there
# are some parameters that should not be sampled when there are fewer accumulators.
# Consider the following model, where we need to set parameters in F=f1 for s3
# to a constant as s3 is never presented in f1. The value of the constant has no
# effect. 

model <- model.dmc(
  p.map=list(A="1",B="1",mean_v=c("F","S","M"),sd_v=c("1"),t0="1",st0="1",N=c("F")),
  constants=c(st0=0,sd_v=1,N.f1=2,N.f2=3,mean_v.f1.s3.true=0,mean_v.f1.s3.false=0),
  match.map=list(M=list(s1="r1",s2="r2",s3="r3")),
  factors=list(S=c("s1","s2","s3"),F=c("f1","f2")),responses=c("r1","r2","r3"),type="norm")

### EXAMPLE 3: NO match.map

# In all of the models examined so far we have used a "match.map" in order to 
# construct the match factor, M, which maps rates to matching and mismatching
# accumulators. However in some situations there is no "correct" response and so
# we want the flexibility to do this in an arbitrary way. The following example
# shows how to do this by creating a factor with level names which can be parsed
# in transform.dmc (by a special purpose function you create there) to set the
# values of parameters for each accumulator (in our example two accumulators 
# called r1 and r2) as required.  

# Here we do this for mean_v in a complicated example that shows the great 
# flexibility of this approach (which comes at the cost of needing some tricky
# functions). In the example each cell of the 6 cells of the design has an 
# associated time and distance covariate for each accumulator (i.e., 4 values
# per cell). The name of each cell has 4 elements that designate these 4 values.
# Rather than just pasting the actual covariate values together (which could get
# messy if they are real numbers with many decimals) we instead paste either the
# numbers 1 or 2, which then serve to lookup the associated covariate values 
# from a table we create in transform.dmc

load_model ("LBA","lba_B_T^-1_D.R")

# First create the factor, which we call DT (distance and time) 
DT <- c("1112","1121","1122","1221","1222","2121")
# Each digit in these levels is used to lookup a numeric value for covariates 
# used to calculate mean_v for each accumulator in the cell. In this case the
# the first two digits (1000s and 100s) refer to T(ime) and D(istance)
# covariates for accumulator r1, and the last two (10s and units) the same for
# accumulator r2.

# Now create a vector to be used as the constants argument.
# First put the levels in as values then make the names appropriate for the model
# that it will be used to create.
consts <- as.numeric(DT); names(consts) <- paste("mean_v",DT,sep=".")
consts
# mean_v.1112 mean_v.1121 mean_v.1122 mean_v.1221 mean_v.1222 mean_v.2121 
#        1112        1121        1122        1221        1222        2121 

# Append on the usual constants to set the scale and turn off non-decision time 
# variability
consts <- c(consts,c(sd_v=1,st0=0))

# Note that no match.map is specified as mapping to the appropriate 
# accumulator (A or B) handled in the model file lba_B_T^-1_D.R. 
model <- model.dmc(p.map = list(A="1",B="R",t0="1",sd_v="1",st0="1",
                                Wt="1",Wd="1",mean_v="DT"), 
                   factors=list(DT=DT),
                   constants = consts, 
                   responses = c("r1","r2"),
                   type="norm")
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"    "B.r1" "B.r2" "t0"   "Wt"   "Wd"  
# 
# Constants are (see attr(,"constants") ):
# mean_v.1112 mean_v.1121 mean_v.1122 mean_v.1221 mean_v.1222 mean_v.2121        sd_v         st0 
#        1112        1121        1122        1221        1222        2121           1           0 

# To see how the mapping is done here is the state of par.df as it is passed into
# transform.dmc for design cell "1112.r1" (i.e., cell 1112 for data where 
# response 1 was made). Note that mean_v contains the numeric values that will
# be used to perform the lookup. The entries in mean_v are redundant over 
# accumulators so only the first row needs to be used.

# To do this we need to fake up some things as if we were inside DMC. This stuff
# is pretty esoteric so you can ignore it and just take on faith it gives us 
# par.df in the right state.

# First create a parameter vector
p.vector <- attr(model,"p.vector")
p.vector[1:length(p.vector)] <- c(1,1,2,.2,1,.1)
# Now run par.df but stop transform.dmc from doing anything
real.transform.dmc <- transform.dmc 
transform.dmc <- function(par.df) par.df
par.df <- p.df.dmc(p.vector,"1112.r1",model,n1order=TRUE)
par.df
#    A B  t0 sd_v st0 Wt  Wd mean_v
# r1 1 1 0.2    1   0  1 0.1   1112
# r2 1 2 0.2    1   0  1 0.1   1112

# Now put the real transform function back.
transform.dmc <- real.transform.dmc


# We will look at transform.dmc from lba_B_T^-1_D.R in pieces. First the head 

# transform.dmc <- function(par.df,
#   AT=c('1'=1,'2'=2),   BT=c('1'=1,'2'=2),   # Could just code T and D as 
#   AD=c('1'=10,'2'=20), BD=c('1'=10,'2'=20)) # A and B assumed the same.

# In this example we have put the lookup table r1 (AT and AD) and r2 (BT and BD)
# (as noted these are identical so we could just use the same for both r1 and r2).
# For example AT=c('1'=1,'2'=2) means that AT["1"] will return a time of 1 and
# AT["2"] will return a time of 2. Similarly AD["1"] will return a distance of 
# 10 and AD["2"] a distance of 20.

# Lets make the default arguments available.packages
AT=c('1'=1,'2'=2)
BT=c('1'=1,'2'=2)
AD=c('1'=10,'2'=20)
BD=c('1'=10,'2'=20)

# The first line of the function splits up the the first entry of mean_v into
# its 4 parts, where 1000s map to AT, 100s to AD, 10s to BT and 1s to BD

lookup <- strsplit(as.character(par.df["r1","mean_v"]),"")[[1]]

# In this case lookup has the values [1] "1" "1" "1" "2"

# Next these values are used to set the mean_v parameter for the first 
# accumulator, which is the sum of 1) the inverse of the product of the
# time covariate AT[lookup[1]]=1 and the time weight parameter par.df["r1","Wt"]
# and 2) the product of the distance covariate AD[lookup[2]]=10 and the 
# distance weight parameter par.df["r1","Wd"]. 
par.df["r1","mean_v"] <- 1/(par.df["r1","Wt"]*AT[lookup[1]]) +
  (par.df["r1","Wd"]*AD[lookup[2]])

# The next line set mean_v for r2 accumulator in the same way 
par.df["r2","mean_v"] <- 1/(par.df["r2","Wt"]*BT[lookup[3]]) +
  (par.df["r2","Wd"]*BD[lookup[4]])

# The remainder of the transform.dmc does the standard things for an lba_B
# type model. 

par.df$b <- par.df$B+par.df$A
par.df[,c("A","b","t0","mean_v","sd_v","st0")]
#    A b  t0 mean_v sd_v st0
# r1 1 2 0.2      2    1   0
# r2 1 3 0.2      3    1   0

# random.dmc and likelihood.dmc are similarly standard. A lesson
# on this model including a parameter recovery study called lba_B_T^-1_D_lesson.R
# is contained in the LBA directory. It provides a model of the sort of checking
# that it is important to do when developing a new model. We close with an 
# important word of caution. Even when you get all of the technical details of
# a new model right it is important to check whether it is possible to recover
# its parameters IN THE SAME SORT OF DESIGN AS YOU WILL APPLY IT TO. Even when
# recovery is good with large samples that is no guarantee of adequate 
# performance in smaller samples.


# EXAMPLE 4:  PDA (Probability Density Approximation)

# If you don’t have an analytic likelihood but can simulate data from a model
# you can use a density estimator to get an approximate likelihood, see:
# 
# Turner, B. M., & Sederberg, P. B. (2014). A generalized, likelihood-free method 
# for posterior estimation. Psychonomic Bulletin & Review, 21(2), 227–250. 
# http://doi.org/10.3758/s13423-013-0530-0
# 
# Holmes, W. R. (2015). A practical guide to the Probability Density Approximation 
# (PDA) with improved implementation and error characterization. Journal of 
# Mathematical Psychology, 68-69, 13–24. http://doi.org/10.1016/j.jmp.2015.08.006
#
# This was used to fit an LBA to data where the stimulus changed from supporting
# one response to supporting the other (the PLBA) in:
#
# Holmes, W. R., Trueblood, J. S., & Heathcote, A. (2016). A new framework for 
# modelling decisions about changing information: The Piecewise Linear Ballistic 
# Accumulator model. Cognitive Psychology, 85, 1–29. 
# http://doi.org/10.1016/j.cogpsych.2015.11.002
#
# Numerical integration can be used to get the PLBA likelihood but PDA ended up
# being quicker (depending on how accurately you approximate the likelihood,
# which depends on how many samples you simulate). We will use the PLBA with
# PDA as an example 

load_model("LBA-pda","plba_B.R")

# Note there is also an lba_B.R model enabling you to fit the standard LBA with
# pda so you can compare its performance with a case with an analytic likelihood.

# dists.R in the models/LBA-pda directory defines a random function for the
# PLBA. 
# 1) The A, B, v (= mean_v), sv (= sd_v) and t0 parameter are as for the LBA. 
# swt is the time at which the stimulus switches.
# 2) At a delay rD after the switch the rates can change to w and sw. 
# 3) At a delay tD the threshold can change from A+B to A+B+C
# 
# Here we all the PLBA random function where the switch happens at 0.5s then 
# after a 0.3s delay the rate and threshold switch. Note that the rate parameters
# are defined assuming the first direction is "correct", i.e., v has a larger 
# rate in the first position, and w (after the switch) has a larger rate in the
# second position.

par.df <- cbind.data.frame(A=c(1.5,1.5),B=c(1.2,1.2),C=c(0.3,0.1),
                           v=c(3.32,2.24),w=c(1.51,3.69),sv=c(1,1),sw=c(1,1),
                           rD=0.3,tD=.3,swt=0.5,t0=0.08)

# Lets simulate and plot some data
srts <- random.dmc(1e5,par.df)

# A plot of the simulated data shows two modes, one for response before the
# stimulus switch (where choice 1 is most common), one after (where choice 2
# is most common).
x=c(1:200)/100
d1 <- pda.dmc(x,srts[srts[,"choice"]==1,"rt"],n=table(srts[,"choice"])[1])
d2 <- pda.dmc(x,srts[srts[,"choice"]==2,"rt"],n=table(srts[,"choice"])[2])
plot(x,d1,type="l")
lines(x,d2,col="red")

# The function pda.dmc is actually the one used to calculate the approximate 
# likelihood from the simulated samples. It is just a wrapper for the standard
# R "density" function to get the standard 512 density estimates at regular
# intervals (i.e., for the simulated samples srts[srts[,"choice"]==1,"rt"] and
# srts[srts[,"choice"]==2,"rt"]). The standard "bw.nrd0" plugin bandwidth
# estimate is used, and can be adjusted with this function's adjust argument, 
# which is set at 0.8 by default (i.e., bw = 80% of the plugin estimate as 
# recommended by Holmes). The intervals are take over the range from the minimum
# value where a likelihood is to be calculated (i.e., min(x) above) minus 3*bw
# to 3*bw above the maximum. Finally the approx function is used to look up the 
# required likelihoods (i.e., at each value of x) by linear interpolation 
# between the 512 estimates produced by density. There is also some protection
# against numerical issues in pda.dmc


# Now lets build a standard DMC model. The response is either a RIGHT or LEFT
# button press. The stimulus is either RL (right motion initially then left 
# motion) or LR (left then right motion). We allow v and w to vary with M(atch)
# so we can specify the true and false rates before (v) and after (w) the switch.
# Note that match.map is defined relative to the first direction. 
# In transform.dmc the order of the elements of the w and sw 2-vectors are 
# flipped so effectively w and sw use the complimentary match.map (i.e., 
# match.map=list(M=list(LR="RIGHT",RL="LEFT"))). Also note the switch time (swt)
# is given as a constant.
model <- model.dmc(
  p.map = list(A="1",B="1",C="1",v="M",w="M",sv="1",sw="1",
               rD="1",tD="1",t0="1",swt="1"), 
  match.map=list(M=list(RL="RIGHT",LR="LEFT")),
  factors=list(S=c("LR","RL")),
  constants=c(swt=.5),
  responses = c("LEFT","RIGHT"),
  type="norm")
# Parameter vector names are: ( see attr(,"p.vector") )
#  [1] "A"       "B"       "C"       "v.true"  "v.false" "w.true"  "w.false"
#  [8] "sv"      "sw"      "rD"      "tD"      "t0"     
# 
# Constants are (see attr(,"constants") ):
# swt 
# 0.5 
# 
# Model type = norm (posdrift= TRUE ) 

# Lets simulate some data where the probability of responding correctly to the
# first stimulus is about 75%.
p.vector=c(A=1.5,B=1.2,C=0.3,v.true=3,v.false=1,w.true=3,w.false=1,
           sv=1,sw=1,rD=.3,tD=.3,t0=.08)

raw.data <- simulate.dmc(p.vector,model,1e4) # 20,000 data points
data.model <- data.model.dmc(raw.data,model,n.pda=1e4)
par(mfrow=c(1,2))
plot.cell.density(data.cell=data.model[data.model$S=="LR",],xlim=c(0,3))
plot.cell.density(data.cell=data.model[data.model$S=="RL",],xlim=c(0,3))

# When making the data.model object we also specify how many samples should be
# taken when using pda to calculate the likelihood (n.pda).
data.model <- data.model.dmc(raw.data,model,n.pda=1e4)

# In this case it takes around a quarter of a second to calculate the 20,000
# likelihoods, i.e., to simulate two samples 1e4 in size then do the 
# interpolation. 
system.time({
  print(head(likelihood.dmc(p.vector,data=data.model,min.like=1e-10)))
})
#    user  system elapsed 
#   0.288   0.000   0.286 

# Note that likelihood.dmc has some special accommodations. In particular it
# does not loop over every design cell and every response, as is standard, but 
# instead has an outer loop for design cells, at which point it simulates, then
# for each response runs pda.dmc to get the corresponding likelihoods.

# The following code use profile plots to check the likelihoods. You can see 
# that the likelihood is a little rough, especially for the parameters that
# rely mostly on things that happen after the switch (as these are rarer). It
# runs reasonably fast; each profile plot has 100 points, so in total the 
# likelihood of the 20,000 data points is evaluated 1200 times in a little
# under 3 minutes.

lo <- p.vector*.25
hi <- p.vector*2
par(mfrow=c(2,6))
system.time({
  A=profile.dmc("A",lo["A"],hi["A"],p.vector,data.model)
  B=profile.dmc("B",lo["B"],hi["B"],p.vector,data.model)
  C=profile.dmc("C",lo["C"],hi["C"],p.vector,data.model)
  v.true=profile.dmc("v.true",lo["v.true"],hi["v.true"],p.vector,data.model)
  v.false=profile.dmc("v.false",lo["v.false"],hi["v.false"],p.vector,data.model)
  w.true=profile.dmc("w.true",lo["w.true"],hi["w.true"],p.vector,data.model)
  w.false=profile.dmc("w.false",lo["w.false"],hi["w.false"],p.vector,data.model)
  sv=profile.dmc("sv",lo["sv"],hi["sv"],p.vector,data.model)
  sw=profile.dmc("sw",lo["sw"],hi["sw"],p.vector,data.model)
  rD=profile.dmc("rD",lo["rD"],hi["rD"],p.vector,data.model)
  tD=profile.dmc("tD",lo["tD"],hi["tD"],p.vector,data.model)
  t0=profile.dmc("t0",lo["t0"],hi["t0"],p.vector,data.model)
})
#    user  system elapsed 
# 163.266   0.237 164.672
est.1e4 <- c(A,B,C,v.true,v.false,w.true,w.false,sv,sw,rD,tD,t0)
print(est.1e4)
#      1.46      1.21      0.32      2.98      0.85      2.98      0.71      0.94 
# -9313.404 -9300.378 -9255.546 -9274.504 -9283.925 -9280.560 -9277.178 -9268.905 
#      1.06       0.3      0.29      0.08 
# -9294.667 -9289.425 -9262.951 -9305.715

# Generally the peaks are not to far off the true values.
p.vector-as.numeric(names(est.1e4))
#       A       B       C  v.true v.false  w.true w.false      sv      sw      rD 
#    0.04   -0.01   -0.02    0.02    0.15    0.02    0.29    0.06   -0.06    0.00 
#      tD      t0 
#    0.01    0.00 

