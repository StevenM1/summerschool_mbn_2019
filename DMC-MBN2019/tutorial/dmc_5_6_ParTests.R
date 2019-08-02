##################  DMC Lesson 5: Advanced Multiple Subjects


### Lesson 5.6:  Testing parameter effects.
rm(list=ls())
source ("dmc/dmc.R")
load_model ("Wald","wald.R")


# Parameter or posterior inference is a useful complement to model selection,
# where differences are allowed in estimates across conditions then tests are
# made to determine the probability of obtaining a difference of the observed
# magnitude (i.e., a posterior p value). This requires calculating the difference
# for each sampled posterior parameter and calculating the probability that the
# distribution of differences is greater than zero.
# 
# The procedure can be applied to individual subject fits, fixed effects fits
# and hierarchical fits, and in the latter case can also be applied to 
# hyper-parameters. For groups of subjects they can be applied to average 
# differences over subjects (i.e., the distribution if of group average 
# differences or to hyper-parameters). However, in the latter case they will not
# reflect correlations inherit in a with-subject contrast as this structure is
# not presently built in to DMC's hierarchical fitting (it is planned for the
# future). So, it is recommended that tests of within-subject factors are 
# done on averages of individual subject parameters, as they do capture the 
# correlation. However, such tests provide only fixed effects inference (about
# the group of subjects you measured not future subjects). Between subject 
# tests can be done on the hypers.

# We will use some fits of a Wald model with Go Failure (see dmc_6_2_Wald_GF.R)
# as an example, with samples (from a hierarchical fit) contained in an 
# object called "hsamples".
load_data("dmc_5_6.RData")

# The experiment had 20 subjects who had to indicate which of two stimuli 
# occurred (factor S=stimulus: low or high). They could do this either while 
# counting backwards by 3s or not (factor L=load: none vs. 3s). Here is the 
# model that was fit.

modelC <- model.dmc(p.map=list(A="1",B=c("L","R"),v=c("S","L","M"),t0="L",gf="L"),
                    match.map=list(M=list(high="HIGH",low="LOW")),
                    factors=list(L = c("none","3s"),S=c("high","low")),
                    constants=c(A=0),
                    responses=c(high="HIGH",low="LOW"),
                    type="wald")

# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "B.none.HIGH"       "B.3s.HIGH"         "B.none.LOW"       
# [4] "B.3s.LOW"          "v.high.none.true"  "v.low.none.true"  
# [7] "v.high.3s.true"    "v.low.3s.true"     "v.high.none.false"
# [10] "v.low.none.false"  "v.high.3s.false"   "v.low.3s.false"   
# [13] "t0.none"           "t0.3s"             "gf.none"          
# [16] "gf.3s"            


# By default the compare.p function performs a test appropriate for a within-
# subjects effect, averaging over participant differences for two parameters
# specified with the pnames argument (second name - first name). It prints out
# a table of the 95% CI for each component of the difference and the difference
# ("contrast") followed by the probability that the difference is greater than
# zero. The names of the components can be controlled with the "pretty"
# argument. The width of the CI can be controlled by the lo.p (= .025) and 
# hi.p (= .975) arguments and significant digits with the digits argument.
# The test of the t0 parameter show participants were faster for 3s than none 
# by 45ms. 

compare.p(hsamples,pnames=c("t0.none","t0.3s"),pretty=c("none","3s"),show.plot=TRUE)
#        none    3s contrast
# 2.5%  0.234 0.188    0.025
# 50%   0.247 0.202    0.045
# 97.5% 0.261 0.216    0.064
# p.gt.0 
#      1

# In this case the number of posterior samples per parameter is not very large
# (48 chains x 200 iterations = 9600) so it does not take to long to compute 
# using all of them but sometimes you may want to use fewer. This can be done 
# with the subsample argument which allows selection of a subset of interations. 
# For example, suppose you only wanted the first 10 (480, not a great idea ...)
compare.p(hsamples,pnames=c("t0.none","t0.3s"),pretty=c("none","3s"),show.plot=TRUE,subsample=10)

# Alternativley you can pick specific interations
compare.p(hsamples,pnames=c("t0.none","t0.3s"),pretty=c("none","3s"),show.plot=TRUE,subsample=c(190:200))


# Here we do the same test on the first subject by passing only their data
# (note the use of [] not [[]] to do this). We also ask for a plot of the
# posterior difference distribution (note that the table output can be turned
# off with show.table=FALSE). Default line0=TRUE shows a vertical line at zero.
# The legend can be turned off (show.legend=FALSE) and its position controlled
# (lpos="topleft" by defualt).
compare.p(hsamples[1],pnames=c("t0.none","t0.3s"),pretty=c("none","3s"),
          show.plot=TRUE,xlab="Non-decison time (s)")
#        none    3s contrast
# 2.5%  0.072 0.072   -0.091
# 50%   0.122 0.135   -0.013
# 97.5% 0.166 0.195    0.066
# p.gt.0 
#  0.373

# From here we will return to focusing on group inference.

# A function can be used to calculate the difference. Here it first converts
# the gf estimates to the probability scale (they are sampled on the probit)
# scale) before taking the difference. Note the function must return a 3-vector
# containing each element of the difference and the difference itself. The 
# results show 2.1% higher go failure in 3s than none.
compare.p(hsamples,show.plot=TRUE,main="Go Failure",pretty=c("none","3s"),
          fun=function(x){x<-pnorm(x[c("gf.none","gf.3s")]);c(x,-diff(x))},xlab="p(GF)")
#        none    3s contrast
# 2.5%  0.010 0.028   -0.027
# 50%   0.013 0.034   -0.021
# 97.5% 0.016 0.040   -0.014
# p.gt.0 
#      0

# Note that further arguments can be passed to fun through ... arguments. For 
# examples suppose we wished to test if go failure in the 3s condition was 
# greater than some given value.
fun=function(x,gf){x<-c(pnorm(x["gf.3s"]),gf);c(x,-diff(x))}

# We can now make the test passing whatever value of the values (gf) we like.
# Here we compare to a 2% rate.
compare.p(hsamples,pretty=c("3s","Fixed"),fun=fun,gf=.02)
#          3s Fixed contrast
# 2.5%  0.028  0.02    0.008
# 50%   0.034  0.02    0.014
# 97.5% 0.040  0.02    0.020
# p.gt.0 
#      1

# We can make the same test at the hyper level using hpar=TRUE. By defualt this
# this test is done on the location hyper. Estimates are a little different and
# CIs wider as between-subject variance is not partialed out of the contrast.
compare.p(hsamples,hyper=TRUE,
          show.plot=TRUE,main="Go Failure",pretty=c("none","3s"),
          fun=function(x){x<-pnorm(x[c("gf.none","gf.3s")]);c(x,-diff(x))},xlab="p(GF)")
#        none    3s contrast
# 2.5%  0.006 0.010   -0.030
# 50%   0.010 0.021   -0.010
# 97.5% 0.016 0.039    0.002
# p.gt.0 
#  0.053

# The same test can be made on the scale hyper (hpar=1 is the default for 
# testing location), showing greater individual differences in go failures under
# load.
compare.p(hsamples,hyper=TRUE,hpar=2,
          show.plot=TRUE,main="Go Failure",pretty=c("none","3s"),
          fun=function(x){x<-pnorm(x[c("gf.none","gf.3s")]);c(x,-diff(x))},xlab="p(GF)")
#        none    3s contrast
# 2.5%  0.592 0.655   -0.177
# 50%   0.636 0.717   -0.080
# 97.5% 0.701 0.802    0.012
# p.gt.0 
#  0.044

# Returning to test of averages over individual parameters, here we use a 
# function on four parameters, first averaging over thresholds for the
# two accumulators (HIGH and LOW) then taking a difference to see if average
# thresholds differ with load, finding an average B higher in 3s than none.

fun <- function(x){
  none <- mean(x[c("B.none.HIGH","B.none.LOW")])             
  s3 <- mean(x[c("B.3s.LOW","B.3s.HIGH")])  
  c(s3,none,s3-none)
}
compare.p(hsamples,show.plot=TRUE,pretty=c("3s","none"),fun=fun,xlab="Average Threshold")
#          3s  none contrast
# 2.5%  2.125 1.785    0.259
# 50%   2.172 1.841    0.331
# 97.5% 2.215 1.898    0.400
# p.gt.0 
#      1

# This function test the difference between 3s and none in the true-false rate
# differences averaged over high and low stimuli, finding them to be larger in 
# none than 3s.
fun <- function(x){
  hi.none.D <- diff(x[c("v.high.none.false","v.high.none.true")])             
  lo.none.D <- diff(x[c("v.low.none.false","v.low.none.true")])             
  hi.3s.D <- diff(x[c("v.high.3s.false","v.high.3s.true")])             
  lo.3s.D <- diff(x[c("v.low.3s.false","v.low.3s.true")])             
  none <- (hi.none.D+lo.none.D)/2
  s3 <- (hi.3s.D+lo.3s.D)/2
  c(none,s3,none-s3)
}
compare.p(hsamples,show.plot=TRUE,pretty=c("none","3s"),fun=fun,xlab="Average Rate")
#        none    3s contrast
# 2.5%  1.160 0.907    0.160
# 50%   1.233 0.973    0.261
# 97.5% 1.308 1.039    0.362
# p.gt.0 
#      1 

# Sometimes you might want to get all of the subject-average medians and 
# credible intervals (e.g., for plotting). That can be done with as:

ci <- subject.average.ci(fits)
round(ci,2)
#       B.none.HIGH B.3s.HIGH B.none.LOW B.3s.LOW v.high.none.true
# 2.5%         1.71      2.10       1.84     2.13             2.34
# 50%          1.78      2.15       1.90     2.19             2.43
# 97.5%        1.85      2.20       1.97     2.25             2.52
#       v.low.none.true v.high.3s.true v.low.3s.true v.high.none.false
# 2.5%             2.43           2.29          2.26              1.19
# 50%              2.52           2.37          2.34              1.29
# 97.5%            2.61           2.44          2.43              1.38
#       v.low.none.false v.high.3s.false v.low.3s.false t0.none t0.3s
# 2.5%              1.09            1.34           1.25    0.23  0.19
# 50%               1.20            1.44           1.33    0.25  0.20
# 97.5%             1.31            1.52           1.42    0.26  0.22
#       gf.none gf.3s
# 2.5%    -2.43 -2.16
# 50%     -2.32 -2.05
# 97.5%   -2.23 -1.94

# Note the default is to do all parameters (but you can do a subset by 
# specifying their names using the pnames argument) and the 95% CI (but you 
# can set lo.p and hi.p to get any interval).

# You can also calculate functions of parametrs. For example, suppose you 
# wanted to plot the difference in the true and false rates. 
fun <- function(x) {
  c(diff(x[c("v.high.none.false","v.high.none.true")]),
    diff(x[c("v.low.none.false","v.low.none.true")]),
    diff(x[c("v.high.3s.false","v.high.3s.true")]),
    diff(x[c("v.low.3s.false","v.low.3s.true")]))
}
# You can use the pretty argument to give the differneces names
pretty <-  c("dv.high.none","dv.low.none","dv.high.3s","dv.low.3s")

ci.dv <- subject.average.ci(hsamples,fun=fun,pretty=pretty)
round(ci.dv,2)
#       dv.high.none dv.low.none dv.high.3s dv.low.3s
# 2.5%          1.01        1.18       0.81      0.89
# 50%           1.15        1.32       0.93      1.01
# 97.5%         1.28        1.46       1.05      1.13


### BETWEEN SUBJECT TESTS

# Here we illustrate testing between subject effects by copying the hsamples
# object as if it were two groups of subjects

hA1 <- hsamples
hA2 <- hsamples
# We first put the two fit objects into a list. If the list has names these
# will be used in output comparisons, otherwise they will be called G1 and G2.
fits=list(A1=hA1,A2=hA2)

# The first test compares the t0 parameter, specified using the within.pnames
# argument, at the hyper level. The arguments controlling outputs for the
# within subject case by compare.p, like "main" in this example, are also used 
# in the between subject version. Of course as they are identical copies we
# expet no difference.
compare.ps(fits=fits,within.pnames="t0.3s",main="t0.3s",hyper=TRUE)
#          A1    A2 contrast
# 2.5%  0.068 0.068        0
# 50%   0.175 0.175        0
# 97.5% 0.232 0.232        0
# p.gt.0 
#      0
     
# As shown the output is also the same as the within subjects version. The test
# can also be done on averages of participant level parameters rather than 
# hyper parameters, but here there is no advantage as there are no within-subject
# correlations to be take account of. As expected the CIs are much tighter for 
# each condition.
compare.ps(fits=fits,within.pnames="t0.3s",main="t0",hyper=FALSE)
#          A1    A2 contrast
# 2.5%  0.188 0.188        0
# 50%   0.202 0.202        0
# 97.5% 0.216 0.216        0
# p.gt.0 
#      0
     

# From here we focus on the hyper parameters, which is the default setting. The
# next example shows how to use a function to combine parameters within each 
# group as input to the between subjects test using the within.fun argument. 
# Note that this function must return a single number, here the average of the
# two t0 parameters.
compare.ps(fits=fits,main="t0",
  within.fun=function(x){mean(x[c("t0.3s","t0.none")])})
#          A1    A2 contrast
# 2.5%  0.155 0.155        0
# 50%   0.209 0.209        0
# 97.5% 0.243 0.243        0
# p.gt.0 
#      0
     
