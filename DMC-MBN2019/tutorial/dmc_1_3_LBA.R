##################  DMC Lesson 1 Models and Data
#
# Note: Before running these scripts, the current working directory must be set 
# to the top-level folder containing the dmc and tutorial subfolders 

### Lesson 1.3 Exploring the LBA model

# This lesson is about familiarising yourself with the LBA model and its 
# parameters, and how you can create different scenarios with different 
# parameter values. See Brown, S. D., & Heathcote, A. (2008). The simplest 
# complete model of choice response time: Linear ballistic accumulation. 
# Cognitive Psychology, 57(3), 153–178. http://doi.org/10.1016/j.cogpsych.2007.12.002



rm(list=ls()); par(mfrow = c(1,1))
source ("dmc/dmc.R")

load_model ("LBA","lba_B.R")

# NB: For this model definiton, and most others supplied, t0 cannot differ
#     between accmulators (if it does the value for the first accumulator is
#     used for all and the others ignored). However, the lba_Bt0.R model 
#     definiton relaxes this restriciton. In this case any difference must be
#     interprited as in the response production part of t0 (as any difference
#     in encoding time will effect which accumulator wins the race, whereas any
#     difference in response production time is simply accomodated by subtracting
#     differnet values from observed RTs before calcualting the likelihood for
#     for each accumulator or adding differnet values when simulating).

# Start with a simple design with a binary stimulus factor, S 
# (e.g, left vs. right motion) and corresponding responses.  
factors=list(S=c("left","right"))
responses=c("LEFT","RIGHT")
match.map=list(M=list(left="LEFT",right="RIGHT"))

# NB1: For the LBA one of B, mean_v or sd_v must be fixed in at least one
#     cell of the design. In this case we fix sd_v.false = 1. For simplicity we 
#     will also fix mean_v.false=0.
#
# NB2: See ?LBA for rtdists help on parameters
consts=c(st0=0,sd_v.false=1,mean_v.false=0)

# Design where there are no effects except greater than chance accuracy due to
# a mean_v.true > mean_v.false
p.map=list(A="1",B="1",mean_v="M",sd_v="M",t0="1",st0="1")
model <- model.dmc(type="norm",constants=consts,p.map=p.map,
                   match.map=match.map,factors=factors,responses=responses)

# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"           "B"           "mean_v.true" "sd_v.true"   "t0"         

# 1) Fast and error prone performance
#    N.B. sd_v.true < sd_v.false as is often seen in practice
p.vector  <- c(A=1,B=0,mean_v.true=1,sd_v.true=0.66,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)

# Score & plot
correct <- data.model$S==tolower(data.model$R)
round(mean(correct),2)                       # Accuracy
round(tapply(data.model$RT,list(correct),mean),2) # Errors/correct similar speed

plot.cell.density(data.model,C=correct,xlim=c(0,5)) # (red represents error)

# 2) Raise the threshold to improve accuracy, makes errors slow
p.vector  <- c(A=1,B=1,mean_v.true=1,sd_v.true=0.66,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
# Score & plot (using a convenience function which does the same as 
# the previous code, useful for simple scoring situations such as this)
plot.score.dmc (data.model)

# 3) Lowering start-point noise has similar effects
p.vector  <- c(A=.25,B=.75,mean_v.true=1,sd_v.true=0.66,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
plot.score.dmc (data.model)

# Further analyses are relative to case (2) above

# 4) Increasing matching rate increases accuracy and speed
p.vector  <- c(A=1,B=1,mean_v.true=1.5,sd_v.true=0.66,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
plot.score.dmc (data.model)

# 5) Decreasing sd_v.true increases accuracy and creates fast errors
p.vector  <- c(A=1,B=1,mean_v.true=1.5,sd_v.true=0.33,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
plot.score.dmc (data.model)

# 6) Increasing t0 slows things down again but has no effect on accuracy
p.vector  <- c(A=1,B=1,mean_v.true=1.5,sd_v.true=0.66,t0=.4)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
plot.score.dmc (data.model)

# NB: t0 has the same effect for all models.

# Let's allow non-decision time to be variable across trials, so do not set st0=0.
# Let's re-run (2) for comparison (but with a larger t0), and add a
# printout of variability (as that is what's affected by st0). We use the inter-
# quartile range to measure variability as it is robust against large outliers.
p.vector  <- c(A=1,B=1,mean_v.true=1,sd_v.true=0.66,t0=.4)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
# Score & plot 
plot.score.dmc (data.model, IQR=TRUE)  # Show inter-quartile range

# 7) We set st0 to 0.2s. As t0 specifies the MINIMUM non-decision time we
#    subtract st0/2 from t0 in order to maintain the same mean non-decision 
#    time (what is usually referred to as ter). Also, any change in non-decision 
#    time parameters has no effect on accuracy (because t0 distributes 
#    uniformly.)
consts=c(sd_v.false=1,mean_v.false=0); consts
p.map=list(A="1",B="1",mean_v="M",sd_v="M",t0="1",st0="1"); p.map
model <- model.dmc(type="norm",constants=consts,p.map=p.map,
                   match.map=match.map,factors=factors,responses=responses)
p.vector  <- c(A=1,B=1,mean_v.true=1,sd_v.true=0.66,t0=.1,st0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
# Score & plot
plot.score.dmc (data.model, IQR=TRUE)

# 8) Response bias can be introduced by allowing B to vary with accumulator (R)
consts=c(st0=0,sd_v.false=1,mean_v.false=0)
p.map=list(A="1",B="R",mean_v="M",sd_v="M",t0="1",st0="1")
model <- model.dmc(type="norm",constants=consts,
                   p.map=p.map,match.map=match.map,
                   factors=factors,responses=responses)
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"           "B.LEFT"      "B.RIGHT"     "mean_v.true" "sd_v.true"   "t0"         

# Bias to respond left
p.vector  <- c(A=1,B.LEFT=.5,B.RIGHT=1.5,mean_v.true=1,sd_v.true=0.66,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)

correct <- data.model$S==tolower(data.model$R)
# Right is now below chance
round(tapply(correct,data.model$S,mean),2)      
# Correct responses to left are faster and to right slower, but the opposite
# patterns shows in errors (and hence slow errors for left and fast for right)
round(tapply(data.model$RT,list(data.model$S,correct),mean),2)   
par(mfrow=c(1,2))
plot.cell.density(data.model[data.model$S=="left",],
                  C=correct[data.model$S=="left"],
                  xlim=c(0,5),ymax=1.0,main="left")
plot.cell.density(data.model[data.model$S=="right",],
                  C=correct[data.model$S=="right"],
                  xlim=c(0,5),ymax=1.0,main="right")
