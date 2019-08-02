##################  DMC Lesson 1 Models and Data
#
# Note: Before running these scripts, the current working directory must be set 
# to the top-level folder containing the dmc and tutorial subfolders 

### Lesson 1.4 Exploring the LNR model

# This lesson is about familiarising yourself with the LNR model and its 
# parameters, and how you can create different scenarios with different 
# parameter values. See Heathcote, A., & Love, J. (2012). Linear deterministic 
# accumulator models of simple choice. Frontiers in Psychology, 3. 
# http://doi.org/10.3389/fpsyg.2012.00292/abstract


rm(list=ls()); par(mfrow = c(1,1))
source ("dmc/dmc.R")

load_model ("LNR","lnr.R")

# Start with a simple design with a binary stimulus factor, S 
# (e.g, left vs. right motion) and corresponding responses.  
factors=list(S=c("left","right"))
responses=c("LEFT","RIGHT")
match.map=list(M=list(left="LEFT",right="RIGHT"))

# NB: The LNR model has parameters meanlog, sdlog (see ?plnorm) and t0. 
#     The LNR model is not implemented in the rtdists package but 
#     "rtdists_extras.R" provides the necessary functions (rlnr,   
#     a random function, and n1PDF.lnr). 

# Only a M(atch) effect
model <- model.dmc(type="lnr",constants=c(st0=0),
                   p.map=list(meanlog="M",sdlog="M",t0="1",st0="1"),
                   match.map=list(M=list(left="LEFT",right="RIGHT")),
                   factors=list(S=c("left","right")),
                   responses=c("LEFT","RIGHT"))
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "meanlog.true"  "meanlog.false" "sdlog.true"    "sdlog.false"   "t0"           

# 1) In this example accuracy turns out to be around 75%
p.vector  <- c(meanlog.true=-1,meanlog.false=0,
               sdlog.true=1,sdlog.false=1,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
plot.score.dmc (data.model)

# In the LNR the threshold and mean rate effects both cause changes in the 
# meanlog parameters. In particular an increase in the threshold causes a
# increase in meanlog and an increase in rate causes a decrease in meanlog 
# (i.e., rates and thresholds operate in opposite directions).

# 2) An increase in meanlog.true decreases accuracy (as it gets closer to 
#    meanlog.false) and speed
p.vector  <- c(meanlog.true=-.5,meanlog.false=0,
               sdlog.true=1,sdlog.false=1,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
plot.score.dmc (data.model)

# 3) A decrease in meanlog.true increases accuracy and speed
p.vector  <- c(meanlog.true=-1.5,meanlog.false=0,
               sdlog.true=1,sdlog.false=1,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
plot.score.dmc (data.model)

# 4) Changing meanlog for both accumulators (true and false) equally has no
#    effect on accuracy (the same is true for the LBA model for mean_v). 
#    In this example, because both are increased, responding is slowed for all 
#    responses compared to (1). 
p.vector  <- c(meanlog.true=-.5,meanlog.false=.5,
               sdlog.true=1,sdlog.false=1,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
plot.score.dmc (data.model)

# 5) A decrease in sdlog.true mainly increases accuracy compared to (1).
# It also speeds up mean error RT, so now errors are fast.
p.vector  <- c(meanlog.true=-0.5,meanlog.false=0.5,
               sdlog.true=0.5,sdlog.false=1,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
plot.score.dmc (data.model)

# 6) In contrast to meanlog, equal changes in sdlog affect both accuracy and 
# average speed. In this example an increase in both decreases accuracy and
# slows mean RT. 
p.vector  <- c(meanlog.true=-1,meanlog.false=0,
               sdlog.true=2,sdlog.false=2,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
plot.score.dmc (data.model)
# N.B. If you decrease both you will see that accuracy increases but mean RT is 
#      not much affected. In general accuracy effects are complicated. 

# You can (with all models) explore a parameter's effects in the same data set
# using a more complicated design. For example, suppose you add a factor, D,
# with levels "easy" and "hard" for stimulus discriminations of different 
# difficulty, with an effect on meanlog.true:
model <- model.dmc(type="lnr",constants=c(st0=0),
  p.map=list(meanlog=c("D","M"),sdlog="M",t0="1",st0="1"),
  match.map=list(M=list(left="LEFT",right="RIGHT")),responses=c("LEFT","RIGHT"),
  factors=list(S=c("left","right"),D=c("easy","hard")))
# [1] "meanlog.easy.true"  "meanlog.hard.true"  "meanlog.easy.false" "meanlog.hard.false"
# [5] "sdlog.true"         "sdlog.false"        "t0"   

p.vector  <- c(meanlog.easy.true=-1.0,meanlog.easy.false=0,
               meanlog.hard.true=-0.5,meanlog.hard.false=0,
               sdlog.true=1,sdlog.false=1,t0=.2)

data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
# Score & plot
plot.score.dmc (data.model)

correct <- data.model$S==tolower(data.model$R)
# Easy more accurate than hard
round(tapply(correct,data.model$D,mean),2)   
# Easy faster than hard
round(tapply(data.model$RT,list(correct,data.model$D),mean),2) 

# Now lets plot RT densities. These are like histograms, but with smooth
# lines, see R help ?density
par(mfrow=c(1,2))
plot.cell.density(data.model[data.model$D=="easy",],
                  C=correct[data.model$D=="easy"],
                  xlim=c(0,5),ymax=2.0,main="easy")
plot.cell.density(data.model[data.model$D!="easy",],
                  C=correct[data.model$D!="easy"],
                  xlim=c(0,5),ymax=2.0,main="hard")

