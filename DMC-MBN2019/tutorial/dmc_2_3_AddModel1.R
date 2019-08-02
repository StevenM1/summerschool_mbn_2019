##################  DMC Lesson 2: Priors, Posteriors and Adding New Models
#
# THIS IS AN ADVANCED LESSON, BEST DONE AFTER WORKING THROUGH TO THE END OF LESSON 4
#
# Note: Before running these scripts, the current working directory must be set 
# to the top-level folder containing the dmc and tutorial subfolders 
#
### Lesson 2.3:  Basics of adding models and modifying transform.dmc

# DMC allows you to add a new class of models and new model parameterizations
# within an existing class of models. Each class of models has its own
# sub-directory in the dmc/models. That directory must contain a file called 
# dists.R, which enables the addition of the essential ingredients for a model
# class, a random function (which provides a way of simulating the model) and a
# function for calculating likelihoods. In some cases no detail needs to be 
# specified because this is all taken care of by a package, usually "rtdists".

# For example dists.R in the DDM directory just loads rtdists. For the LBA it is
# a combination of rtdists and some custom functions (actually modified rtdists
# functions adding some functionality not available in early versions of the 
# package). For LNR dists.R defines everything required.

# The functions made available by dists.R are then used in a model file to 
# define the three standard functions required by DMC. The first "likelihood.dmc"
# has already been discussed. The other two, "random.dmc" and "transform.dmc" are
# used internally, with the user calling the former through simulate.dmc (as
# also described in earlier lessons) and transform.dmc, which will be discussed
# in detail here. Note it is this re-use of these functions under different 
# models that precludes DMC being a package. That would require implementing
# types, which might be done in the future.

# To make a new model file you need to define these three functions. We will use
# the lnr.R function in the LNR directory as a base from which to do this.

rm(list=ls())
source ("dmc/dmc.R")

# Usually we would use a call to load_model ("LNR","lnr.R"). When you add your 
# own model class you just need to call load_model with the name of the 
# directory and model file you created. 

# All that load model does is source the dists.R and lnr.R files in the 
# dmc/models/LNR subdirectory. Lets do this direclty. 
source("dmc/models/LNR/dists.R")
# Rather than running the following line to source the model file 
# source("dmc/models/LNR/lnr.R")
# we will create its components here:

# First consider likelihood.dmc

likelihood.dmc <- function(p.vector,data,min.like=1e-10) 
  # Returns vector of likelihoods for each RT in data (in same order)
{
  
  likelihood <- numeric(dim(data)[1])
  for ( i in row.names(attr(data,"model")) ) if ( !attr(data,"cell.empty")[i] )
  {
    p.df <- p.df.dmc(p.vector,i,attributes(data)$model,n1order=TRUE)
    likelihood[ attr(data,"cell.index")[[i]] ] <-
      n1PDF.lnr(data$RT[attr(data,"cell.index")[[i]]],
                t0=p.df$t0, 
                meanlog=p.df$meanlog,
                sdlog=p.df$sdlog,
                st0=p.df$st0[1]
      )
  }
  pmax(likelihood,min.like)
}

# NB: The final line just ensures that a zero likelihood is never returned to
#     avoid numerical issues when its logarithm is calculated.

# To make your own model class all you need do is substitute for n1PDF.lnr the 
# name of the function you made in dists.R to calculate likelihoods and provide
# it with appropriate arguments from the data frame p.df which will have one
# row for the parameter of each accumulator and column names corresponding to 
# the parameter names used by that function (created by p.df.dmc, an internal
# dmc function) from p.vector given an index, i, of the relevant design cell. You 
# can just leave this and the other components which just use internal DMC 
# methods to loop over the cells of the design in the data object. Note that the 
# "n1PDF.lnr" function gives the likelihood for the winning accumulator with 
# parameters given in the first row of p.df ("n1" is "node 1"), so your new 
# function must do the same.

# Second consider the random function

random.dmc<- function(n,p.df,model)
  # Returns a data frame with columns RT and R 
{
  rlnr(n,meanlog=p.df$meanlog,sdlog=p.df$sdlog,
       t0=p.df$t0,st0=p.df$st0[1])
}

# Again you can substitute here your own random function name from rtdists, 
# providing the appropriate parameter names from p.df, which is created by a
# call to p.df.dmc in the simulate.dmc function.

# Finally the transform.dmc function is:

transform.dmc <- function(par.df) 
  # Transforms parameters to a form suitable for the model being used. 
{
  par.df[,c("meanlog","sdlog","t0","st0")]
  
}

# Transform.dmc is called inside p.df.dmc. In this case it effectively does 
# nothing, except that it filters its input data frame, par.df, to only output
# the parameter recognized by random.dmc and likelihood.dmc. If these parameters
# are not present then when you first use model.dmc (as described in
# earlier lessons) an error will be given to indicate something is wrong. 

# To understand the flexibility given by transform.dmc to reparameterize models
# that all use the same dists.R definitions consider transform.dmc from the 
# lnrPP.R model file:

# "par.df" is a data frame of parameters types, some of which may need to be 
#   transformed, or new columns created, so that the full set of internal 
#   parameter types, specified in "type.par.names", required by the type of 
#   evidence accumulation model being used is present.

transform.dmc <- function(par.df) 
  # Transforms parameters to a form suitable for the model being used. 
{
  par.df$t0 <- exp(par.df$t0)  
  par.df$sdlog <- exp(par.df$sdlog)
  par.df[,c("meanlog","sdlog","t0","st0")]
  
}

# In the first line we see that the t0 parameter sampled on an unbounded log 
# seconds scale is exponentiated to the seconds scale expected by the dists.R
# functions. On the second line the same is done for the sdlog parameter.

# If you look at the top of the lnrPP.R model file you will see a comment 
# describing the transform in terms of "External parameters" (created when you
# define a model with model.dmc) and "Internal parameters", which are
# expected by the dists.R functions. It is a good idea to put in such comments
# when you make your own model files.

# transform.dmc is not limited to rescaling existing parameters, it can also
# create new ones. Consider the lba_B.R model file in the LBA directory. There

transform.dmc <- function(par.df) 
{
  par.df$b <- par.df$B+par.df$A
  
  par.df[,c("A","b","t0","mean_v","sd_v","st0")]
}

# Here the internal parameter "b" is created by adding the external B and A
# parameters. The function of the final line is now evident, it filters the
# output to only contain the internal parameter expected by the dists.R 
# functions. This also enables the check made in models.dmc; if the appropriate
# parameters expected by dists.R are not returned an error occurs. 


### EXAMPLE 1

# transform.dmc can also be used to replace particular rows of the par.df data
# frame. As an example, consider a one-dimensional LBA where the rate for the 
# mismatching accumulator is always one minus the rate of the matching
# accumulator as defined in lba_B1v.R in LBA. Here we use a trick where the
# mismatch rate is set to a constant (-Inf) that cannot occur in sampling the
# match rate, and so it allows the mismatch rates to be uniquely identified
# and replaced. Suppose we had the following model definition.


load_model ("LBA","lba_B1v.R")

model <- model.dmc(
  p.map <- list(A="1",B="1",d="1",mean_v="M",sd_v="S",t0="1",st0="1"),
  constants=c(st0=0,mean_v.false=-Inf,sd_v.s1=0.4),
  match.map=list(M=list(s1=1,s2=2)),
  factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="norm")

# Note that for the LBA one dimensional model there is no need to fix any 
# parameters for scaling purposes, but in order to enforce make 
# v.false = 1-v.true need to set all v.false to a constant -Inf (for more 
# complicated designs this is true for EVERY v.false parameter)

# In transform.dmc all rows in the mean_v column with the value -Inf are 
# replaced with one minus the value in the remaining row. Note that this trick
# would work for an N choice LBA, but only where all mismatching accumulators
# have the same rate.
transform.dmc <- function(par.df) 
{
  par.df$b <- par.df$B+par.df$A
  par.df$mean_v[par.df$mean_v==-Inf] <- 1-par.df$mean_v[par.df$mean_v!=-Inf] 
  
  par.df[,c("A","b","t0","mean_v","sd_v","st0")]
}

# The LBA model directory contains a lesson for a more general form lba_Bdv.R which 
# (lba_Bdv_lesson.R) which sets the mismatch rate at d - match rate, where d is 
# the sum of match and mismatch rates, which determines the overall "drive" (i.e. 
# strength of accumulation). When d is estimated something will need to be set 
# to a constant for scaling purposes. Note that the "1v" model is just a special 
# case of the "dv" model with a constant d=1.

### EXAMPLE 2

# As another example of this sort of approach consider the lba_B2v.R model file
# in LBA. This allows the flexibility to replace selected mean_v values with a 
# separately estimated parameter. To do so we create a new external parameter 
# mean_v1 (any name can be used as long as it is appropriately dealt with by 
# transform.dmc).

load_model ("LBA","lba_B2v.R")
model <- model.dmc(
  p.map=list(A="1",B="1",mean_v=c("S","M"),mean_v1="1",sd_v=c("M"),t0="1",st0="1"),
  constants=c(st0=0,sd_v.true=1,mean_v.s1.false=Inf),
  match.map=list(M=list(s1="r1",s2="r2")),
  factors=list(S=c("s1","s2")),responses=c("r1","r2"),type="norm")

# This model allows us to set the same value for mismatch rates for both 
# levels of the factor S but different values for the match rates. Note that
# we again set the rows in mean_v to be replaced with distinctive numeric values 
# (Inf, in general both -Inf and Inf are easy options, or you could set a prior
# to exclude some range of values and use any value in that range, giving you
# unlimited flexibility). The transform.dmc in lba_B2v.R is

transform.dmc <- function(par.df) 
{
  par.df$b <- par.df$B+par.df$A
  
  # Replace mismatch in mean_v (set to Inf as a constant) with mean_v1
  par.df$mean_v[par.df[,"mean_v"]==Inf] <- par.df[1,"mean_v1"]
  
  par.df[,c("A","b","t0","mean_v","sd_v","st0")]
}

# Note that you could use a similar approach to have separate parameters for 
# B for each accumulator, but now using B="R" in the model.dmc call and 
# then setting e.g., B.r1 to a distinctive constant. 

# These approaches do become somewhat unwieldy when mean_v or sd_v (which can
# be treated like mean_v) or B or A (which can be treated like B) are functions
# of many factors. However, this method does provide a work around for DMC that
# means you do not have to change DMCs internal methods of dealing with 
# the accumulator factors (M and R).

