# Distance and Time approach avoidacne model
# Hyperbolic time + linear Distance

# Template setup for n-choice LBA, B=b-A parameterization
#   External parameters types: A, B, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all applications

# source("rtdists_extras.R")


# This function transfroms parameters to a form suitbale for the model 
#   being used.
# "par.df" is a data frame of parameters types , some of which may need to be 
#   transformed, or new columns created, so that the full set of internal 
#   parameter types, specified in "type.par.names", required by the type of 
#   evidence accumulation model being used ("norm" etc.) is present.

# FOR EXAMPLE par.df might look like
#   A B  t0 sd_v st0 Wt .Wd mean_v
# A 1 1 0.2    1   0  1  .1   1112
# B 1 2 0.2    1   0  1  .1   1112
#
# Note that row names are responses, and mean_v is the paraemter to be 
# mapped as a function of the cell coordinates in mean_v, which are assumed to
# index a look table for A.Time (1000s), A.distance (100s), B.time (10s)
# and B.distance (units). Note must use the response names A,B for this to work.

# Pass in lookup values as constants
transform.dmc <- function(par.df,
  #AT=c('1'=1,'2'=2,'3'=4,'4'=8),   BT=c('1'=1,'2'=2,'3'=4,'4'=8),   # Could just code T and D as 
  #AD=c('1'=10,'2'=20,'3'=40,'4'=80), BD=c('1'=10,'2'=20,'3'=40,'4'=80)) # A and B assumed the same.

  AT=c('1'=18,'2'=20,'3'=24,'4'=30,'5'=40,'6'=60,'7'=90),
  BT=c('1'=18,'2'=20,'3'=24,'4'=30,'5'=40,'6'=60,'7'=90),   # Could just code T and D as 
  AD=c('1'=30,'2'=45,'3'=67.5,'4'=90,'5'=112.5,'6'=135,'7'=150), 
  BD=c('1'=30,'2'=45,'3'=67.5,'4'=90,'5'=112.5,'6'=135,'7'=150))
  
{
  
  # Maps 1000s to AT, 100s to AD, 10s to BT and 1s to BD
  lookup <- strsplit(as.character(par.df["r1","mean_v"]),"")[[1]]
  
  # Set mean_v for A accumulator 
  par.df["r1","mean_v"] <- par.df["r1","Wt"]*(AT[lookup[1]])^par.df["r1","tau"]
  
  # Set mean_v for B accumulator 
  par.df["r2","mean_v"] <- par.df["r2","Wt"]*(BT[lookup[3]])^par.df["r2","tau"]
  
  
  # Usual transfrom for _B form of LBA
  par.df$b <- par.df$B+par.df$A
  
  # Filter output to internal parameters
  par.df[,c("A","b","t0","mean_v","sd_v","st0")]
}

random.dmc<- function(n,p.df,model)
{
  rlba.norm(n,A=p.df$A,b=p.df$b,t0=p.df$t0, 
            mean_v=p.df$mean_v,sd_v=p.df$sd_v,st0=p.df$st0[1],
            posdrift = attr(model,"posdrift"))
}



likelihood.dmc <- function(p.vector,data,ok.types=c("norm"),min.like=1e-10) 
# Returns vector of likelihoods for each RT in data (in same order)
{

  likelihood <- numeric(dim(data)[1])
  for ( i in row.names(attr(data,"model")) ) if ( !attr(data,"cell.empty")[i] )
  {
    p.df <- p.df.dmc(p.vector,i,attributes(data)$model,n1order=TRUE)
    likelihood[ attr(data,"cell.index")[[i]] ] <-
      switch(attr(attributes(data)$model,"type"),
        norm=n1PDF(data$RT[attr(data,"cell.index")[[i]]],
          A=as.list(p.df$A),
          b=as.list(p.df$b),
          t0=p.df$t0[1], 
          mean_v=p.df$mean_v,
          sd_v=p.df$sd_v,
          st0=p.df$st0[1],
          distribution=attr(attr(data,"model"),"type"),
          args.list = list(posdrift=attr(attr(data,"model"),"posdrift")),
          silent=TRUE)
      )
 }
 pmax(likelihood,min.like)
}


