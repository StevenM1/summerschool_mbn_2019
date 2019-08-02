# Template setup for n-choice LBA
#   External parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all applications

#source("rtdists_extras.R")

# This function transfroms parameters to a form suitbale for the model 
#   being used. Called inside of get.par.mat. 
# "par.df" is a data frame of parameters types , some of which may need to be 
#   transformed, or new columns created, so that the full set of internal 
#   parameter types, specified in "type.par.names", required by the type of 
#   evidence accumulation model being used ("norm" etc.) is present.
transform.dmc <- function(par.df) 
{
  # User supplied tranforms go here

  par.df[,c("A","b","t0","mean_v","sd_v","st0")]
}


random.dmc <- function(n,p.df,model)
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


