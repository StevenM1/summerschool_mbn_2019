# n-choice Wald, with probit go failure and positve normal rate.
# External parameters types: a, v, sv, t0, probit gf
# Internal parameters types: a, v, sv, t0, gf

# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all applications

# source("rtdists_extras.R")


# This function transfroms parameters to a form suitbale for the model 
#   being used. Called inside of get.par.mat. 
# "par.df" is a data frame of parameters types , some of which may need to be 
#   transformed, or new columns created, so that the full set of internal 
#   parameter types, specified in "type.par.names", required by the type of 
#   evidence accumulation model being used ("norm" etc.) is present.
transform.dmc <- function(par.df) 
{
  # User supplied tranforms go here
  par.df$gf <- pnorm(par.df$gf)  
  
#   # COMMENT OUT this check for speed after debugging
#   if ( !all(type.par.names %in% names(par.df)) )
#     stop("Trasform has not created parameter(s) required by the model.")
  
  par.df[,c("a","v","sv","t0","gf")]
}


random.dmc<- function(n,p.df,model)
{
  rPNRWaldRace(n,v=p.df$v,sv=p.df$sv,a=p.df$a,t0=p.df$t0,gf=p.df$gf[1])
}


likelihood.dmc <- function(p.vector,data,ok.types=c("pnrwald"),min.like=1e-10) 
# Returns vector of likelihoods for each RT in data (in same order)
{

#   # COMMENT OUT this check for speed after debugging
#   if ( !all(attr(model,"type") %in% ok.types) )
#     stop("Distribution funciton type not supported by likelihood.dmc")
  
  likelihood <- numeric(dim(data)[1])
  for ( i in row.names(attr(data,"model")) ) if ( !attr(data,"cell.empty")[i] )
  {
    p.df <- p.df.dmc(p.vector,i,attributes(data)$model,n1order=TRUE)
    likelihood[ attr(data,"cell.index")[[i]] ] <-
      switch(attr(attributes(data)$model,"type"),
        pnrwald=n1PNRWald(data$RT[attr(data,"cell.index")[[i]]],
          a=p.df$a,
          t0=p.df$t0[1], 
          v=p.df$v,
          sv=p.df$sv,
          gf=p.df$gf[1])
      )
 }
 pmax(likelihood,min.like)
}


