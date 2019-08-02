# Template setup for n-choice Wald, B=b-A parameterization with probit go failure
# Note that B is internal.
# External parameters types: v, B, A, t0, probit gf
# Internal parameters types: v, B, A, t0, gf

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
  
  par.df[,c("A","B","t0","v","gf")]
}

random.dmc<- function(n,p.df,model)
{
  rWaldRace(n,v=p.df$v,B=p.df$B,A=p.df$A,t0=p.df$t0,gf=p.df$gf[1])
}
  



likelihood.dmc <- function(p.vector,data,ok.types=c("wald"),min.like=1e-10) 
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
        wald=n1Wald(data$RT[attr(data,"cell.index")[[i]]],
          A=p.df$A,
          B=p.df$B,
          t0=p.df$t0[1], 
          v=p.df$v,
          gf=p.df$gf[1])
      )
 }
 pmax(likelihood,min.like)
}


