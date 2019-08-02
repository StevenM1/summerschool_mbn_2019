# Template setup for n-choice LBA, B=b-A parameterization, with Go-NoGo
#   External parameters types: A, B, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# User edited funcitons for the DMC (Dynamic Models of Choice)
#    Source in all appliciatons

my.integrate <- function(f,lower,...) {
   out <- try(integral(fun=f,xmin=lower,xmax=Inf,method="Kronrod",...) ,silent=TRUE)
   if (class(out)=="try-error") 0 else out
}

# This function transfroms parameters to a form suitbale for the model 
#   being used. Called inside of get.par.mat. 
# "par.df" is a data frame of parameters types , some of which may need to be 
#   transformed, or new columns created, so that the full set of internal 
#   parameter types, specified in "type.par.names", required by the type of 
#   evidence accumulation model being used ("norm" etc.) is present.
transform.dmc <- function(par.df) 
{
  # User supplied tranforms go here
  par.df$b <- par.df$B+par.df$A
  
#   # COMMENT OUT this check for speed after debugging
#   if ( !all(type.par.names %in% names(par.df)) )
#     stop("Trasform has not created parameter(s) required by the model.")
  
  par.df[,c("A","b","t0","mean_v","sd_v","st0")]
}

random.dmc<- function(n,pdf,model)
{
  rlba.norm(n,A=p.df$A,b=p.df$b,t0=p.df$t0, 
            mean_v=p.df$mean_v,sd_v=p.df$sd_v,st0=p.df$st0[1],
            posdrift = attr(model,"posdrift")) 
}

# Returns vector of likelihoods for each RT in data (in same order)
# !!! TO DO !!!
# !!! types other than norm
# !!! posdrift has dropped out !!!

likelihood.dmc <- function(p.vector,data,ok.types=c("normgng"),min.like=1e-10) 
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
        normgng=n1PDF.normgng(data$RT[attr(data,"cell.index")[[i]]],
          A=p.df$A,
          b=p.df$b,
          t0=p.df$t0, 
          mean_v=p.df$mean_v,
          sd_v=p.df$sd_v,
          st0=p.df$st0[1])
      )
  }
  pmax(likelihood,min.like)
}


