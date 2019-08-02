# n-choice LNR, t0, mu, sigma parameterization, with Go-NoGo
#    External parameters types: meanlog, sdlog, t0, st0 (optional)
#    Internal parameters types: meanlog, sdlog, t0, st0 (optional)

# User edited funcitons for the DMC (Dynamic Models of Choice)
#    Source in all appliciatons

my.integrate <- function(f,lower,...) {
   out <- try(integral(fun=f,xmin=lower,xmax=Inf,method="Kronrod",...) ,silent=TRUE)
   if (class(out)=="try-error") 0 else out
}

# source("rtdists_extras.R")

transform.dmc <- function(par.df) 
# This function transfroms parameters to a form suitbale for the model 
#   being used. Called inside of get.par.mat. 
# "par.df" is a data frame of parameters types , some of which may need to be 
#   transformed, or new columns created, so that the full set of internal 
#   parameter types, specified in "type.par.names", required by the type of 
#   evidence accumulation model being used is present.
{
  par.df[,c("meanlog","sdlog","t0","st0")]
  
}

random.dmc<- function(n,p.df,model)
{rlnr(n,meanlog=p.df$meanlog,sdlog=p.df$sdlog,
      t0=p.df$t0,st0=p.df$st0[1])
}

likelihood.dmc <- function(p.vector,data,min.like=1e-10) 
# Returns vector of likelihoods for each RT in data (in same order)
{
  
  likelihood <- numeric(dim(data)[1])
  for ( i in row.names(attr(data,"model")) ) if ( !attr(data,"cell.empty")[i] )
  {
    p.df <- p.df.dmc(p.vector,i,attributes(data)$model,n1order=TRUE)
    likelihood[ attr(data,"cell.index")[[i]] ] <-
      n1PDF.lnrgng(data$RT[attr(data,"cell.index")[[i]]],
          t0=p.df$t0, 
          meanlog=p.df$meanlog,
          sdlog=p.df$sdlog,
          st0=p.df$st0[1]
      )
 }
 pmax(likelihood,min.like)
}


