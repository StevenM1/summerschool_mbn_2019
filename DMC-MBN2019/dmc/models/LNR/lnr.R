# n-choice LNR, t0, mu, sigma parameterization
#    External parameters types: meanlog, sdlog, t0, st0 (optional)
#    Internal parameters types: meanlog, sdlog, t0, st0 (optional)

# User edited funcitons for the DMC (Dynamic Models of Choice)
#    Source in all appliciatons

# source("rtdists_extras.R")

transform.dmc <- function(par.df) 
# Transfroms parameters to a form suitbale for the model being used. 
{
  par.df[,c("meanlog","sdlog","t0","st0")]

}

random.dmc<- function(n,p.df,model)
# Retruns a data frame with columns RT and R 
{
  rlnr(n,meanlog=p.df$meanlog,sdlog=p.df$sdlog,
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
      n1PDF.lnr(data$RT[attr(data,"cell.index")[[i]]],
          t0=p.df$t0, 
          meanlog=p.df$meanlog,
          sdlog=p.df$sdlog,
          st0=p.df$st0[1]
      )
 }
 pmax(likelihood,min.like)
}


