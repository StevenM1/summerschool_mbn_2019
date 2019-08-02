# Template setup for n-choice LBA, B=b-A parameterization, with Go-NoGo
#   External parameters types: A, B, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# User edited funcitons for the DMC (Dynamic Models of Choice)
#    Source in all appliciatons

my.integrate <- function(...,big=10)
# Avoids but in integrate upper=Inf that uses only 1  subdivision
# Use of  big=10 is arbitary ...
{
  out <- try(integrate(...,upper=Inf),silent=TRUE)
  if (class(out)=="try-error") 0 else 
  {
    if (out$subdivisions==1) 
    {
      out <- try(integrate(...,upper=big),silent=TRUE)
      if (class(out)=="try-error") 0 else
      {
         if (out$subdivisions==1) 0 else out$value   
      }
    } else out$value
  }
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

random.dmc<- function(n,p.df,model)
{
  rlba.norm(n,A=p.df$A,b=p.df$b,t0=p.df$t0, 
          mean_v=p.df$mean_v,sd_v=p.df$sd_v,st0=p.df$st0[1],
          posdrift = attr(model,"posdrift"))
}

likelihood.dmc <- function(p.vector,data,ok.types=c("normgng"),min.like=1e-10) 
{
  
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


