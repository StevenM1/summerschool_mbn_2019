# Standard DDM model
#   External parameters types: a, v, t0, z = 0.5, d = 0, sz = 0, sv = 0, st0 = 0
#   Internal parameters types: a, v, t0, z = 0.5, d = 0, sz = 0, sv = 0, st0 = 0

# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all appliciatons

# This function transfroms parameters to a form suitbale for the model 
#   being used. Called inside of get.par.mat. 
# "par.df" is a data frame of parameters types , some of which may need to be 
#   transformed, or new columns created, so that the full set of internal 
#   parameter types, specified in "type.par.names", required by the type of 
#   evidence accumulation model being used ("norm" etc.) is present.
transform.dmc <- function(par.df) 
{

  # User supplied tranforms go here
  par.df$sv <- 1/par.df$tau

#   # COMMENT OUT this check for speed after debugging
#   if ( !all(type.par.names %in% names(par.df)) )
#     stop("Trasform has not created parameter(s) required by the model.")
  
  par.df[,c("a","v","t0","z","d","sz","sv","st0")]
}

random.dmc(n,p.df,model)
{
rdiffusion(n,a=p.df$a[1],v=p.df$v[1],t0=p.df$t0[1],
           sz=p.df$sz[1]*p.df$a[1],z=p.df$z[1]*p.df$a[1], # convert to absolute
           d=p.df$d[1],sv=p.df$sv[1],st0=p.df$st0[1])
}


likelihood.dmc <- function(p.vector,data,ok.types=NULL,min.like=1e-10,
                           precision=2.5) 
# Returns vector of likelihoods for each RT in data (in same order). Precision 
# usually in range 2 (fast but less accruate) to 3 (slower but more accurate)
{  
  bad <- function(p) 
    # Stops drd crashing if given bad values, may need to add extra protection
    # for v and upper bounds for a, sv and t0.
  {
    (p$a[1]<0)      | (p$z[1] <1e-6) | (p$z[1] >.999999) | (p$t0[1]<1e-6)  | 
      (p$sz[1]<0) | (p$st0[1]<0)    | (p$sv[1]<0) |
      (p$sz[1]>1.999999*min(c(p$z[1],1-p$z[1]))) 
  }
  
  bound <- rep("lower",dim(attr(data,"model"))[1])
  bound[as.vector(sapply(
    paste("",names(attr(attributes(data)$model,"match.map")$M),
    attr(attributes(data)$model,"match.map")$M,sep="*"),
    function(x){grep(glob2rx(x),row.names(attr(data,"model")))}))] <- "upper"
  names(bound) <- row.names(attr(data,"model"))
  
  likelihood <- numeric(dim(data)[1])
  for ( i in row.names(attr(data,"model")) ) if ( !attr(data,"cell.empty")[i] )
  {
    p.df <- p.df.dmc(p.vector,i,attributes(data)$model,n1order=TRUE) 

    # convert to absolute z. 
    p.df$z <- p.df$z*p.df$a

    # Avoid numerical errors
    if ( bad(p.df) ) likelihood[ attr(data,"cell.index")[[i]] ] <- 0 else 
    {
      likelihood[ attr(data,"cell.index")[[i]] ] <-
        abs(ddiffusion(data$RT[attr(data,"cell.index")[[i]]],
          a=p.df$a[1],
          v=p.df$v[1],
          t0=p.df$t0[1], 
          z=p.df$z[1]*p.df$a[1],   # convert to absolute z.
          d=p.df$d[1],
          sz=p.df$sz[1]*p.df$a[1], # convert to absolute sz. 
          sv=p.df$sv[1],
          st0=p.df$st0[1],
          response=bound[i],
          precision=precision
        ))
    }
 }
 pmax(likelihood,min.like)
}



