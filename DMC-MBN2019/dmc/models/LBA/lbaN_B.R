# Template setup for n-choice LBA, B=b-A parameterization
#   External parameters types: A, B, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# Allows for a variable number of accumulators by passing parameter for number
# of accumulators in p.df$N[1]


# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all applications

transform.dmc <- function(par.df) 
{
  par.df$b <- par.df$B+par.df$A
  
  par.df[,c("A","b","t0","mean_v","sd_v","st0","N")]
}

random.dmc<- function(n,p.df,model)
{
  rlba.norm(n,
    A=p.df$A[1:p.df$N[1]],
    b=p.df$b[1:p.df$N[1]],
    mean_v=p.df$mean_v[1:p.df$N[1]],
    sd_v=p.df$sd_v[1:p.df$N[1]],
    t0=p.df$t0[1], 
    st0=p.df$st0[1],
    posdrift = attr(model,"posdrift"))
}



likelihood.dmc <- function(p.vector,data,min.like=1e-10) 
# Returns vector of likelihoods for each RT in data (in same order)
{

  likelihood <- numeric(dim(data)[1])
  for ( i in row.names(attr(data,"model")) ) if ( !attr(data,"cell.empty")[i] )
  {
    p.df <- p.df.dmc(p.vector,i,attributes(data)$model,n1order=TRUE)
    likelihood[ attr(data,"cell.index")[[i]] ] <-
      n1PDF(data$RT[attr(data,"cell.index")[[i]]],
          A=as.list(p.df$A[1:p.df$N[1]]),
          b=as.list(p.df$b[1:p.df$N[1]]),
          mean_v=p.df$mean_v[1:p.df$N[1]],
          sd_v=p.df$sd_v[1:p.df$N[1]],
          t0=p.df$t0[1], 
          st0=p.df$st0[1],
          distribution="norm",
          args.list = list(posdrift=attr(attr(data,"model"),"posdrift")),
          silent=TRUE)
      
  }
  pmax(likelihood,min.like)
}


