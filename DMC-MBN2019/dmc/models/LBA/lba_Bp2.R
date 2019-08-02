# Template setup for n-choice LBA, B=b-A parameterization, with probability bp
# of incrementing using B1+B2 instead of B1 on any trial
#   External parameters types: A, B1, B2, pb2, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all applications

# source("rtdists_extras.R")

transform.dmc <- function(par.df) 
{
  # User supplied tranforms go here
  par.df$b1 <- par.df$B1+par.df$A
  par.df$b2 <- par.df$B2+par.df$b1
  
  par.df[,c("A","b1","b2","pb2","t0","mean_v","sd_v","st0")]
}

random.dmc <- function(n,p.df,model)
{
  n2 <- rbinom(1,n,p.df$pb2[1])
  out <- matrix(nrow=n,ncol=2)
  dimnames(out) <- list(NULL,c("rt","response"))
  out <- data.frame(out)
  if (n2<n) out[1:(n-n2),] <- rlba.norm(n-n2,A=p.df$A,b=p.df$b1,t0=p.df$t0[1],
    mean_v=p.df$mean_v,sd_v=p.df$sd_v,st0=p.df$st0[1],
    posdrift=attr(model,"posdrift"))
  if (n2>0) out[(n-n2+1):n,] <- rlba.norm(n2,A=p.df$A,b=p.df$b2,t0=p.df$t0[1],
    mean_v=p.df$mean_v,sd_v=p.df$sd_v,st0=p.df$st0[1],
    posdrift = attr(model,"posdrift"))
  out
}



likelihood.dmc <- function(p.vector,data,min.like=1e-10) 
# Returns vector of likelihoods for each RT in data (in same order)
{

  if (attr(attr(data,"model"),"type")=="normBp2") dist <- "norm" else dist <- "lnorm"
  likelihood <- numeric(dim(data)[1])
  for ( i in row.names(attr(data,"model")) ) if ( !attr(data,"cell.empty")[i] )
  {
    p.df <- p.df.dmc(p.vector,i,attributes(data)$model,n1order=TRUE)
    likelihood[ attr(data,"cell.index")[[i]] ] <-
      (1-p.df$pb2[1])*n1PDF(data$RT[attr(data,"cell.index")[[i]]],
            A=as.list(p.df$A),b=as.list(p.df$b1),
            mean_v=p.df$mean_v,sd_v=p.df$sd_v,
            t0=p.df$t0[1],st0=p.df$st0[1], 
            distribution=dist,
            args.list = list(posdrift=attr(attr(data,"model"),"posdrift")),
            silent=TRUE) + 
      p.df$pb2[1]*n1PDF(data$RT[attr(data,"cell.index")[[i]]],
            A=as.list(p.df$A),b=as.list(p.df$b2),
            mean_v=p.df$mean_v,sd_v=p.df$sd_v,
            t0=p.df$t0[1],st0=p.df$st0[1], 
            distribution=dist,
            args.list = list(posdrift=attr(attr(data,"model"),"posdrift")),
            silent=TRUE) 
        }

 pmax(likelihood,min.like)
}


