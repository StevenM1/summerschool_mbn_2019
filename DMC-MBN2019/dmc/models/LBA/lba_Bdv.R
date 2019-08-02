# n-choice LBA, B=b-A, mean_v(mismatch)=d-mean_v(match) parameterization
#   External parameters types: A, B, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all applications

transform.dmc <- function(par.df) 
{
  par.df$b <- par.df$B+par.df$A
  par.df$mean_v[par.df$mean_v==-Inf] <- 
    par.df$d[par.df$mean_v==-Inf]-par.df$mean_v[par.df$mean_v!=-Inf] 
  
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


