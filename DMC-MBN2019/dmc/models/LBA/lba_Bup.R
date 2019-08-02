# Template setup for n-choice LBA, B=b-A unbounded parameterization (except st0 so can fix at zero)
#   External parameters types: log(A), log(B), log(t0), mean_v, log(sd_v), st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all applications

transform.dmc <- function(par.df) 
{
  par.df$A <- exp(par.df$A)
  par.df$B <- exp(par.df$B)
  par.df$b <- par.df$B+par.df$A
  par.df$t0 <- exp(par.df$t0)
  par.df$sd_v <- exp(par.df$sd_v)
  par.df[,c("A","b","t0","mean_v","sd_v","st0")]
}


random.dmc <- function(n,p.df,model)
{
  rlba.norm(n,A=p.df$A,b=p.df$b,t0=p.df$t0, 
            mean_v=p.df$mean_v,sd_v=p.df$sd_v,st0=p.df$st0[1],
            posdrift = attr(model,"posdrift"))
}


likelihood.dmc <- function(p.vector,data,ok.types=c("norm"),min.like=1e-10)   
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


