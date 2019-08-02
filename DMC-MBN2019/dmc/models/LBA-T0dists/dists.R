require(rtdists) 

####### LBA Gamma non-decision, MR linear angle effects (AS) on shape  ----

{
  rlba.norm.gammaMR <- function(n,A,b,t0,mean_v,sd_v,shape,scale,AS,posdrift=TRUE) 
  {
    shape <- shape*AS  
    if (any(b < A)) 
      stop("b cannot be smaller than A!")
    n_v <- length(mean_v)
    if (posdrift) 
      drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v, 
                              lower = 0), nrow = n_v) else 
                                drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v), nrow = n_v)
    out <- make.r(drifts=drifts,b=b,A=A,n_v=n_v,t0=t0,n = n)
    if (shape>0 & scale>0) 
      cbind.data.frame(RT=out$rt + rgamma(n,shape=shape[1],scale=scale[1]),
                       R=factor(out$response)) else
                         cbind.data.frame(RT=out$rt,R=factor(out$response))
  }
  
  
  n1PDF.norm.gammaMR <- function (t,A,b,mean_v,sd_v,t0,shape,scale,AS, 
                                  posdrift=TRUE,robust = FALSE)
  {
    shape <- shape*AS  
    dt <- pmax(t-t0[1], 0)
    if (shape[1] == 0 | scale[1] == 0) 
      return(n1PDFfixedt0.norm.old(dt=dt,A=A,b=b,mean_v=mean_v,sd_v=sd_v, 
                                   posdrift=posdrift,robust = robust))
    else {
      tmpf <- function(tau,ti0,A,b,mean_v,sd_v,scale,shape) 
        n1PDFfixedt0.norm.old(tau,A,b,mean_v,sd_v)*
        dgamma(ti0-tau,shape=shape,scale=scale)
      outs <- numeric(length(dt))
      for (i in c(1:length(outs))) {
        tmp <- try(integrate(f=tmpf,
                             lower=0,upper=dt[i],ti0=dt[i],
                             A=A,b=b,mean_v=mean_v,sd_v=sd_v,
                             scale=scale[1],shape=shape[1])$value,silent=T)
        if (is.numeric(tmp)) 
          outs[i] <- tmp else
            outs[i] <- 0
      }
      return(outs)
    }
  }
  
  
  # # Check
  # n=1e6
  # shape=1; scale=1
  # AS=0
  # AS=1
  # t0=.2; A=c(.5,.5); b=c(1,1); mean_v=c(1,.5); sd_v=c(1,1)
  # sim <- rlba.norm.gammaMR(n=n,A=A,b=b,mean_v=mean_v,sd_v=sd_v,
  #                        shape=shape,scale=scale,t0=t0,AS=AS)
  # dns <- plot.cell.density(sim,C=1,xlim=c(0,7),save.density=TRUE)
  # # n1PDF check corret
  # d <- n1PDF.norm.gammaMR(dns$correct$x,A=A,b=b,mean_v=mean_v,sd_v=sd_v,t0=t0,
  #                       shape=shape,scale=scale,AS=AS)
  # # red=black?
  # plot(dns$correct$x,dns$correct$y,type="l",ylab="density",xlab="RT")
  # lines(dns$correct$x,d,col="red")
  # # n1PDF check, error
  # d <- n1PDF.norm.gammaMR(dns$error$x,A=A[2:1],b=b[2:1],
  #                       mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0,
  #                       shape=shape,scale=scale,AS=AS)
  # lines(dns$error$x,dns$error$y,lty=2)
  # lines(dns$error$x,d,col="red",lty=2)

}

####################### LBA Lognormal non-decision ----
{
  rlba.norm.lnorm <- function(n,A,b,t0,mean_v,sd_v,meanlog,sdlog,posdrift=TRUE) 
  {
    if (any(b < A)) 
      stop("b cannot be smaller than A!")
    n_v <- length(mean_v)
    if (posdrift) 
      drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v, 
                              lower = 0), nrow = n_v)
    else drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v), 
                          nrow = n_v)
    out <- make.r(drifts=drifts,b=b,A=A,n_v=n_v,t0=t0,n = n)
    if (sdlog>0) 
      cbind.data.frame(RT=out$rt + rlnorm(n,meanlog=meanlog,sdlog=sdlog),
                       R=factor(out$response)) else
                         cbind.data.frame(RT=out$rt,R=factor(out$response))
  }
  
  n1PDF.norm.lnorm <- function (t,A,b,mean_v,sd_v,t0,meanlog,sdlog, 
                                posdrift=TRUE,robust = FALSE)
  {
    dt <- pmax(t-t0[1], 0)
    if (sdlog[1] == 0) 
      return(n1PDFfixedt0.norm.old(dt=dt,A=A,b=b,mean_v=mean_v,sd_v=sd_v, 
                                   posdrift=posdrift,robust = robust))
    else {
      tmpf <- function(tau,ti0,A,b,mean_v,sd_v,meanlog,sdlog) 
        n1PDFfixedt0.norm.old(tau,A,b,mean_v,sd_v)*
        dlnorm(ti0-tau,meanlog=meanlog,sdlog=sdlog)
      outs <- numeric(length(dt))
      for (i in c(1:length(outs))) {
        tmp <- try(integrate(f=tmpf,
                             lower=0,upper=dt[i],ti0=dt[i],
                             A=A,b=b,mean_v=mean_v,sd_v=sd_v,
                             meanlog=meanlog[1],sdlog=sdlog[1])$value,silent=T)
        if (is.numeric(tmp)) 
          outs[i] <- tmp else
            outs[i] <- 0
      }
      return(outs)
    }
  }
  
  # # Check
  # n=1e6
  # sdlog=1; meanlog=0
  # t0=.2; A=c(.5,.5); b=c(1,1); mean_v=c(1,.5); sd_v=c(1,1)
  # sim <- rlba.norm.lnorm(n=n,A=A,b=b,mean_v=mean_v,sd_v=sd_v,
  #                        meanlog=meanlog,sdlog=sdlog,t0=t0)
  # dns <- plot.cell.density(sim,C=1,xlim=c(0,7),save.density=TRUE)
  # # n1PDF check, correct
  # d <- n1PDF.norm.lnorm(dns$correct$x,A=A,b=b,mean_v=mean_v,sd_v=sd_v,t0=t0,
  #                       meanlog=meanlog,sdlog=sdlog)
  # # red=black?
  # plot(dns$correct$x,dns$correct$y,type="l",ylab="density",xlab="RT")
  # lines(dns$correct$x,d,col="red")
  # # n1PDF check, error
  # d <- n1PDF.norm.lnorm(dns$error$x,A=A[2:1],b=b[2:1],
  #                       mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0,
  #                       meanlog=meanlog,sdlog=sdlog)
  # lines(dns$error$x,dns$error$y,lty=2)
  # lines(dns$error$x,d,col="red",lty=2)
}