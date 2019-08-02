#####################  LBA with go failure
require(rtdists) 

#######  Standard model ----

make.r <- function (drifts, b, A, n_v, t0, st0 = 0, n, return.ttf=FALSE) 
{
    drifts[drifts < 0] <- 0
    starts <- matrix(runif(min = 0, max = A, n = n * n_v), nrow = n_v)
    ttf <- t0 + (b - starts)/drifts
    if (return.ttf) return(ttf)
    resp <- apply(ttf, 2, which.min)
    rt <- ttf[cbind(resp,1:n)]
    if (st0[1]>0) rt <- rt + runif(min = 0, max = st0[1], n = n)
    bad <- !is.finite(rt)
    if (any(bad)) {
        warning(paste(sum(bad), "infinite RTs removed and less than", 
            n, "rts returned"))
        resp <- resp[!bad]
        rt <- rt[!bad]
    }
    data.frame(rt = rt, response = resp)
}


rlba.norm <- function (n,A,b,t0,mean_v,sd_v,st0=0,posdrift=TRUE,return.ttf=FALSE) 
{
    if (any(b < A)) 
        stop("b cannot be smaller than A!")
    n_v  <- ifelse(is.null(dim(mean_v)), length(mean_v), dim(mean_v)[1])
    if (posdrift) 
        drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v, 
            lower = 0), nrow = n_v)
    else drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v), 
        nrow = n_v)
    make.r(drifts = drifts, b = b, A = A, n_v = n_v, t0 = t0, st0 = st0, n = n,
      return.ttf=return.ttf)
}



#####################  LBA with go failure ----

rlba.normGF <- function(n,A,b,t0,mean_v,sd_v,gf=0,st0=0,posdrift=TRUE,return.ttf=FALSE) 
  # random function for LBA race with go failure.
{
  out <- rlba.norm (n=n,A=A,b=b,t0=t0,mean_v=mean_v,sd_v=sd_v,st0=st0,
                    posdrift=posdrift,return.ttf=return.ttf)     
  if (return.ttf) return(out)
  names(out) <- c("RT","R")
  if (gf > 0) {
    is.gf <- as.logical(rbinom(dim(out)[1],1,gf))
    out$RT[is.gf] <- NA
    out$R[is.gf] <- 1
  }
  
  out
}


n1PDFGF <- function(dt,A,b,mean_v,sd_v,t0,gf=0,st0=0,
                    distribution="norm",args.list=list(posdrift=TRUE),silent=TRUE)
  # Generates defective PDF for responses on node=1, with go failure 
  # dt (decison time) is a vector of times
{
  out <- numeric(length(dt))
  is.go <- !is.na(dt)
  
  out[is.go] <- (1-gf[1])*n1PDF(dt[is.go],A=A,b=b,t0=t0,mean_v=mean_v,sd_v=sd_v,st0=st0,
                                distribution=distribution,args.list=args.list,silent=silent)
  
  out[!is.go] <- gf[1]
  
  out
}


# # Check go failure 
# n=1e5
# v=c(2,1); sd_v = c(1,1); B=c(1,1); A=c(2,2);t0=1; gf=.2
# b=A+B
# sim <- rlba.normGF(n=n,A=A,b=b,mean_v=v,sd_v=sd_v,t0=t0,gf=gf)
# par(mfrow=c(1,3))
# dns <- plot.cell.density(sim,C=1,xlim=c(0,4),save.density=TRUE)
# pc <- integrate(n1PDFGF,lower=0,upper=Inf,A=A,b=b,mean_v=v,sd_v=sd_v,t0=t0,gf=gf)$value/(1-gf)
# print(pc)
# dt=dns$correct$x
# d <- n1PDFGF(dt,A=A,b=b,mean_v=v,sd_v=sd_v,gf=gf,t0=t0)
# # red=black?
# plot(dns$correct$x,dns$correct$y,type="l")
# lines(dns$correct$x,d,col="red")
# dt=dns$error$x
# d <- n1PDFGF(dt,A=A[2:1],b=b[2:1],mean_v=v[2:1],sd_v=sd_v[2:1],gf=gf,t0=t0)
# plot(dns$error$x,dns$error$y,lty=2,type="l")
# lines(dns$error$x,d,col="red",lty=2)
# 
# # Check (effectively) single accumulator case 
# n=1e5
# v=c(2,-100); sd_v = c(1,1); B=c(1,1); A=c(2,2);t0=1; gf=.2
# b=A+B
# sim <- rlba.normGF(n=n,A=A,b=b,mean_v=v,sd_v=sd_v,t0=t0,gf=gf)
# par(mfrow=c(1,3))
# dns <- plot.cell.density(sim,C=1,xlim=c(0,7),save.density=TRUE)
# pc <- integrate(n1PDFGF,lower=0,upper=Inf,A=A,b=b,mean_v=v,sd_v=sd_v,t0=t0,gf=gf)$value/(1-gf)
# print(pc)
# dt=dns$correct$x
# d <- n1PDFGF(dt,A=A,b=b,mean_v=v,sd_v=sd_v,gf=gf,t0=t0)
# # red=black?
# plot(dns$correct$x,dns$correct$y,type="l")
# lines(dns$correct$x,d,col="red")
# 
# v=c(-100,2)
# sim <- rlba.normGF(n=n,A=A,b=b,mean_v=v,sd_v=sd_v,t0=t0,gf=gf)
# par(mfrow=c(1,3))
# dns <- plot.cell.density(sim,C=1,xlim=c(0,7),save.density=TRUE)
# pe <- integrate(n1PDFGF,lower=0,upper=Inf,A=A,b=b,mean_v=v,sd_v=sd_v,t0=t0,gf=gf)$value/(1-gf)
# print(pe)
# dt=dns$error$x
# d <- n1PDFGF(dt,A=A[2:1],b=b[2:1],mean_v=v[2:1],sd_v=sd_v[2:1],gf=gf,t0=t0)
# plot(dns$error$x,dns$error$y,lty=2,type="l")
# lines(dns$error$x,d,col="red",lty=2)
