####################### LBA go-nogo ----
require(rtdists) 


# The following functions are from the rtdists package but with some tweaks 
# to vectorize t0.

## Simulate LBA trials

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

### LBA likelihood

dlba.norm.core <- function (t,A,b,mean_v,sd_v,posdrift=TRUE,robust=FALSE,nn) 
# like dlba_norm_core but t0 dealt with outside (removed from t), and
# default added to nn
{
  
  pnormP <- function (x, mean = 0, sd = 1, lower.tail = TRUE) 
    ifelse(abs(x) < 7, pnorm(x, mean = mean, sd = sd, lower.tail = lower.tail), 
      ifelse(x < 0, 0, 1))
  
  dnormP <- function (x, mean = 0, sd = 1) 
    ifelse(abs(x) < 7, dnorm(x, mean = mean, sd = sd), 0)

    if (robust) {
        pnorm1 <- pnormP
        dnorm1 <- dnormP
    }
    else {
        pnorm1 <- pnorm
        dnorm1 <- dnorm
    }
    if (posdrift) 
        denom <- pmax(pnorm1(mean_v/sd_v), 1e-10)
    else denom <- rep(1, nn)
    if (any(A < 1e-10)) {
        A_small <- A < 1e-10
        out_A <- pmax(0, ((b[A_small]/t[A_small]^2) * 
                            dnorm1(b[A_small]/t[A_small], 
            mean_v[A_small], sd = sd_v[A_small]))/denom[A_small])
        zs <- t[!A_small] * sd_v[!A_small]
        zu <- t[!A_small] * mean_v[!A_small]
        chiminuszu <- b[!A_small] - zu
        chizu <- chiminuszu/zs
        chizumax <- (chiminuszu - A[!A_small])/zs
        out_o <- pmax(0, (mean_v[!A_small] * (pnorm1(chizu) - 
            pnorm1(chizumax)) + sd_v[!A_small] * (dnorm1(chizumax) - 
            dnorm1(chizu)))/(A[!A_small] * denom[!A_small]))
        out <- numeric(nn)
        out[!A_small] <- out_o
        out[A_small] <- out_A
        return(out)
    }
    else {
        zs <- t * sd_v
        zu <- t * mean_v
        chiminuszu <- b - zu
        chizu <- chiminuszu/zs
        chizumax <- (chiminuszu - A)/zs
        return(pmax(0, (mean_v * (pnorm1(chizu) - pnorm1(chizumax)) + 
            sd_v * (dnorm1(chizumax) - dnorm1(chizu)))/(A * denom)))
    }
}


plba.norm.core <- function (t,A,b,mean_v,sd_v,posdrift=TRUE,robust=FALSE,nn) 
# like plba_norm_core but t0 dealt with outside (removed from t), and
# default added to nn
{
  
    pnormP <- function (x, mean = 0, sd = 1, lower.tail = TRUE) 
      ifelse(abs(x) < 7, pnorm(x, mean = mean, sd = sd, lower.tail=lower.tail), 
        ifelse(x < 0, 0, 1))
  
    dnormP <- function (x, mean = 0, sd = 1) 
      ifelse(abs(x) < 7, dnorm(x, mean = mean, sd = sd), 0)

    if (robust) {
        pnorm1 <- pnormP
        dnorm1 <- dnormP
    }
    else {
        pnorm1 <- pnorm
        dnorm1 <- dnorm
    }
    if (posdrift) 
        denom <- pmax(pnorm1(mean_v/sd_v), 1e-10)
    else denom <- 1
    if (any(A < 1e-10)) {
        A_small <- A < 1e-10
        out_A <- pmin(1, pmax(0, (pnorm1(b[A_small]/t[A_small], 
            mean = mean_v[A_small], sd = sd_v[A_small], 
            lower.tail = FALSE))/denom[A_small]))
        zs <- t[!A_small] * sd_v[!A_small]
        zu <- t[!A_small] * mean_v[!A_small]
        chiminuszu <- b[!A_small] - zu
        xx <- chiminuszu - A[!A_small]
        chizu <- chiminuszu/zs
        chizumax <- xx/zs
        tmp1 <- zs * (dnorm1(chizumax) - dnorm1(chizu))
        tmp2 <- xx * pnorm1(chizumax) - chiminuszu * pnorm1(chizu)
        out_o <- pmin(pmax(0, (1 + (tmp1 + tmp2)/A[!A_small])/denom[!A_small]), 
            1)
        out <- numeric(length(mean_v))
        out[!A_small] <- out_o
        out[A_small] <- out_A
        return(out)
    }
    else {
        zs <- t * sd_v
        zu <- t * mean_v
        chiminuszu <- b - zu
        xx <- chiminuszu - A
        chizu <- chiminuszu/zs
        chizumax <- xx/zs
        tmp1 <- zs * (dnorm1(chizumax) - dnorm1(chizu))
        tmp2 <- xx * pnorm1(chizumax) - chiminuszu * pnorm1(chizu)
        return(pmin(pmax(0, (1 + (tmp1 + tmp2)/A)/denom), 1))
    }
}


dlba.norm <- function (t, A, b, mean_v, sd_v, posdrift = TRUE, robust = FALSE) 
{
    nn <- length(t)
    tpos <- (t>0) & (b >= A)
    out <- numeric(length(t))
    out[tpos] <- dlba.norm.core(t = t[tpos], A = A[tpos], b = b[tpos], 
      mean_v = mean_v[tpos], sd_v = sd_v[tpos], 
      posdrift = posdrift, robust = robust, nn = nn)
    out
}


plba.norm <- function (t, A, b, mean_v, sd_v, posdrift = TRUE, robust = FALSE) 
{
    nn <- length(t)
    tpos <- (t>0) & (b >= A)
    out <- numeric(length(t))
    out[tpos] <- plba.norm.core(t = t[tpos], A = A[tpos], b = b[tpos], 
      mean_v = mean_v[tpos], sd_v = sd_v[tpos], 
      posdrift = posdrift, robust = robust, nn = nn)
    out
}


n1PDFfixedt0.norm=function(dt,A,b,mean_v,sd_v, 
                           posdrift=TRUE,robust = FALSE) 
# Generates defective PDF for responses on node= 1
# dt (decison time) is a matrix with length(mean_v) rows, one row for
# each accumulator to allow for different start times
{
  
  n_acc <- length(mean_v)
  n <- dim(dt)[2]
  A <- matrix(rep(A,n),nrow=n_acc)
  b <- matrix(rep(b,n),nrow=n_acc)
  mean_v <- matrix(rep(mean_v,n),nrow=n_acc)
  sd_v <- matrix(rep(sd_v,n),nrow=n_acc)
    
  dt[1,] <- dlba.norm(dt[1,],A=A[1,],b=b[1,],mean_v=mean_v[1,],sd_v=sd_v[1,],
                      posdrift=posdrift,robust=robust)
  if (n_acc>1) for (i in 2:n_acc)
    dt[1,] <- dt[1,]*(1-plba.norm(dt[i,],
      A=A[i,],b=b[i,],mean_v=mean_v[i,],sd_v=sd_v[i,],
      posdrift=posdrift,robust=robust))
  dt[1,]
}
 

n1PDFfixedt0.normgng=function(dt,A,b,mean_v,sd_v,st0=0,posdrift=TRUE) 
    # Same as n1PDFfixedt0 except dt=NA done by integration
{
    
    stopfn <- function(t,A,b,mean_v,sd_v,st0=0,posdrift=TRUE) 
    {
      n1PDFfixedt0.norm(
        matrix(rep(t,each=length(mean_v)),nrow=length(mean_v)),
        A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift
      )
    }
    
    is.stop <- is.na(dt[1,])
    if (any(!is.stop)) dt[1,!is.stop] <- n1PDFfixedt0.norm(dt[,!is.stop],
                  A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift)
    # tmp <- try(integrate(f=stopfn,lower=0,upper=Inf,
    #     A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift)$value,silent=TRUE)
    # if (!is.numeric(tmp)) tmp <- 0
    if (any(is.stop)) {
      tmp <- my.integrate(f=stopfn,lower=0,
                        A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift)
      dt[1,is.na(dt[1,])] <- tmp 
    }
    dt[1,]
}

# dt=data$RT[attr(data,"cell.index")[[i]]]
# A=p.df$A;b=p.df$b;t0=p.df$t0;mean_v=p.df$mean_v;sd_v=p.df$sd_v;st0=p.df$st0[1]
# posdrift=TRUE
  
n1PDF.normgng <- function(dt,A,b,t0,mean_v,sd_v,st0=0,posdrift=TRUE)
    # Same as n1PDF except NAs dealt with by integration
{
    
    if ( st0 < 1e-3 ) # smaller values can cause numerical integration problems
      return(n1PDFfixedt0.normgng(
        A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift,
        dt=matrix(pmax(rep(dt,each=length(mean_v))-t0,0),nrow=length(mean_v))
      )) else
      {    
        
        integrate.f <- function(A,b,t0,mean_v,sd_v,st0,posdrift) 
          n1PDFfixedt0.normgng(
            A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift,
            dt=matrix(pmax(rep(dt,each=length(mean_v))-t0,0),
                      nrow=length(mean_v)))/st0
        
        outs <- numeric(length(dt))
        for (i in 1:length(outs)) {
          tmp <- try(integrate(f=integrate.f,
                               lower=dt[i]-st0[1],upper=dt[i],A=A,b=b,t0=t0,mean_v=mean_v,sd_v=sd_v,
                               st0=st0[1],posdrift=posdrift)$value,silent=TRUE)
          if (is.numeric(tmp)) 
            outs[i] <- tmp else
              outs[i] <- 0
        }
        return(outs)
      }
}
  
  
  # # Check
  # 
  # # 2 Choice case
  # n=1e5
  # A=c(1,1);b=c(2,2);t0=.2;mean_v=c(1,0);sd_v=c(1,1);st0=0;posdrift=TRUE
  # sim <- rlba.normgng(n,A,b,t0,mean_v,sd_v,st0,posdrift)
  # names(sim) <- c("RT","R")
  # dns <- plot.cell.density(sim,xlim=c(0,7),save.density=TRUE)
  # d <- n1PDF.normgng(dt=dns$'2'$x,A=A[2:1],b=b[2:1],mean_v=mean_v[2:1],
  #   sd_v=sd_v[2:1],t0=t0,st0=st0,posdrift=posdrift)
  # # d <- n1PDF.normgng(dt=dns$'2'$x,A=A,b=b,mean_v=mean_v,
  # #   sd_v=sd_v,t0=t0,st0=st0,posdrift=posdrift)
  # 
  # lines(dns$'2'$x,d,col="red")
  # 
  # # p(Stop check)
  # mean(is.na(sim$RT))
  # n1PDFfixedt0.normgng(dt=matrix(rep(NA,2),ncol=1),A=A,b=b,
  #                      mean_v=mean_v,sd_v=sd_v)
  # 
  # # 3 Choice case
  # n=1e5
  # A=c(1,1,1);b=c(2,2,2);t0=.2;mean_v=c(1,.5,0);sd_v=c(1,1,1);st0=0;posdrift=TRUE
  # sim <- rlba.normgng(n,A,b,t0,mean_v,sd_v,st0,posdrift)
  # names(sim) <- c("RT","R")
  # dns <- plot.cell.density(sim,xlim=c(0,3),save.density=TRUE)
  # d2 <- n1PDF.normgng(dt=dns$'2'$x,A=A[c(2,3,1)],b=b[c(2,3,1)],mean_v=mean_v[c(2,3,1)],
  #   sd_v=sd_v[c(2,3,1)],t0=t0,st0=st0,posdrift=posdrift)
  # lines(dns$'2'$x,d2,col="red")
  # d3 <- n1PDF.normgng(dt=dns$'3'$x,A=A[3:1],b=b[3:1],mean_v=mean_v[3:1],
  #   sd_v=sd_v[3:1],t0=t0,st0=st0,posdrift=posdrift)
  # lines(dns$'3'$x,d3,col="red",lty=2)
  # 
  # 
  # # p(Stop check)
  # mean(is.na(sim$RT))
  # n1PDFfixedt0.normgng(dt=matrix(rep(NA,3),ncol=1),A=A,b=b,
  #                      mean_v=mean_v,sd_v=sd_v)
  
