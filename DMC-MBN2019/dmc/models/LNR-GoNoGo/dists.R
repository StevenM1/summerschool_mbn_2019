####################### LNR go-nogo 


####################### LNR n-choice ----

  rlnr <- function (n, meanlog, sdlog, t0, st0 = 0) 
    # Race among n_acc accumulators, mealnlog and sdlog can be n_acc length
    # vectors or n_acc x n matrices. t0 can be
    # a) a scalar, b) a vector of length number of accumulators or
    # c) a matrix with 1 row per accumulator, when start times differ on each trial
    # st0, range of non-decison time variability, must be a scalar, as the same
    # variability is assumed in a common encoding/production stage
    
  {
    n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])
    dt <- matrix(rlnorm(n = n*n_acc, meanlog = meanlog, sdlog = sdlog), 
                 nrow = n_acc) + t0
    winner <- apply(dt,2,which.min)    
    if (st0[1]==0) data.frame(RT=dt[cbind(winner,1:n)],R=winner) else
      data.frame(RT=dt[cbind(winner,1:n)]+runif(n,0,st0[1]),R=winner)
  }
  
  n1PDFfixedt0.lnr=function(dt,meanlog,sdlog) 
    # Generates defective PDF for responses first among n_acc accumulator at 
    # dt (decison time), a matrix with one row for each accumulator (allowing for
    # different start times per accumulator) 
  {
    
    n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])
    if (!is.matrix(meanlog)) meanlog <- matrix(rep(meanlog,dim(dt)[2]),nrow=n_acc)
    if (!is.matrix(sdlog))     sdlog <- matrix(rep(  sdlog,dim(dt)[2]),nrow=n_acc)
    # winner
    dt[1,] <- dlnorm(dt[1,],meanlog[1,],sdlog[1,])
    # loosers
    if (dim(meanlog)[1]>1) for (i in 2:dim(meanlog)[1])
      dt[1,] <- dt[1,]*plnorm(dt[i,],meanlog[i,],sdlog[i,],lower.tail=FALSE)
    dt[1,]
  }

  
####################### LNR go-nogo ----

  
  rlnrgng <- function (n, meanlog, sdlog, t0, st0 = 0) 
    # Race among length(meanlog) accumulators, first of which is a no-go 
    # accumulator. For trials with winning first accumulator RT set to NA.   
  {
    out <- rlnr(n,meanlog,sdlog,t0,st0)
    out[out$R==1,"RT"] <- NA
    out$R <- factor(as.numeric(out$R))
    data.frame(out)
  }
  
  
  n1PDFfixedt0.lnrgng=function(dt,meanlog,sdlog) 
    # Same as n1PDFfixedt0.lnr except dt=NA done by integration
  {
    
    stopfn <- function(t,meanlog,sdlog) 
    {
      n1PDFfixedt0.lnr(
        matrix(rep(t,each=length(meanlog)),nrow=length(meanlog)),
        meanlog,sdlog
      )
    }
    
    n.trials <- dim(dt)[2]
    out <- numeric(n.trials)
    is.stop <- is.na(dt[1,])
    dt[1,!is.stop] <- n1PDFfixedt0.lnr(dt[,!is.stop,drop=F],meanlog,sdlog)
    # tmp <- try(integrate(f=stopfn,lower=0,upper=Inf,
    #     meanlog=meanlog,sdlog=sdlog)$value,silent=TRUE)
    # if (!is.numeric(tmp)) tmp <- 0
    tmp <- my.integrate(f=stopfn,lower=0,meanlog=meanlog,sdlog=sdlog)
    dt[1,is.na(dt[1,])] <- tmp 
    dt[1,]
  }
  
  n1PDF.lnrgng <- function(dt,meanlog,sdlog,t0,st0=0)
    # Same as n1PDF.lnr except NAs dealt with by integration
  {
    
    if ( st0 < 1e-3 ) # smaller values can cause numerical integration problems
      return(n1PDFfixedt0.lnrgng(meanlog=meanlog,sdlog=sdlog,dt=matrix(
        pmax(rep(dt,each=length(meanlog))-t0,0),nrow=length(meanlog)))
      ) else
      {    
        
        integrate.f <- function(dt,meanlog,sdlog,t0,st0) 
          n1PDFfixedt0.lnrgng(meanlog=meanlog,sdlog=sdlog,dt=matrix(pmax(
            rep(dt,each=length(meanlog))-t0,0),nrow=length(meanlog)))/st0
        
        outs <- numeric(length(dt))
        for (i in 1:length(outs)) {
          tmp <- try(integrate(f=integrate.f,
                               lower=dt[i]-st0[1],upper=dt[i],
                               meanlog=meanlog,sdlog=sdlog,t0=t0,st0=st0[1])$value,silent=TRUE)
          if (is.numeric(tmp)) 
            outs[i] <- tmp else
              outs[i] <- 0
        }
        return(outs)
      }
  }
  
  # # Check
  # 
  # # 2 choice 
  # n=1e5
  # meanlog=c(.5,.75); sdlog=c(1,1); t0=.2
  # sim <- rlnrgng(n=n,meanlog,sdlog,t0,st0=0)
  # dns <- plot.cell.density(sim,xlim=c(0,7),save.density=TRUE)
  # d <- n1PDF.lnrgng(dns$'2'$x,meanlog[c(2,1)],sdlog[c(2,1)],t0,st0=0)
  # lines(dns$'2'$x,d,col="red")
  # 
  # # p(Stop check)
  # mean(is.na(sim$RT))
  # n1PDFfixedt0.lnrgng(dt=matrix(rep(NA,2),ncol=1),meanlog,sdlog)
  # 
  # # 3 choice
  # n=1e5
  # meanlog=c(.5,.75,1); sdlog=c(1,1,1); t0=.2
  # sim <- rlnrgng(n=n,meanlog,sdlog,t0,st0=0)
  # dns <- plot.cell.density(sim,xlim=c(0,7),save.density=TRUE)
  # d <- n1PDF.lnrgng(dns$'2'$x,meanlog[c(2,1,3)],sdlog[c(2,1,3)],t0,st0=0)
  # lines(dns$'2'$x,d,col="red")
  # d <- n1PDF.lnrgng(dns$'3'$x,meanlog[c(3,2,1)],sdlog[c(3,2,1)],t0,st0=0)
  # lines(dns$'3'$x,d,col="red",lty=2)
  # 
  # # p(Stop check)
  # mean(is.na(sim$RT))
  # n1PDFfixedt0.lnrgng(dt=matrix(rep(NA,3),ncol=1),meanlog,sdlog)
  
  
