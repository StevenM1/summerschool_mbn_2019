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
  
  n1PDF.lnr <- function(dt,meanlog,sdlog,t0,st0=0)
    # dt (decision time) is a vector, meanlog and sdlog have same length = number of 
    # accumulators, t0 is the lower bound of non-decision time, it can be:
    # a) a scalar, b) a vector of length number of accumulators or
    # c) a matrix with 1 row per accumulator, when start times differ on each trial
    # st0, range of non-decison time variability, must be a scalar, as the same
    # variability is assumed in a common encoding/production stage
  {
    
    # NOTE: No need to flag negative dt in dt-t0 as plnorm/dlnorm return 0
    
    n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])
    if (!is.matrix(meanlog)) meanlog <- matrix(rep(meanlog,length(dt)),nrow=n_acc)
    if (!is.matrix(sdlog))     sdlog <- matrix(rep(  sdlog,length(dt)),nrow=n_acc)
    if ( st0 < 1e-3 ) # smaller values can cause numerical integration problems
      return(n1PDFfixedt0.lnr(meanlog=meanlog,sdlog=sdlog,dt=matrix(
        pmax(rep(dt,each=n_acc)-t0,0),nrow=n_acc))) else
        {    
          
          integrate.f <- function(dt,meanlog,sdlog,t0,st0) 
            n1PDFfixedt0.lnr(meanlog=meanlog,sdlog=sdlog,dt=matrix(
              pmax(rep(dt,each=length(meanlog))-t0,0),nrow=length(meanlog)))/st0
          
          outs <- numeric(length(dt))
          for (i in 1:length(outs)) {
            tmp <- try(integrate(f=integrate.f,
                                 lower=dt[i]-st0[1],upper=dt[i],
                                 meanlog=meanlog[,i],sdlog=sdlog[,i],t0=t0,st0=st0[1])$value,silent=TRUE)
            if (is.numeric(tmp)) 
              outs[i] <- tmp else
                outs[i] <- 0
          }
          return(outs)
        }
  }
  
  
  # # Check
  # n=1e5
  # meanlog=c(.5,.75,1); sdlog=c(1,1,1)
  # n_acc <- length(meanlog)
  # 
  # # check scalar t0
  # t0=.2
  # t0=c(.2,1,1)
  # 
  # # check t0 noise
  # st0=1
  # st0=0
  # 
  # sim <- rlnr(n=n,meanlog,sdlog,t0,st0)
  # dns <- plot.cell.density(sim,xlim=c(0,7),save.density=TRUE)
  # 
  # # Check matrix forms work
  # meanlog <- matrix(rep(meanlog,length(dns[[ichar]]$x)),nrow=n_acc)
  # sdlog <- matrix(rep(sdlog,length(dns[[ichar]]$x)),nrow=n_acc)
  # t0 <- matrix(rep(t0,length.out=n_acc*length(dns[[ichar]]$x)),nrow=n_acc)
  # 
  # check_n1 <- TRUE   # n1PDF check
  # check_n1 <- FALSE  # n1PDFfixedt0.lnr check (ignores st0)
  # 
  # for (i in 1:n_acc) {
  #   ichar <- as.character(i)
  #   indx <- c(i,c(1:n_acc)[-i])
  #   if (!is.matrix(meanlog)) meanlogi <- meanlog[indx] else meanlogi <- meanlog[indx,]
  #   if (!is.matrix(sdlog))     sdlogi <- sdlog[indx] else     sdlogi <- sdlog[indx,]
  #   if (!is.matrix(t0)) {
  #     if (length(t0)!=1) t0i <- t0[indx] else t0i <- t0 
  #   } else t0i <- t0[indx,]
  #   if (check_n1) d <- n1PDF.lnr(dt=dns[[ichar]]$x,meanlogi,sdlogi,t0i,st0=st0) else
  #     d <- n1PDFfixedt0.lnr(meanlog=meanlogi,sdlog=sdlogi,
  #       dt=matrix(pmax(rep(dns[[ichar]]$x,each=n_acc)-t0i,0),nrow=n_acc))
  #   if (i!=1) lines(dns[[ichar]]$x,dns[[ichar]]$y,lty=i) else
  #     plot(dns[[ichar]]$x,dns[[ichar]]$y,type="l",xlab="RT",ylab="Density") 
  #   lines(dns[[ichar]]$x,d,col="red",lty=i)
  # }
