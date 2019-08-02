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
    if (st0[1]==0) data.frame(rt=dt[cbind(winner,1:n)],response=winner) else
      data.frame(rt=dt[cbind(winner,1:n)]+runif(n,0,st0[1]),response=winner)
  }
