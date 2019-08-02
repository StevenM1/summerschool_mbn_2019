require(rtdists) 


##These functions could be replaced by rlba_norm in rtdists.

##Simulate LBA trials

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


##Simulate PLBA trials

# A=c(1.5,1.5); B=c(1.2,1.2); C=c(0.3,0.1)
# v=c(3.32,2.24); w=c(1.51,3.69); sv=c(1,1); sw=c(1,1)
# rD=0.3; tD=.3; swt=0.5; t0=0.08

rplba <- function(n, A=c(1.5,1.5), B=c(1.2,1.2), C=c(0.3,0.3),
  v=c(3.32,2.24), w=c(1.51,3.69), sv=c(1,1),sw=c(1,1), 
  rD=0.3, tD=.3, swt=0.5, t0=0.08) 
  
{
  # Stage 1 LBA
  v1 <- rtnorm(n, v[1], sv[1], 0)
  v2 <- rtnorm(n, v[2], sv[2], 0)
  sp <- matrix(runif(2*n,0,A),nrow=2)
  
  # Calcualte thresholds
  b1 <- A[1]+B[1]
  b2 <- A[2]+B[2]
  c1 <- b1 + C[1]
  c2 <- b2 + C[2]
  
  # Race
  dt <- rbind((b1-sp[1,])/v1,(b2-sp[2,])/v2)
  # dt[dt<0] <- Inf
  choice <- apply(dt,2,which.min)
  rt <- dt[cbind(choice,1:n)]
  
  # Calculate effective switch times
  swt_b <- swt[1] + tD[1] # Threshold delay
  swt_r <- swt[1] + rD[1] # Rate delay
  
  # Which switch is first
  swt <- pmin(swt_b,swt_r)
  if (swt_b==swt_r) {
    change <- "both"
  } else if (swt_r < swt_b) {
    change <- "rate"
  } else { 
    change <- "threshold"
  }
  
  # Which are finished?
  done <- rt <= swt
  n2 <- sum(!done)
  
  # Stage 2 LBA
  
  # Distance left to travel for those not finished
  # threshold - distance already travelled
  if ( change=="rate" ) {
    B1 <- b1 - (sp[1,!done] + swt*v1[!done]) 
    B2 <- b2 - (sp[2,!done] + swt*v2[!done]) 
  } else {
    B1 <- c1 - (sp[1,!done] + swt*v1[!done]) 
    B2 <- c2 - (sp[2,!done] + swt*v2[!done]) 
  }
  
  # Change rates?
  if ( change=="threshold" ) {
    w1 <- v1[!done]; w2 <- v2[!done]
  } else {
    w1 <- rtnorm(n2,w[1],sw[1],0)
    w2 <- rtnorm(n2,w[2],sw[2],0)
  } 
  
  # Race  
  dt <- rbind(B1/w1,B2/w2)
  # dt[dt<0] <- Inf
  choice[!done] <- apply(dt,2,which.min)
  rt[!done] <- swt+dt[cbind(choice[!done],1:n2)]
  
  if ( change != "both" ) { # Stage 3 LBA
    
    if ( change=="threshold" ) swt1 <- swt_r else swt1 <- swt_b
    t2 <- swt1-swt
    
    # Which are finished?
    done1 <- rt[!done] < swt1
    n2 <- sum(!done1)
    
    if ( !all(done1) ) {
      
      # Distance left to travel for those not finished
      # Distance left at end of stage 1 - further travel
      B1 <- B1[!done1] - t2*w1[!done1]
      B2 <- B2[!done1] - t2*w2[!done1]
      
      if ( change=="threshold" ) {
        w1 <- rtnorm(n2,w[1],sw[1],0)
        w2 <- rtnorm(n2,w[2],sw[2],0)
      }  else {
        w1 <- w1[!done1]; w2 <- w2[!done1]
        B1 <- B1 + C[1]
        B2 <- B2 + C[2]
      }
      
      # Race  
      dt <- rbind(B1/w1,B2/w2)
      # dt[dt<0] <- Inf
      choice[!done][!done1] <- apply(dt,2,which.min)
      rt[!done][!done1] <- swt1+dt[cbind(choice[!done][!done1],1:n2)]
    }
    
  }
  
  # save results NOTE: columns must be in this order!
  cbind(rt=t0[1]+rt,choice=choice)
}

