######################## PNRWald ----

# n-choice postivie (zero truncated) normal rate trial to trial ~N(v,sv) Wald 
#    race, with t0, v, sv, a (boundary) parameterixaiton

# require("statmod") invgauss random, pdf, cdf etc.
# This has parameterizaiton mean (expected value) and dispersion = 1/shape. In 
# the accumulator interpritaiton dispersion is scaled relative to threshold (a)
# so without loss of generality we can set diffusion coefficeint to 1 and so
# shape = a^2. Using the notation v = mean rate, expected value (mean) = a/v 
# and varaince = mean^3 / a^2 = a / v^3

### Single accumulator model

rPNRWald <- function(n,a,v,sv,t0)
  # random function for single acumulator
{
  drifts <- matrix(rtnorm(n = n , mean = v, sd = sv, lower = 0), 
                   nrow = length(v))
  t0+rinvgauss(n,mean=a/drifts,shape=a^2)        
}


dPNRWald <- function(t,a,v,sv,posdrift=TRUE)
  # density for single accumulator
{
  # Convert to Desmond and Yang's notation
  d <- v/a # mean of delta
  l <- a^2 # lambda assuming diffusive variance is 1
  v <- sv^2 # variance of delta
  sqrt(l/(2*pi*t^3*(v*t+1)))*exp(-(l*(d*t-1)^2)/(2*t*(v*t+1)))*
    pnorm((v+d)*sqrt(l/(t*v^2+v)))/pnorm(d*sqrt(l/v)) # normalize
}


# # Check density and normalization
# n=1e6; a=3; v=2; sv=1
# integrate(dPNRWald,lower=0,upper=Inf,a=a,v=v,sv=sv)
# sim <- rPNRWald(n=n,a=a,v=v,sv=sv,t0=0)
# bad <- sim>10; mean(bad) # if this is large the check isnt valid
# dns <- density(sim[!bad])
# x <- dns$x[dns$x>0]
# d <- dPNRWald(x,a=a,v=v,sv=sv)
# plot(dns)
# lines(x,d,col="red")


pPNRWald <- function(t,a,v,sv)
  # cumulative density for single accumulator, NB: direct integration of
  # dPNRWald fails! This method using the analytic for the un-mixed
  # Wald CDF works much better
{
  if ( length(c(a,v,sv))>3 )
    stop("pPNRWald can only handle scalar parameters")
  ifun <- function(x,a,v,sv,t) {
    out <- dtnorm(x,v,sv,lower=0)*pinvgauss(t,mean=a/x,shape=a^2)
    ifelse(is.finite(out),out,0) 
  }
  
  #   for (i in 1:length(t)) {
  #     int <- try(integrate(ifun,lower=0,upper=Inf,a=a,v=v,sv=sv,t=t[i]),silent=TRUE)
  #     if (class(int)=="try-error") t[i] <- 0 else t[i] <- int$value
  #   }     
  
  uniquet <- unique(t) # scrimp some speed by calculating only for unqiue
  for (i in 1:length(uniquet)) {
    int <- try(integrate(ifun,lower=0,upper=Inf,a=a,v=v,sv=sv,t=uniquet[i]),silent=TRUE)
    if (class(int)=="try-error") 
      t[t==uniquet[i]] <- 0 else t[t==uniquet[i]] <- int$value
  }     
  
  t[t<0] <- 0; t[t>1] <- 1 
  t
}


# # Check cumulative density
# n=1e6; a=.1; v=.2; sv=1
# sim <- rPNRWald(n=n,a=a,v=v,sv=sv,t0=0)
# probs=1:990/1000
# qs <- quantile(sim,probs=probs)
# cd <- pPNRWald(qs,a=a,v=v,sv=sv)
# plot(qs,probs,type="l")
# lines(qs,cd,col="red")
# 
# # Speed up by using unique t
# sim <- rPNRWald(n=1e4,a=a,v=v,sv=sv,t0=0)
# system.time(pPNRWald(sim,a=a,v=v,sv=sv))
# #    user  system elapsed 
# #   7.720   0.009   7.731
# sim <- round(sim*1000)/1000
# system.time(pPNRWald(sim,a=a,v=v,sv=sv))
# #    user  system elapsed 
# #   0.875   0.000   0.875 

### Race model

rPNRWaldRace <- function(n,a,v,sv,t0,gf=0) 
  # random function for PNRWald race, if all accumualtors have non-poistive 
  # rates race deleted and warning returned. If only one has positve rate and sv
  # then do single accumulator special case, so can be accessed by setting v=0
  # and sv=0 for all other accumulators.
{
  n_v <- length(v)
  single_v <- c(1:length(v))[v>0] # Single accumulator case 
  if ( (length(single_v)==1) && 
       ( ( (length(sv)==1) && (sv<=0) ) || (all(sv[-single_v]<=0)) ) ) { 
    out <- data.frame(RT=rPNRWald(n,a[single_v],v[single_v],sv[single_v],t0[1]),
                      R=factor(rep(single_v,n),levels=1:n_v))  
  } else {
    drifts <- matrix(rtnorm(n=n*n_v,mean=v,sd=sv,lower=0),nrow=n_v) 
    ttf <- matrix(t0 + rinvgauss(n*n_v,mean=a/drifts,shape=a^2),nrow=n_v)
    resp <- apply(ttf, 2, which.min)
    out <- data.frame(RT = ttf[cbind(resp,1:n)], 
                      R = factor(apply(ttf, 2, which.min),levels=1:n_v))
  }
  if (gf[1] > 0) {
    is.gf <- as.logical(rbinom(dim(out)[1],1,gf))
    out$RT[is.gf] <- NA
    out$R[is.gf] <- 1
  }
  
  out
  
}


n1PNRWald <- function(dt,a,v,sv,t0,gf=0,posdrift=TRUE)
  # Generates defective PDF for responses on node=1dt (decison time) is a vector of times
  # If only one has positve rate and sv then do single accumulator special case, 
  # so can be accessed by setting v=0 and sv=0 for all other accumulators.
{
  # Some protection for negative dt values, remove if not needed
  dt <- dt-t0
  good <- is.go <- !is.na(dt)
  good[good] <- dt[good] > 0
  
  d <- numeric(length(dt))
  d[good] <- (1-gf[1])*dPNRWald(dt[good],a=a[1],v=v[1],sv=sv[1])
  single_v <- c(1:length(v))[v>0] 
  if ( !( (length(single_v)==1) && 
          ( ( (length(sv)==1) && (sv<=0) ) || (all(sv[-single_v]<=0)) )) )
    for (i in 2:length(v))
      d[good] <- d[good]*(1-pPNRWald(dt[good],a=a[i],v=v[i],sv=sv[i]))
  
  d[!is.go] <- gf[1]
  
  d
}


# # Check 
# n=1e5
# v=c(4,0.5); sv=c(1,1); a=c(2,2); t0=0.5
# sim <- rPNRWaldRace(n=n,a=a,v=v,sv=sv,t0=t0)
# 
# par(mfrow=c(1,2))
# dns <- plot.cell.density(sim,C=1,xlim=c(0,4),save.density=TRUE)
# # NB: Direclty integrating again appears very inaccurate
# integrate(n1PNRWald,lower=t0,upper=Inf,a=a,v=v,sv=sv,t0=t0)$value
# dt=dns$correct$x
# d <- n1PNRWald(dt,a=a,v=v,sv=sv,t0=t0)
# # red=black?
# plot(dns$correct$x,dns$correct$y,type="l")
# lines(dns$correct$x,d,col="red")
# 
# # With go failure
# n=1e5
# v=c(2,0.5); sv=c(1,1); a=c(2,2); t0=0.5; gf=.2
# sim <- rPNRWaldRace(n=n,a=a,v=v,sv=sv,t0=t0,gf=gf)
# par(mfrow=c(1,3))
# dns <- plot.cell.density(sim,C=1,xlim=c(0,4),save.density=TRUE)
# # NB: Direclty integrating again appears very inaccurate
# integrate(n1PNRWald,lower=t0,upper=Inf,a=a,v=v,sv=sv,t0=t0,gf=gf)$value/(1-gf)
# dt=dns$correct$x
# d <- n1PNRWald(dt,a=a,v=v,sv=sv,t0=t0,gf=gf)
# # red=black?
# plot(dns$correct$x,dns$correct$y,type="l")
# lines(dns$correct$x,d,col="red")
# dt=dns$error$x
# d <- n1PNRWald(dt,a=a[2:1],v=v[2:1],sv=sv[2:1],t0=t0,gf=gf)
# # red=black?
# plot(dns$error$x,dns$error$y,type="l")
# lines(dns$error$x,d,col="red")

# # Single accumlator with go failure
# n=1e5
# v=c(3,0); sv=c(1,0); a=c(2,2); t0=0.5; gf=.2
# sim <- rPNRWaldRace(n=n,a=a,v=v,sv=sv,t0=t0,gf=gf)
# par(mfrow=c(1,3))
# dns <- plot.cell.density(sim,C=1,xlim=c(0,4),save.density=TRUE)
# # NB: Direclty integrating again appears very inaccurate
# integrate(n1PNRWald,lower=t0,upper=Inf,a=a,v=v,sv=sv,t0=t0,gf=gf)$value/(1-gf)
# dt=dns$correct$x
# d <- n1PNRWald(dt,a=a,v=v,sv=sv,t0=t0,gf=gf)
# # red=black?
# plot(dns$correct$x,dns$correct$y,type="l")
# lines(dns$correct$x,d,col="red")
# 
# n=1e5
# v=c(0,3); sv=c(0,1); a=c(2,2); t0=0.5; gf=.2
# sim <- rPNRWaldRace(n=n,a=a,v=v,sv=sv,t0=t0,gf=gf)
# par(mfrow=c(1,3))
# dns <- plot.cell.density(sim,C=1,xlim=c(0,4),save.density=TRUE)
# # NB: Direclty integrating again appears very inaccurate
# 1-integrate(n1PNRWald,lower=t0,upper=Inf,a=a[2:1],v=v[2:1],sv=sv[2:1],t0=t0,gf=gf)$value/(1-gf)
# dt=dns$error$x
# d <- n1PNRWald(dt,a=a[2:1],v=v[2:1],sv=sv[2:1],t0=t0,gf=gf)
# # red=black?
# plot(dns$error$x,dns$error$y,type="l")
# lines(dns$error$x,d,col="red")

######################## Wald ----

# n-choice uniformly varying start point (0-A) Wald 
#    race, with t0, v, A, b (boundary) parameterixaiton

### Single accumulator model

# pigt, digt, rwaldt Copyright (C) 2013  Trisha Van Zandt distributed with: 
# Logan, Van Zandt, Verbruggen, and Wagenmakers (2014).  On the ability to 
# inhibit thought and action: General and special theories of an act of control.
# Psychological Review. Comments and changes added by Andrew Heathcote. Trish's
# code is for k = threshold, a = half width of uniform threshold variability,
# l = rate of accumulation. Note that Wald mean = k/l and shape = k^2.

# Following functions use a different parameterization in terms of v=l (rate),
# uniform start point variability from 0-A (A>=0), threshold b (>0) and hence 
# B=b-A (>=0) as a threshold gap. Hence k = b-A/2 = B + A/2 and a=A 

rWald <- function(n,B,v,A)
  # random function for single acumulator
{
  
  rwaldt <- function(n,k,l,tiny=1e-6) {
    # random sample of n from a Wald (or Inverse Gaussian)
    # k = criterion, l = rate, assumes sigma=1 Browninan motion
    # about same speed as statmod rinvgauss
    
    rlevy <- function(n=1, m=0, c=1) {
      if (any(c<0)) stop("c must be positive")
      c/qnorm(1-runif(n)/2)^2+m
    }
    
    flag <- l>tiny
    x <- rep(NA,times=n)
    
    x[!flag] <- rlevy(sum(!flag),0,k[!flag]^2)
    mu <- k/l
    lambda <- k^2
    
    y <- rnorm(sum(flag))^2
    mu.0 <- mu[flag]
    lambda.0 <- lambda[flag]
    
    x.0 <- mu.0 + mu.0^2*y/(2*lambda.0) -
      sqrt(4*mu.0*lambda.0*y + mu.0^2*y^2)*mu.0/(2*lambda.0)
    
    z <- runif(length(x.0))
    test <- mu.0/(mu.0+x.0)
    x.0[z>test] <- mu.0[z>test]^2/x.0[z>test]
    x[flag] <- x.0
    x[x<0] <- max(x)
    x
  }
  
  # Act as if negative v never terminates, cluge to do single accumulator
  # case by passing negative v
  if (length(v)!=n) v <- rep(v,length.out=n)
  if (length(B)!=n) B <- rep(B,length.out=n)
  if (length(A)!=n) A <- rep(A,length.out=n)
  
  # Kluge to return -Inf for negative rates, so can implment one accumulator case
  out <- numeric(n)
  ok <- v>0  
  nok <- sum(ok)
  bs <- B[ok]+runif(nok,0,A[ok])
  out[ok] <- rwaldt(nok,k=bs,l=v[ok])
  out[!ok] <- Inf
  out
}


dWald <- function(t,v,B,A)
  # density for single accumulator
{
  
  digt <- function(t,k=1,l=1,a=.1,tiny=1e-10) {
    # pdf of inverse gaussian at t with k +/- a/2 uniform variability
    # returns digt.0 if a<1e-10 
    
    digt.0 <- function(t,k=1,l=1) {
      # pdf of inverse gaussian at t with no k variability
      # much faster than statmod's dinvgauss funciton
      
      lambda <- k^2
      l0 <- l==0
      e <- numeric(length(t))
      if ( any(!l0) ) {
        mu <- k[!l0]/l[!l0]
        e[!l0] <- -(lambda[!l0]/(2*t[!l0])) * (t[!l0]^2/mu^2 - 2*t[!l0]/mu  + 1)
      }
      if ( any(l0) )  e[l0] <- -.5*lambda[l0]/t[l0]
      x <- exp(e + .5*log(lambda) - .5*log(2*t^3*pi))
      x[t<=0] <- 0
      x
    }
    
    options(warn=-1)
    if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
    if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
    if(length(a)!=length(t)) a <- rep(a,length.out=length(t))
    
    tpos <- t<=0
    
    atiny <- a<=tiny & !tpos
    a[atiny] <- 0
    
    ltiny <- (l<=tiny) & !atiny & !tpos
    notltiny <- (l>tiny) & !atiny & !tpos
    l[l<=tiny] <- 0
    
    x <- numeric(length(t))
    
    # No threshold variability
    if ( any(atiny) )
      x[atiny] <- digt.0(t=t[atiny],k=k[atiny],l=l[atiny])
    
    # Threshold variability
    if ( any(!atiny) ) {
      
      if ( any(notltiny) ) { # rate non-zero
        
        sqr.t <- sqrt(t[notltiny])
        
        term.1a <- -(a[notltiny]-k[notltiny]+t[notltiny]*l[notltiny])^2/(2*t[notltiny])
        term.1b <- -(a[notltiny]+k[notltiny]-t[notltiny]*l[notltiny])^2/(2*t[notltiny])
        term.1 <- (exp(term.1a) - exp(term.1b))/sqrt(2*pi*t[notltiny])
        
        term.2a <- log(.5)+log(l[notltiny])
        term.2b <- 2*pnorm((-k[notltiny]+a[notltiny])/sqr.t+sqr.t*l[notltiny])-1
        term.2c <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
        term.2d <- term.2b+term.2c
        term.2 <- exp(term.2a)*term.2d
        
        term.3 <- term.1+term.2
        term.4 <- log(term.3)-log(2)-log(a[notltiny])
        x[notltiny] <- exp(term.4)
      }
      
      if ( any(ltiny) ) {  # rate zero
        log.t <- log(t[ltiny])
        term.1 <- -.5*(log(2)+log(pi)+log.t)
        term.2 <- (k[ltiny]-a[ltiny])^2/(2*t[ltiny])
        term.3 <- (k[ltiny]+a[ltiny])^2/(2*t[ltiny])
        term.4 <- (exp(-term.2)-exp(-term.3))
        term.5 <- term.1+log(term.4) - log(2) - log(a[ltiny])
        x[ltiny] <- exp(term.5)
      }
      
    }
    
    x[x<0 | is.nan(x) ] <- 0
    x
  }
  
  out <- numeric(length(t))
  ok <- v>0
  out[ok] <- digt(t[ok],k=B[ok]+A[ok]/2,l=v[ok],a=A[ok]/2)
  out[!ok] <- 0
  out
}


# Check density and normalization
# par(mfrow=c(1,2))
# n=1e6; A=0; v=2; B=1
# integrate(dWald,lower=0,upper=Inf,A=A,v=v,B=B)
# sim <- rWald(n=n,A=A,v=v,B=B)
# bad <- sim>10; mean(bad) # if this is large the check isnt valid
# dns <- density(sim[!bad])
# x <- dns$x[dns$x>0]
# d <- dWald(x,A=A,v=v,B=B)
# plot(dns)
# lines(x,d,col="red")
# 
# n=1e6; A=1; v=2; B=1
# integrate(dWald,lower=0,upper=Inf,A=A,v=v,B=B)
# sim <- rWald(n=n,A=A,v=v,B=B)
# bad <- sim>10; mean(bad) # if this is large the check isnt valid
# dns <- density(sim[!bad])
# x <- dns$x[dns$x>0]
# d <- dWald(x,A=A,v=v,B=B)
# plot(dns)
# lines(x,d,col="red")


pWald <- function(t,v,B,A)
  # cumulative density for single accumulator
{
  pigt <- function(t,k=1,l=1,a=.1,tiny=1e-10) {
    # cdf of inverse gaussian at t with k +/- a/2 uniform variability
    # returns pigt.0 if a<=0
    
    pigt.0 <- function(t,k=1,l=1) {
      # cdf of inverse gaussian at t with no k variability
      # much faster than statmod's pinvgauss funciton
      
      mu <- k/l
      lambda <- k^2
      
      e <- exp(log(2*lambda) - log(mu))
      add <- sqrt(lambda/t) * (1 + t/mu)
      sub <- sqrt(lambda/t) * (1 - t/mu)
      
      p.1 <- 1 - pnorm(add)
      p.2 <- 1 - pnorm(sub)
      x <- exp(e + log(p.1)) + p.2
      
      x[t<0] <- 0
      x
    }
    
    options(warn=-1)
    if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
    if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
    if(length(a)!=length(t)) a <- rep(a,length.out=length(t))
    
    tpos <- t<=0
    
    atiny <- a<=tiny & !tpos
    a[atiny] <- 0
    
    ltiny <- (l<=tiny) & !atiny & !tpos
    notltiny <- (l>tiny) & !atiny & !tpos
    l[l<=tiny] <- 0
    
    x <- numeric(length(t))
    
    # No threshold variability
    if ( any(atiny) )
      x[atiny] <- pigt.0(t[atiny],k[atiny],l[atiny])
    
    # Threshold variability
    if ( any(!atiny) ) {
      
      if ( any(notltiny) ) { # rate non-zero
        
        log.t <- log(t[notltiny])
        sqr.t <- sqrt(t[notltiny])
        
        term.1a <- .5*log.t-.5*log(2*pi)
        term.1b <- exp(-((k[notltiny]-a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
        term.1c <- exp(-((k[notltiny]+a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
        term.1 <- exp(term.1a)*(term.1b-term.1c)
        
        term.2a <- exp(2*l[notltiny]*(k[notltiny]-a[notltiny]) + 
                         log(pnorm(-(k[notltiny]-a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
        term.2b <- exp(2*l[notltiny]*(k[notltiny]+a[notltiny]) + 
                         log(pnorm(-(k[notltiny]+a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
        term.2 <- a[notltiny] + (term.2b-term.2a)/(2*l[notltiny])
        
        term.4a <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
        term.4b <- 2*pnorm((k[notltiny]-a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
        term.4c <- .5*(t[notltiny]*l[notltiny] - a[notltiny] - k[notltiny] + .5/l[notltiny])
        term.4d <- .5*(k[notltiny] - a[notltiny] - t[notltiny]*l[notltiny] - .5/l[notltiny])
        term.4 <- term.4c*term.4a + term.4d*term.4b
        
        x[notltiny] <- (term.4 + term.2 + term.1)/(2*a[notltiny])
      }
      
      if ( any(ltiny) ) {  # rate zero
        sqr.t <- sqrt(t[ltiny])
        log.t <- log(t[ltiny])
        term.5a <- 2*pnorm((k[ltiny]+a[ltiny])/sqr.t)-1
        term.5b <- 2*pnorm(-(k[ltiny]-a[ltiny])/sqr.t)-1
        term.5 <- (-(k[ltiny]+a[ltiny])*term.5a - (k[ltiny]-a[ltiny])*term.5b)/(2*a[ltiny])
        
        term.6a <- -.5*(k[ltiny]+a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
        term.6b <- -.5*(k[ltiny]-a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
        term.6 <- 1 + exp(term.6b) - exp(term.6a)
        
        x[ltiny] <- term.5 + term.6
      }
      
    }
    
    x[x<0 | is.nan(x) ] <- 0
    x
  }
  
  out <- numeric(length(t))
  ok <- v>0
  out[ok] <- pigt(t[ok],k=B[ok]+A[ok]/2,l=v[ok],a=A[ok]/2)
  out[!ok] <- 0
  out
  
}


# # Check cumulative density
# par(mfrow=c(1,2))
# 
# n=1e6; A=0; v=2; B=1
# sim <- rWald(n=n,A=A,v=v,B=B)
# probs=1:990/1000
# qs <- quantile(sim,probs=probs)
# cd <- pWald(qs,A=A,v=v,B=B)
# plot(qs,probs,type="l")
# lines(qs,cd,col="red")
# 
# n=1e6; A=1; v=2; B=1
# sim <- rWald(n=n,A=A,v=v,B=B)
# probs=1:990/1000
# qs <- quantile(sim,probs=probs)
# cd <- pWald(qs,A=A,v=v,B=B)
# plot(qs,probs,type="l")
# lines(qs,cd,col="red")

### Race model

rWaldRace <- function(n,v,B,A,t0,gf=0,return.ttf=FALSE) 
  # random function for Wald race.
{
  B[B<0] <- 0 # Protection for negatives 
  A[A<0] <- 0
  n_v  <- ifelse(is.null(dim(v)), length(v), dim(v)[1])
  ttf <- matrix(t0 + rWald(n*n_v,B=B,v=v,A=A),nrow=n_v)
  if (return.ttf) return(ttf)
  resp <- apply(ttf, 2, which.min)
  out <- data.frame(RT = ttf[cbind(resp,1:n)], R = apply(ttf, 2, which.min))
  
  if (gf[1] > 0) {
    is.gf <- as.logical(rbinom(dim(out)[1],1,gf))
    out$RT[is.gf] <- NA
    out$R[is.gf] <- 1
  }
  
  out
}


n1Wald <- function(dt,v,B,A,t0,gf=0)
  # Generates defective PDF for responses on node=1, dt (decison time) is a vector of times
{
  B[B<0] <- 0 # Protection for negatives 
  A[A<0] <- 0
  n_acc <- ifelse(is.null(dim(v)),length(v),dim(v)[1])
  if (is.null(dim(dt))) dt <- matrix(rep(dt,each=n_acc),nrow=n_acc)
  dt <- dt-t0
  
  is.go <- !is.na(dt[1,])
  n.go <- sum(is.go)
  
  if (!is.matrix(v)) v <- matrix(rep(v,n.go),nrow=n_acc)
  if (!is.matrix(B)) B <- matrix(rep(B,n.go),nrow=n_acc)
  if (!is.matrix(A)) A <- matrix(rep(A,n.go),nrow=n_acc)
  
  # Winner
  dt[1,is.go] <- (1-gf[1])*dWald(dt[1,is.go],A=A[1,],v=v[1,],B=B[1,])
  if (n_acc > 1) for (i in 2:n_acc)
    dt[1,is.go] <- dt[1,is.go]*(1-pWald(dt[i,is.go],A=A[i,],v=v[i,],B=B[i,]))
  
  dt[1,!is.go] <- gf[1]
  
  dt[1,]
}


# # Check no go failure 
# n=1e5
# v=c(2,1); B=c(1,1); A=c(2,2); t0=1
# sim <- rWaldRace(n=n,A=A,v=v,B=B,t0=t0)
# par(mfrow=c(1,3))
# dns <- plot.cell.density(sim,C=1,xlim=c(0,4),save.density=TRUE)
# pc <- integrate(n1Wald,lower=0,upper=Inf,A=A,v=v,B=B,t0=t0)$value
# print(pc)
# dt=dns$correct$x
# d <- n1Wald(dt,A=A,v=v,B=B,t0=t0)
# # red=black?
# plot(dns$correct$x,dns$correct$y,type="l")
# lines(dns$correct$x,d,col="red")
# dt=dns$error$x
# d <- n1Wald(dt,A=A[2:1],v=v[2:1],B=B[2:1],t0=t0)
# plot(dns$error$x,dns$error$y,lty=2,type="l")
# lines(dns$error$x,d,col="red",lty=2)
# 
# # Check go failure 
# n=1e5
# v=c(2,1); B=c(1,1); A=c(2,2); t0=1; gf=.2
# sim <- rWaldRace(n=n,A=A,v=v,B=B,t0=t0,gf=gf)
# par(mfrow=c(1,3))
# dns <- plot.cell.density(sim,C=1,xlim=c(0,4),save.density=TRUE)
# pc <- integrate(n1Wald,lower=0,upper=Inf,A=A,v=v,B=B,t0=t0)$value
# print(pc)
# dt=dns$correct$x
# d <- n1Wald(dt,A=A,v=v,B=B,gf=gf,t0=t0)
# # red=black?
# plot(dns$correct$x,dns$correct$y,type="l")
# lines(dns$correct$x,d,col="red")
# dt=dns$error$x
# d <- n1Wald(dt,A=A[2:1],v=v[2:1],B=B[2:1],gf=gf,t0=t0)
# plot(dns$error$x,dns$error$y,lty=2,type="l")
# lines(dns$error$x,d,col="red",lty=2)

# # Check (effectively) single accumulator case 
# n=1e5
# v=c(2,-1); B=c(1,1); A=c(2,2); t0=1
# sim <- rWaldRace(n=n,A=A,v=v,B=B,t0=t0)
# par(mfrow=c(1,4))
# dns <- plot.cell.density(sim,C=1,xlim=c(0,4),save.density=TRUE)
# pc <- integrate(n1Wald,lower=0,upper=Inf,A=A,v=v,B=B,t0=t0)$value
# print(pc)
# dt=dns$correct$x
# d <- n1Wald(dt,A=A,v=v,B=B,t0=t0)
# # red=black?
# plot(dns$correct$x,dns$correct$y,type="l")
# lines(dns$correct$x,d,col="red")
# 
# v=c(-1,2); B=c(1,1); A=c(2,2); t0=1
# sim <- rWaldRace(n=n,A=A,v=v,B=B,t0=t0)
# dns <- plot.cell.density(sim,C=1,xlim=c(0,4),save.density=TRUE)
# pe <- integrate(n1Wald,lower=0,upper=Inf,A=A,v=v,B=B,t0=t0)$value
# print(pe)
# dt=dns$error$x
# d <- n1Wald(dt,A=A[2:1],v=v[2:1],B=B[2:1],gf=gf,t0=t0)
# # red=black?
# plot(dns$error$x,dns$error$y,type="l")
# lines(dns$error$x,d,col="red")