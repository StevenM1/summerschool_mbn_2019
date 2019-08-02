####################### ExGaussian stop signal----

  # WHY TRIAL SCALING IS BAD IN EXG, mean not proportional to sd
  # TRIALS=c(0:10)/10
  # mu=.5; sigma=.025; tau=.05; ts=.25
  # sigma=sigma+ts*TRIALS
  # tau=tau + ts*TRIALS
  # mn= mu + tau 
  # sd=sqrt( (sigma)^2 + (tau)^2 )
  # ratio=mn/sd
  # round(rbind(sigma,tau,mn,sd,ratio),3)
  # ratio[1]/ratio[10]
  
  
  # Modified from gamlss.dist to make cdf in nu > 0.05 * sigma case robust,
  # and robust to -Inf and Inf inputs, returns NA for bad sigma or tau and
  # robust against small sigma cases.
  pexGAUS <- function (q, mu = 5, sigma = 1, nu = 1, lower.tail = TRUE, log.p = FALSE) 
  {
    if (sigma <= 0) return(rep(NA,length(q)))
    if (nu <= 0) return(rep(NA,length(q)))
    
    # if (sigma < 0.05*nu) 
    if (sigma < 1e-4) 
      return(pexp(q-mu,1/nu,log=log.p,lower.tail=lower.tail)) # shfited exponential
    
    ly <- length(q)
    sigma <- rep(sigma, length = ly)
    mu <- rep(mu, length = ly)
    nu <- rep(nu, length = ly)
    index <- seq(along = q)
    z <- q - mu - ((sigma^2)/nu)
    cdf <- ifelse(is.finite(q), 
                  ifelse(nu > 0.05 * sigma, 
                         pnorm((q - mu)/sigma) - 
                           exp(log(pnorm(z/sigma)) + ((mu + (sigma^2/nu))^2 - (mu^2) - 
                                                        2 * q * ((sigma^2)/nu))/(2 * sigma^2)), 
                         pnorm(q, mean = mu, sd = sigma)),
                  ifelse(q<0,0,1)   
    )
    if (lower.tail == TRUE) 
      cdf <- cdf
    else cdf <- 1 - cdf
    if (log.p == FALSE) 
      cdf <- cdf
    else cdf <- log(cdf)
    cdf
  }
  
  # gamlss.dist function, but returns NA for bad sigma or tau, and
  # robust against small sigma cases.
  dexGAUS <- function (x, mu = 5, sigma = 1, nu = 1, log = FALSE) 
  {
    if (sigma <= 0) return(rep(NA,length(x)))
    if (nu <= 0) return(rep(NA,length(x)))
    
    # if (sigma < 0.05*nu) 
    if (sigma < 1e-4) 
      return(dexp(x-mu,1/nu,log=log)) # shfited exponential
    
    ly <- length(x)
    sigma <- rep(sigma, length = ly)
    mu <- rep(mu, length = ly)
    nu <- rep(nu, length = ly)
    z <- x - mu - ((sigma^2)/nu)
    logfy <- ifelse(nu > 0.05 * sigma, -log(nu) - (z + (sigma^2/(2 * 
                                                                   nu)))/nu + log(pnorm(z/sigma)), dnorm(x, mean = mu, sd = sigma, 
                                                                                                         log = TRUE))
    if (log == FALSE) 
      fy <- exp(logfy)
    else fy <- logfy
    fy
  }
  
  rexg <- function (mu,sigma,tau) 
    # NOTE: gamlss.dist rexGAUS is very slow so didnt use.
    # Ex-gaussian race, mu is a matrix with 1 row per accumulator, columns
    # are trials, ouput is a two column (RT, R) data frame, on row for each trial.
    # Sigma and tau can be vectors with one value, or one value for each accumulator 
    # or a matrix in the same form as mu. NB: CAN PRODUCE NEGATIVE RTS!
  {
    dt <- matrix(rnorm(n=length(mu), mean=mu,sd=sigma) +
                   rexp(n=length(mu),rate=1/tau),nrow = dim(mu)[1])
    winner <- apply(dt,2,which.min)    
    data.frame(RT=dt[cbind(winner,1:dim(mu)[2])],R=winner)
  }
  
  n1PDF.exg <- function(dt,mu,sigma,tau)
    # Generates defective PDF for responses on node= 1
    # dt (decison time) is a matrix with length(mu) rows, one row for
    # each accumulator to allow for different start times
  {
    
    dt[1,] <- dexGAUS(dt[1,],mu[1],sigma[1],tau[1])
    if (length(mu)>1) for (i in 2:length(mu))
      dt[1,] <- dt[1,]*pexGAUS(q=dt[i,],mu=mu[i],sigma=sigma[i],nu=tau[i],lower.tail=FALSE)
    dt[1,]
  }
  
  
  rexgss <- function (n, mu, sigma, tau, tf=0, gf=0, 
                      SSD=Inf, staircase=NA) 
    # Race among n accumulators, first of which is a stop accumulator.
    # NA returned for RT when winner = 1. Optional SSD arguement can be used to 
    # adjust mu for first accumulator to mu+SSD. SSD can be a scalar or vector 
    # length n. For trials with winning first accumulator RT and R set to NA. 
    # Adds SSD column to output. 
    # tf = trigger failure probability, gf = go failure probability
  {
    if ( length(SSD)==1 ) SSD <- rep(SSD,n)
    if ( any(is.na(SSD)) || length(SSD) != n )
      stop("SSD cannot have NAs and must be a scalar or same length as n")
    if ( !is.matrix(mu) ) 
      mu <- matrix(rep(mu,times=n),nrow=length(mu)) 
    if (gf > 0) # Setup for GO failure
      is.gf <- as.logical(rbinom(length(SSD),1,gf)) else
        is.gf <- logical(length(SSD))
      
      if ( all(!is.finite(SSD)) ) {              # ALL GO
        mu[1,] <- mu[1,] + SSD
        out <- rexg(mu,sigma,tau)
      } else {                                   # SOME STOP
        if ( any(is.na(staircase)) ) {           # STOP fixed SSD
          mu[1,] <- mu[1,] + SSD
          out <- rexg(mu,sigma,tau)
          if ( tf>0 ) {
            is.tf <- logical(length(SSD))
            is.tf[is.finite(SSD)][as.logical(rbinom(sum(is.finite(SSD)),1,tf))] <- TRUE  
            if ( any(is.tf) ) { 
              out[is.tf,] <- 
                rexg(mu[-1,is.tf,drop=FALSE],sigma[-1],tau[-1])
              out[is.tf,"R"] <- out[is.tf,"R"]+1
            }
          }
        } else {                                 # STOP, staircase
          if ( !is.numeric(staircase) | length(staircase)!=1 )
            stop("Staircase must be a numeric vector of length 1 specifying the step.")
          n_acc <- dim(mu)[1]
          dt <- matrix(rnorm(n=length(mu), mean=mu,sd=sigma) +
                         rexp(n=length(mu),rate=1/tau),nrow = dim(mu)[1])
          winner <- numeric(n)
          for (i in c(1:n)) {
            dt[1,i] <- dt[1,i] + SSD[i]
            if ( runif(1)<tf ) # Trigger failure
              winner[i] <- which.min(dt[2:n_acc,i])+1 else
                winner[i] <- which.min(dt[,i])
              if (is.gf[i]) winner[i] <- 1
              if (i!=n) if (winner[i]==1) 
                SSD[i+1] <- SSD[i] + staircase else
                  SSD[i+1] <- pmax(SSD[i] - staircase,0)
          }
          out <- data.frame(RT=dt[cbind(winner,1:n)],R=winner)
        }
      }
      out[out$R==1,"RT"] <- NA
      if (gf > 0) {
        out$RT[is.gf] <- NA
        out$R[is.gf] <- 1
      }
      cbind.data.frame(out,SSD=SSD)
  }
  
  
  n1PDF.exgss=function(dt,mu,sigma,tau,tf=0,gf=0,
                       SSD,Si) 
    # SSD is either a scalar or vector of length(dt)
    # stop accumulator must have name "NR"
    # SSD is subtracted from dt[Si,]
    # (i.e., stop accumulator RT) and dt=NA done by integration.
    # ts = slope of slowing (speeding if negative) over TRIALS, linear affect on 
    # mean and sd:  sigma + ts*TRIALS, tau + ts*trials (truncated to be positive)
    #
    # tf= probabiliy of trigger failure, where
    # L = trigger fail & respond + trigger and respond + trigger and no-response
    #   = tf*L(N-1)+(1-tf)[L(N)+p(S)],
    # L(N-1) = choice race likelihood (no stop accumulator), 
  # L(N) = full N unit race likelihood given they did respond, 
  # p(S) probability of stop winning
  #
  # gf = probabiliy of go failure. 
  # On go trials:   L = go fail (so no response) + go and L above 
  # L = gf + (1-gf)*[tf*L(N-1)+(1-tf)[L(N)+p(S)]]  or similarly
  # L =    [ p(non-response) ]    +           [ p(response) ] 
  #   = [ pf + (1-gf)(1-tf)p(S) ] + [ (1-gf){(tf*Ln(n-1) + (1-tf)*L(N))} ]
  {
    
    stopfn <- function(t,mu,sigma,tau,SSD,Si) 
    {
      dt <- matrix(rep(t+SSD,each=length(mu)),nrow=length(mu))
      dt[Si,] <- dt[Si,]-SSD
      i <- c(Si,c(1:length(mu))[-Si])
      n1PDF.exg(dt[i,,drop=FALSE],mu[i],sigma[i],tau[i])
    }
    
    dt=matrix(rep(dt,each=length(mu)),nrow=length(mu))
    
    is.stop <- is.na(dt[1,])  
    dt[Si,] <- dt[Si,] - SSD
    if ( any(!is.stop) ) 
    {
      if ( tf > 0 ) 
      {
        dt[1,!is.stop] <- (1-gf)*(
          tf*n1PDF.exg(dt[-Si,!is.stop,drop=FALSE],mu[-Si],sigma[-Si],tau[-Si]) +
            (1-tf)*n1PDF.exg(dt[,!is.stop,drop=FALSE],mu,sigma,tau))
      } else 
        dt[1,!is.stop] <- (1-gf)*n1PDF.exg(dt[,!is.stop,drop=FALSE],mu,sigma,tau)
    }
    if ( any(is.stop) ) for ( i in unique(SSD[is.stop]) ) {
      # tmp <- ifelse(!is.finite(i),0,
      #   try(integrate(f=stopfn,lower=-Inf,upper=Inf,
      #   mu=mu,sigma=sigma,tau=tau,SSD=i,Si=Si)$value,silent=TRUE))
      # if (!is.numeric(tmp)) tmp <- 0
      tmp <- ifelse(!is.finite(i),0,
                    my.integrate(f=stopfn,lower=-Inf,
                                 mu=mu,sigma=sigma,tau=tau,SSD=i,Si=Si))
      dt[1,is.stop & (SSD==i)] <- gf + (1-gf)*(1-tf)*tmp 
    }
    dt[1,]
  }
  
  
  # ########### TWO ACCUMULATOR CASE
  # 
  # # Check, two different SSDs 
  # n <- 1e5
  # SSD <- rep(c(1,10)/10,each=n/2)
  # TRIALS <- NA
  # 
  # # Stop and one go accumulator
  # mu=c(.75,1); sigma=c(.25,.25); tau=c(.25,.5)
  # 
  # ### RUN ONE OF THE FOLLOWING FOUR LINES
  # # Without trigger failure or go failure
  # tf=0; gf=0
  # # With trigger failure, no go failure
  # tf=.1;gf=0
  # # Without trigger failure, with go failure
  # tf=0; gf=.1
  # # With trigger failure and go failure
  # tf=.1;gf=.1
  # 
  # # Simulate data and get simlate stop rate
  # sim <- rexgss(n=n,mu,sigma,tau,SSD=SSD,tf=tf,gf=gf,TRIALS=TRIALS,ts=ts)
  # # Plot data
  # par(mfrow=c(1,2))
  # dns1 <- plot.cell.density(sim[sim$SSD==.1,],xlim=c(0,7),save.density=TRUE,main="SSD=.1")
  # dns2 <- plot.cell.density(sim[sim$SSD!=.1,],xlim=c(0,7),save.density=TRUE,main="SSD=1")
  # 
  # # Signal respond RT
  # dat <- sim; dat <- dat[!is.na(dat$RT),]; dat$R <- factor(as.character(dat$R))
  # round(tapply(dat$RT,dat[,c("R","SSD")],mean),2)
  # 
  # # Save simulated densities
  # x1 <- dns1$'2'$x; x2 <- dns2$'2'$x
  # SSD <- c(rep(.1,length(x1)),rep(1,length(x2)))
  # r1 <- c(2,1)
  # d.r1 <- n1PDF.exgss(dt=c(x1,x2),mu=mu[r1],sigma=sigma[r1],tau=tau[r1],
  #                     SSD=SSD,Si=2,tf=tf,gf=gf)
  # 
  # # Plot simulated (black) and theoretical (red) densities
  # par(mfrow=c(1,2))
  # # red=black?
  # plot(x1,dns1$'2'$y,type="l",main="SSD=.1",ylab="Density",xlab="RT",
  #      ylim=c(0,max(dns1$'2'$y)))
  # lines(x1,d.r1[1:length(x1)],col="red")
  # # red=black?
  # plot(x2,dns2$'2'$y,type="l",main="SSD=1",ylab="Density",xlab="RT",
  #      ylim=c(0,max(dns2$'2'$y)))
  # lines(x2,d.r1[(length(x2)+1):(2*length(x2))],col="red")
  # 
  # # p(Stop check)
  # tapply(is.na(sim$RT),sim$SSD,mean) # empirical
  # # Theoretical
  # n1PDF.exgss(NA,mu,sigma,tau,SSD=.1,Si=1,tf=tf,gf=gf)
  # n1PDF.exgss(NA,mu,sigma,tau,SSD=1,Si=1,tf=tf,gf=gf)
  # 
#   ########### THREE ACCUMULATOR CASE
#   
#   # Check, two different SSDs 
#   n=1e5
#   SSD = rep(c(1,10)/10,each=n/2)
#   
#   # Stop, first go accumulator correct, second error
#   mu=c(.75,1,1.25); sigma=c(.25,.25,.25); tau=c(.25,.5,.5)
#   
#   ### RUN ONE OF THE FOLLOWING FOUR LINES
#   # Without trigger failure or go failure
#   tf=0; gf=0
#   # With trigger failure, no go failure
#   tf=.1;gf=0
#   # Without trigger failure, with go failure
#   tf=0; gf=.1
#   # With trigger failure and go failure
#   tf=.1;gf=.1
#   
#   # Simulate data and get simlate stop rate
#   sim <- rexgss(n=n,mu,sigma,tau,SSD=SSD,tf=tf,gf=gf)
#   # Plot data
#   par(mfrow=c(1,2))
#   dns1 <- plot.cell.density(sim[sim$SSD==.1,],xlim=c(0,7),save.density=TRUE,main="SSD=.1")
#   dns2 <- plot.cell.density(sim[sim$SSD!=.1,],xlim=c(0,7),save.density=TRUE,main="SSD=1")
#   
#   # Signal respond RT
#   dat <- sim; dat <- dat[!is.na(dat$RT),]; dat$R <- factor(as.character(dat$R))
#   round(tapply(dat$RT,dat[,c("R","SSD")],mean),2)
#   
#   # Save simulated densities
#   x1c <- dns1$'2'$x; x2c <- dns2$'2'$x
#   x1e <- dns1$'3'$x; x2e <- dns2$'3'$x
#   
#   SSD <- c(rep(.1,length(x1c)),rep(1,length(x2c)))
#   r1 <- c(2,1,3)
#   d.r1 <- n1PDF.exgss(dt=c(x1c,x2c),mu=mu[r1],sigma=sigma[r1],tau=tau[r1],
#                       SSD=SSD,Si=2,tf=tf,gf=gf)
#   SSD <- c(rep(.1,length(x1e)),rep(1,length(x2e)))
#   r2 <- c(3,1,2)
#   d.r2 <- n1PDF.exgss(dt=c(x1e,x2e),mu=mu[r2],sigma=sigma[r2],tau=tau[r2],
#                       SSD=SSD,Si=2,tf=tf,gf=gf)
#   
#   # Plot simulated (black) and theoretical (red) densities
#   par(mfrow=c(1,2))
#   # red=black?
#   plot(x1c,dns1$'2'$y,type="l",main="SSD=.1",ylab="Density",xlab="RT",
#        ylim=c(0,max(dns1$'2'$y)))
#   lines(x1c,d.r1[1:length(x1c)],col="red")
#   lines(x1e,dns1$'3'$y,lty=2)
#   lines(x1e,d.r2[1:length(x1e)],col="red",lty=2)
#   # red=black?
#   plot(x2c,dns2$'2'$y,type="l",main="SSD=1",ylab="Density",xlab="RT",
#        ylim=c(0,max(dns2$'2'$y)))
#   lines(x2c,d.r1[(length(x2c)+1):(2*length(x2c))],col="red")
#   lines(x2e,dns2$'3'$y,lty=2)
#   lines(x2e,d.r2[(length(x2e)+1):(2*length(x2e))],col="red",lty=2)
#   
#   # p(Stop check)
#   tapply(is.na(sim$RT),sim$SSD,mean) # empirical
#   # Theoretical
#   n1PDF.exgss(NA,mu,sigma,tau,SSD=.1,Si=1,tf=tf,gf=gf)
#   n1PDF.exgss(NA,mu,sigma,tau,SSD=1,Si=1,tf=tf,gf=gf)


#   ########### FOUR ACCUMULATOR CASE
#   
#   # Check, two different SSDs 
#   n=1e5
#   SSD = rep(c(1,10)/10,each=n/2)
#   
#   # Stop, first go accumulator correct, second and third error
#   mu=c(.75,1,1.25,1.25); sigma=c(.25,.25,.25,.25); tau=c(.25,.5,.5,.5)
#   
#   ### RUN ONE OF THE FOLLOWING FOUR LINES
#   # Without trigger failure or go failure
#   tf=0; gf=0
#   # With trigger failure, no go failure
#   tf=.1;gf=0
#   # Without trigger failure, with go failure
#   tf=0; gf=.1
#   # With trigger failure and go failure
#   tf=.1;gf=.1
#   
#   # Simulate data and get simlate stop rate
#   sim <- rexgss(n=n,mu,sigma,tau,SSD=SSD,tf=tf,gf=gf)
#   # Plot data
#   par(mfrow=c(1,2))
#   dns1 <- plot.cell.density(sim[sim$SSD==.1,],xlim=c(0,7),save.density=TRUE,main="SSD=.1")
#   dns2 <- plot.cell.density(sim[sim$SSD!=.1,],xlim=c(0,7),save.density=TRUE,main="SSD=1")
#   
#   # Signal respond RT
#   dat <- sim; dat <- dat[!is.na(dat$RT),]; dat$R <- factor(as.character(dat$R))
#   round(tapply(dat$RT,dat[,c("R","SSD")],mean),2)
#   
#   # Save simulated densities
#   x1c <- dns1$'2'$x; x2c <- dns2$'2'$x
#   x1e1 <- dns1$'3'$x; x2e1 <- dns2$'3'$x
#   x1e2 <- dns1$'4'$x; x2e2 <- dns2$'4'$x
#   
#   SSD <- c(rep(.1,length(x1c)),rep(1,length(x2c)))
#   r1 <- c(2,1,3,4)
#   d.r1 <- n1PDF.exgss(dt=c(x1c,x2c),mu=mu[r1],sigma=sigma[r1],tau=tau[r1],
#                       SSD=SSD,Si=2,tf=tf,gf=gf)
#   SSD <- c(rep(.1,length(x1e1)),rep(1,length(x2e1)))
#   r2 <- c(3,1,2,4)
#   d.r2 <- n1PDF.exgss(dt=c(x1e1,x2e1),mu=mu[r2],sigma=sigma[r2],tau=tau[r2],
#                       SSD=SSD,Si=2,tf=tf,gf=gf)
#   SSD <- c(rep(.1,length(x1e2)),rep(1,length(x2e2)))
#   r2 <- c(4,1,2,3)
#   d.r3 <- n1PDF.exgss(dt=c(x1e2,x2e2),mu=mu[r2],sigma=sigma[r2],tau=tau[r2],
#                       SSD=SSD,Si=2,tf=tf,gf=gf)
#   
#   # Plot simulated (black) and theoretical (red) densities
#   par(mfrow=c(1,2))
#   # red=black?
#   plot(x1c,dns1$'2'$y,type="l",main="SSD=.1",ylab="Density",xlab="RT",
#        ylim=c(0,max(dns1$'2'$y)))
#   lines(x1c,d.r1[1:length(x1c)],col="red")
#   lines(x1e1,dns1$'3'$y,lty=2)
#   lines(x1e1,d.r2[1:length(x1e1)],col="red",lty=2)
#   lines(x1e2,dns1$'4'$y,lty=3)
#   lines(x1e2,d.r3[1:length(x1e2)],col="red",lty=3)
#   # red=black?
#   plot(x2c,dns2$'2'$y,type="l",main="SSD=1",ylab="Density",xlab="RT",
#        ylim=c(0,max(dns2$'2'$y)))
#   lines(x2c,d.r1[(length(x2c)+1):(2*length(x2c))],col="red")
#   lines(x2e1,dns2$'3'$y,lty=2)
#   lines(x2e1,d.r2[(length(x2e1)+1):(2*length(x2e1))],col="red",lty=2)
#   lines(x2e2,dns2$'4'$y,lty=3)
#   lines(x2e2,d.r3[(length(x2e2)+1):(2*length(x2e2))],col="red",lty=3)
#   
#   # p(Stop check)
#   tapply(is.na(sim$RT),sim$SSD,mean) # empirical
#   # Theoretical
#   n1PDF.exgss(NA,mu,sigma,tau,SSD=.1,Si=1,tf=tf,gf=gf)
#   n1PDF.exgss(NA,mu,sigma,tau,SSD=1,Si=1,tf=tf,gf=gf)
