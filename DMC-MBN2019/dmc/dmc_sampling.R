# System functions for the DMC (Dynamic Models of Choice)
#    Functions specific to sampling and related operations
#    Usually user does not need to edit

######### PRIOR and POSTERIOR ----
 
dbeta_lu <- function(x,shape1,shape2,lower,upper,log=FALSE) 
  # Used with beta prior
{
  if (!log) dbeta((x-lower)/(upper-lower),shape1,shape2,log=FALSE)/(upper-lower) else 
    dbeta((x-lower)/(upper-lower),shape1,shape2,log=TRUE) - log(upper-lower) 
}

rbeta_lu <- function(n,shape1,shape2,lower,upper) 
  # Used with beta prior
{
  lower + rbeta(n,shape1,shape2)*(upper-lower)    
}

dgamma_l <- function(x,shape,scale,lower,log=FALSE) 
  # Used with gamma prior
{
  dgamma(x-lower,shape=shape,scale=scale,log=log)  
}

rgamma_l <- function(n,shape,scale,lower) 
  # Used with gamma prior
{
  lower + rgamma(n,shape=shape,scale=scale)    
}

dlnorm_l <- function(x,meanlog,sdlog,lower,log=FALSE) 
  # Used with lognormal prior
{
  dlnorm(x-lower,meanlog,sdlog,log=log)  
}

rlnorm_l <- function(n,meanlog,sdlog,lower) 
  # Used with lognormal prior
{
  lower + rlnorm(n,meanlog,sdlog)    
}


dcauchy <- function(x,location,scale,log=FALSE) 
  # Used with cauchy prior
{
  dcauchy(x,location,scale,log=log)  
}

rcauchy <- function(n,location,scale) 
  # Used with cauchy prior
{
  rcauchy(n,location,scale)    
}


dtt <- function(x, location, scale, df, lower, upper, log = FALSE)
  # Used with truncated t prior
{
  logdens <- ifelse(x >= lower & x <= upper,
                    dt((x - location) / scale, df = df, log = TRUE) - log(scale),
                    -Inf)
  logconstant <- log(pt((upper - location) / scale, df = df) -
                       pt((lower - location) / scale, df = df))
  out <- logdens - logconstant
  
  if ( ! log) {
    out <- exp(out)
  }
  
  return(out)
}


rtt <- function(n,location,scale,df,lower,upper,log=FALSE)
  # Used with truncated t prior
{
  u <- runif(n)
  p <- pt((lower - location) / scale, df = df) +
    u * (pt((upper - location) / scale, df = df) -
         pt((lower - location) / scale, df = df))
  out <- qt(p = p, df = df) * scale + location
  return(out)
}


dconstant <- function(x,constant,log=FALSE)
  # Used with constant prior
{
  if (log) rep(0,length(constant)) else
           rep(1,length(constant))
}


rconstant <- function(n,constant)
  # Used with constant prior
{
  rep(constant,n)
}

  
prior.p.dmc <- function(p1,p2,p3=rep(NA,length(p1)),
                        lower=rep(NA,length(p1)),
                        upper=rep(NA,length(p1)),
                        dists=rep("tnorm",length(p1)),
                        untrans=rep("identity",length(p1)),
                        dist.types=c("tnorm","beta","gamma","lnorm","cauchy","constant","tt")
)
  # Makes a list of prior distribution parameters.
{
  
  if (length(p2)==1)
    p2 <- rep(p2,length(p1))
  if ( length(p1)!=length(p2) )
    stop("p1 and p2 must have the same length")
  if ( length(p1)!=length(lower) )
    stop("p1 and lower must have the same length")
  if ( length(p1)!=length(upper) )
    stop("p1 and upper must have the same length")
  both.not.na <- !is.na(upper) & !is.na(lower)
  if ( any(upper[both.not.na]<=lower[both.not.na]) )
    stop("All elements of upper must be greater than lower")
  if ( length(p1)!=length(dists) )
    stop("p1 and dists must have the same length")
  if ( !all(dists %in% dist.types) )
    stop(paste("Unsupported distribution, allowable types are:",
               paste(dist.types,collapse=", ")))
  is3 <- dists %in% "tt"
  name.untrans <- length(untrans) != length(p1)
  if (name.untrans & (is.null(names(untrans)) | is.null(names(p1))))
    stop("If untrans vector is not the same length as p1 it must have p1 names")
  if (!(all(names(untrans) %in% names(p1))))
    stop("untrans vector has names not in p1 names")
  is.tt <- (names(dists)=="tt")
  if (any(is.tt) && any(is.na(p3[is.tt])))
     stop("For tt distribution must specify df in corresponding slot of p3")
  prior <- vector(mode="list",length=length(p1))
  names(prior) <- names(p1)
  for ( i in 1:length(p1) ) {
    prior[[i]] <- switch(dists[i],
                         tnorm={
                           if (is.na(lower[i])) lower[i] <- -Inf
                           if (is.na(upper[i])) upper[i] <- Inf             
                           p <- c(p1[i],p2[i],lower[i],upper[i])
                           names(p) <- c("mean","sd","lower","upper")
                           p <- as.list(p)
                           attr(p,"dist") <- "tnorm"
                           p
                         },
                         beta={
                           if (is.na(lower[i])) lower[i] <- 0
                           if (is.na(upper[i])) upper[i] <- 1             
                           p <- c(p1[i],p2[i],lower[i],upper[i])
                           names(p) <- c("shape1","shape2","lower","upper")
                           p <- as.list(p)
                           attr(p,"dist") <- "beta_lu"
                           p
                         },
                         gamma={
                           if (is.na(lower[i])) lower[i] <- 0
                           p <- c(p1[i],p2[i],lower[i])
                           names(p) <- c("shape","scale","lower")
                           p <- as.list(p)
                           attr(p,"dist") <- "gamma_l"
                           p
                         },
                         lnorm={
                           if (is.na(lower[i])) lower[i] <- 0
                           p <- c(p1[i],p2[i],lower[i])
                           names(p) <- c("meanlog","sdlog","lower")
                           p <- as.list(p)
                           attr(p,"dist") <- "lnorm_l"
                           p
                         },
                         cauchy={
                           p <- c(p1[i],p2[i])
                           names(p) <- c("location","scale")
                           p <- as.list(p)
                           attr(p,"dist") <- "cauchy"
                           p
                         },
                         tt={
                           if (is.na(lower[i])) lower[i] <- -Inf
                           if (is.na(upper[i])) upper[i] <- Inf             
                           p <- c(p1[i],p2[i],p3[i],lower[i],upper[i])
                           names(p) <- c("location","scale","df","lower","upper")
                           p <- as.list(p)
                           attr(p,"dist") <- "tt"
                           p
                         },
                         {
                           p <- p1[i]
                           names(p) <- c("constant")
                           p <- as.list(p)
                           attr(p,"dist") <- "constant"
                           p
                         }
    )
    prior[[i]]$log <- TRUE
    if (!name.untrans) attr(prior[[i]],"untrans") <- untrans[i] else
      if (is.na(untrans[names(p1)[i]])) 
        attr(prior[[i]],"untrans") <- "identity" else
        attr(prior[[i]],"untrans") <- untrans[names(p1)[i]]
  }
  prior
}





log.prior.dmc=function(p.vector,p.prior)
  # log of prior density for p.vector  
{
  prior <- numeric(length(p.vector))
  names(prior) <- names(p.vector)
  for ( i in names(p.vector) ) {
    p.prior[[i]]$x <- p.vector[i]
    prior[i] <- do.call(paste("d",attr(p.prior[[i]],"dist"),sep=""),
                        p.prior[[i]])
  }
  prior
}


rprior.dmc=function(p.prior,n=1)
  # sample from prior  
{
  prior <- matrix(nrow=n,ncol=length(p.prior))
  dimnames(prior) <- list(NULL,names(p.prior))
  for ( i in 1:length(names(p.prior)) ) {
    p.prior[[i]]$n <- n
    prior[,i] <- do.call(paste("r",attr(p.prior[[i]],"dist"),sep=""),
                         p.prior[[i]][names(p.prior[[i]])!="log"])
  }
  prior
}


log.posterior.dmc=function(p.vector,p.prior,data,min.like=1e-10)
  # Summed log posterior likelihood
{
  sum (log.likelihood(p.vector,data,min.like=min.like)) + 
    summed.log.prior (p.vector, p.prior,
      p.names=names(attr(attributes(data)$model,"p.vector")))
  
}


log.likelihood=function (p.vector, data, min.like=1e-10)
  # Get log likelihood
{
  names(p.vector)=names(attr(attributes(data)$model,"p.vector")) 
  suppressWarnings( log(likelihood.dmc(p.vector,data,min.like=min.like) ) )
}


summed.log.prior=function (p.vector, p.prior, p.names)
{
  names(p.vector)=p.names 
  suppressWarnings( sum(log.prior.dmc(p.vector,p.prior)) )
}  


assign.pp <- function(pp,p.prior)
  # Slot pp values into p.prior
{
  for (i in 1:length(p.prior)) 
    p.prior[[i]][1:2] <- c(pp[[1]][i],pp[[2]][i])
  p.prior
}



######### Sampling ----

# thin=1;samples=NULL;theta1=NULL; data=NULL; p.prior=NULL
# restart=TRUE;add=FALSE;remove=NA;start.from=NA
# start.prior=NULL;rp=.001;verbose=TRUE; replace.bad.chains=NULL
# n.chains=length(attr(attr(data,"model"),"p.vector"))*3

samples.dmc <- function(nmc,p.prior=NULL,data=NULL,thin=1,samples=NULL,theta1=NULL,
                        restart=TRUE,add=FALSE,remove=NA,start.from=NA,
                        start.prior=NULL,rp=.001,verbose=TRUE,replace.bad.chains=NULL,cut=5, 
                        n.chains=length(attr(attr(data,"model"),"p.vector"))*3)
  # Setup for sampling  
{
  if ( !is.null(data) && !is.data.frame(data) )       
    stop("For single subject data must be a data frame")
  if ( !is.null(samples) && is.null(samples$theta) )  
    stop("Use h.samples.dmc for a list of subjects")

  if ( is.null(samples) )   # setup new samples
  { 
    model <- attributes(data)$model
    
    if ( is.null(model) )
      stop("Must specify a model")
    if ( is.null(p.prior) )
      stop("Must specify a p.prior argument")

    p.names <- names(attr(model,"p.vector"))
    n.pars <- length(p.names)
    theta <- array(NA,c(n.chains,n.pars,nmc))   # sampled parameters
    dimnames(theta)[[2]] <- p.names
    if (!is.null(theta1) && !is.matrix(theta1) || 
        (!all(dim(theta1)==c(n.chains,n.pars))) )
      stop("If theta1 is supplied it must be a n.chains x n.pars matrix")
    
    summed_log_prior   <- array(-Inf,c(nmc,n.chains))      
    log_likelihoods    <- array(-Inf,c(nmc,n.chains))
    
    # Generate random start points from prior
    if (verbose) cat("Generating start points for each chain: ")
      
    for( i in 1:n.chains ) 
    {
      if (verbose) cat(".")
      j <- 1

      while ( any(is.infinite(summed_log_prior[1,i])) ||
              any(is.infinite(log_likelihoods[1,i])) )  # Make sure parameters are valid.
      {
        if ( is.null(theta1) ) {           # Generate parameters from prior
          if ( is.null(start.prior) ) {   
            theta[i,,1] <- rprior.dmc(p.prior)[,p.names]  # sample from prior
          } else {
            theta[i,,1] <- rprior.dmc(start.prior)[,p.names] # sample start.prior 
          }
        } else { # Copy parameters from theta 1 
          theta[i,p.names,1] <- theta1[i,p.names]
        }
        
        log_likelihoods[1,i]   <- sum(log.likelihood (theta[i,,1],data))
        summed_log_prior[1,i]  <- summed.log.prior (theta[i,,1],p.prior,
          p.names=names(attr(attributes(data)$model,"p.vector")))

        if ( !is.null(theta1) && (any(any(is.infinite(summed_log_prior[1,i])) ||
              any(is.infinite(log_likelihoods[1,i])))) ) {
          print(theta[i,,1])
          stop("Bad theta1")
        }

        j <- j + 1
        if (j > 1e4) 
          stop("Oops, could not sample valid data-level parameter from prior after 10,000 tries")
      }
    }

    if (verbose) cat("\n")
    
    samples <- list(theta=theta,
                    summed_log_prior=summed_log_prior,
                    log_likelihoods=log_likelihoods,
                    data=data,p.prior=p.prior,start=1,
                    n.pars=n.pars,p.names=p.names,rp=rp,
                    nmc=nmc,thin=thin,n.chains=n.chains)
    
  } else { # restart from previous sampling

    # Clean up remove
    if (!all(is.na(remove))) {
      remove <- remove[(remove>0) & (remove<=dim(samples$theta)[3])]
      if (length(remove)==0) remove <- NA
    }
    # remove some samples 
    if ( !all(is.na(remove)) ) {
      if ( !all(remove %in% 1:dim(samples$theta)[3]) )
        stop(paste("Remove range must be a sequence in 1:",
                   dim(samples$theta)[3],sep=""))
      if (length(remove)==length(1:dim(samples$theta)[3]))
        stop("Cannot remove all samples")
      
      samples$theta <- samples$theta[,,-remove]
      samples$nmc <- dim(samples$theta)[3]
      
      samples$summed_log_prior <- samples$summed_log_prior[-remove,]
      samples$log_likelihoods  <- samples$log_likelihoods[-remove,]
      
    }
    
    if ( add ) { # add on
      theta <- array(NA,c(samples$n.chains,samples$n.pars,samples$nmc+nmc))   
      dimnames(theta)[[2]] <- samples$p.names
      theta[,,1:samples$nmc] <- samples$theta
      samples$theta <- theta
      
      summed_log_prior <- array(-Inf,c(samples$nmc+nmc,samples$n.chains))         
      summed_log_prior[1:samples$nmc,] <- samples$summed_log_prior
      samples$summed_log_prior <- summed_log_prior
      
      log_likelihoods  <- array(-Inf,c(samples$nmc+nmc,samples$n.chains))
      log_likelihoods[1:samples$nmc,] <- samples$log_likelihoods
      samples$log_likelihoods  <- log_likelihoods 
      
      samples$start <- samples$nmc
      samples$nmc <- samples$nmc+nmc
    }
    
    if ( !add & restart ) { # start afresh
      old.nmc <- dim(samples$summed_log_prior)[1]
      if (!is.na(start.from)) {
        if (start.from > old.nmc) 
          stop(paste("start.from must not be >",old.nmc,"\n"))
        old.nmc <- start.from
      }
      theta1 <- samples$theta[,,old.nmc]
      summed_log_prior1 <- samples$summed_log_prior[old.nmc,]
      log_likelihoods1  <- samples$log_likelihoods[old.nmc,]
     
      if (is.logical(replace.bad.chains)) {
        bad <- pick.stuck.pars.dmc(samples,cut=cut,start=1,end=dim(samples$theta)[3])
        replace.bad.chains <- unique(unlist(bad[unlist(lapply(bad,function(x){length(x)>0}))]))
        if (length(replace.bad.chains)==0) replace.bad.chains <- NULL else { 
          cat("Replacing bad chains\n")
          print(replace.bad.chains)
        }
      }
      if (!is.null(replace.bad.chains)) {
        if (!all(replace.bad.chains %in% 1:samples$n.chains))
          stop(paste("replace.bad.chains must be in the range 1 to",samples$n.chains))
        good.chains <- c(1:samples$n.chains)[-replace.bad.chains]
        if (length(good.chains)<2)
          stop("Must be at least two good chains")
        good.chains <- sample(good.chains,length(replace.bad.chains),replace=TRUE)
        theta1[replace.bad.chains,] <- theta1[good.chains,] 
        summed_log_prior1[replace.bad.chains] <- summed_log_prior1[good.chains]
        log_likelihoods1[replace.bad.chains] <- log_likelihoods1[good.chains]
      }
  
      samples$theta=array(NA,c(samples$n.chains,samples$n.pars,nmc))
      dimnames(samples$theta)[[2]] <- samples$p.names
      samples$theta[,,1] <- theta1
      
      samples$summed_log_prior <- array(-Inf,c(nmc,samples$n.chains))         
      samples$summed_log_prior[1,]  <- summed_log_prior1

      samples$log_likelihoods  <- array(-Inf,c(nmc,samples$n.chains))
      samples$log_likelihoods[1,]  <- log_likelihoods1

      samples$nmc <- nmc
      samples$thin <- thin
      samples$start <- 1

    }
  }
  samples
}


crossover <- function(k,pars,use.theta,use.logprior,use.loglike,p.prior,data,
                      rp,gamma.mult=2.38,force=FALSE)
  # DEMCMC crossover update of one chain, data level
{
  # step size
  if (is.na(gamma.mult)) gamma <- runif(1,0.5,1) else
    gamma <-  gamma.mult/sqrt(2*length(pars))
  # pick two other chains
  index <- sample(c(1:dim(use.theta)[1])[-k],2,replace=F)

  # DE step
  theta <- use.theta[k,]
  names(theta)=names(attr(attributes(data)$model,"p.vector")) 

  theta[pars] <- 
    use.theta[k,pars] + 
    gamma*(use.theta[index[1],pars]-use.theta[index[2],pars]) + 
    runif(1,-rp,rp)  
  
  # Old post
  if (force) {
    use.logprior[k] <- summed.log.prior (p.vector=use.theta[k,],p.prior=p.prior,
      p.names=names(attr(attributes(data)$model,"p.vector")))
    use.loglike[k] <- sum(log.likelihood (p.vector=use.theta[k,],data=data))
  }
  summed.use.post <-  use.loglike[k] + use.logprior[k]
  
  # new post
  log.prior <- summed.log.prior (p.vector=theta,p.prior=p.prior,
    p.names=names(attr(attributes(data)$model,"p.vector")))
  loglike  <- sum(log.likelihood (p.vector=theta,data=data))
  post <- log.prior + loglike
  if ( is.na(post) ) post <- -Inf
  
  # Metropolis step
  if ( runif(1) < exp(post-summed.use.post) ) {
    use.theta[k,]   <- theta
    use.logprior[k] <- log.prior
    use.loglike[k]  <- loglike
  }
 
  c(use.logprior[k], use.loglike[k], use.theta[k,])
}


# x is data and theta is n.chain x n.par matrix
migrate.yishin <- function(theta, x) {
  out    <- theta         ## clone a theta copy, so I can change theta internally
  nchain <- nrow(theta)
  npar   <- ncol(theta)
  subchains <- 1 + ggdmc::GetSubchains(nchain)[,1]
  nsubchain <- length(subchains)
  rj <- numeric(nchain)   ## record reject rates; default is 0, accepting.
  theta_star <- numeric(npar)

  for (i in 1:nsubchain) {
    if (i == nsubchain) {
      next_chain <- subchains[1]
    } else {
      next_chain <- subchains[i+1]
    }

    k <- subchains[i]
    theta_cur  <- theta[next_chain, ]
    for(j in 1:npar) { # normal perturbation
      theta_star[j] <- rnorm(1, theta[k, j]) 
    } 
    ## for(j in 1:npar) { theta_star[j] <- theta_star[j] + runif(1, -0.01, 0.01)}

    ## Proposal posterior
    tmp_logprior <- dnorm(theta_star[1], 0, 8, log = TRUE) +
      ggdmc::dtnorm(theta_star[2], 0, 8, 0, Inf, log = TRUE)[,1]
    tmp_loglike  <- sumloglike_(theta_star, x)
    tmp_logpos   <- tmp_logprior + tmp_loglike
    if (is.na(tmp_logpos)) tmp_logpos <- -Inf

    ## Current posterior
    cur_logprior <- dnorm(theta_cur[1], 0, 8, log = TRUE)  +
      ggdmc::dtnorm(theta_cur[2], 0, 8, 0, Inf, log = TRUE)[,1]
    cur_loglike <- sumloglike_(theta_cur, x)
    cur_logpos <- cur_logprior + cur_loglike
    if (is.na(cur_logpos)) cur_logpos <- -Inf

    ratio <- (tmp_logpos - cur_logpos)
    if (is.nan(ratio)) ratio <- -Inf

    ## Decision step
    if (exp(ratio) > runif(1)) {
      ## accept the proposal
      out[next_chain, ] <- theta_star
      rj[next_chain] <- FALSE
    } else {
      ## continue with theta_cur; no reject
      out[next_chain, ]  <- theta_cur
      rj[next_chain] <- TRUE
    }
  }

  return(list(out, rj))
}


migrate <- function(use.theta,use.logprior,use.loglike,
                    p.prior,data,rp)
  # DEMCMC migrate set, all chains, data level
{
  n.chains <- dim(use.theta)[1]
  pars  <- dim(use.theta)[2]
  n.groups <- sample(c(1:n.chains),1)  # determine how many groups to work with
  groups <- sort(sample(c(1:n.chains),n.groups,replace=F))	# get groups 
  
  thetaset <- matrix(NA,n.groups,pars)									   # initialize
  currentset.logprior <- propset.logprior <- numeric(n.groups)
  propw.logprior      <- currw.logprior   <- numeric(n.groups)
  currentset.loglike  <- propset.loglike  <- numeric(n.groups)
  propw.loglike       <- currw.loglike    <- numeric(n.groups)

  for (i in 1:n.groups) {
    # create a set of these particles to swap
    thetaset[i,] <- use.theta[groups[i],] + runif(1,-rp,rp)	

    currentset.logprior[i] <- use.logprior[groups[i]]
    currentset.loglike[i]  <- use.loglike[groups[i]]

    propset.logprior[i] <- summed.log.prior (p.vector=thetaset[i,], 
      p.prior=p.prior,p.names=names(attr(attributes(data)$model,"p.vector")))
    propset.loglike[i]  <- sum(log.likelihood (p.vector=thetaset[i,], 
                                                 data=data))
    propw.logprior[i]      <- propset.logprior[i]
    propw.loglike[i]       <- propset.loglike[i]
    
    currw.logprior[i]      <- currentset.logprior[i]
    currw.loglike[i]       <- currentset.loglike[i]
    
  }

  mh <- exp( (propset.loglike[n.groups] + propset.logprior[n.groups]) - 
                      (currentset.loglike[1] + currentset.logprior[1]))
  if ( !is.na(mh) && (runif(1) < mh) ) {
    use.theta[groups[1],]   <- thetaset[n.groups,]	# swap the 1st with last 
    use.logprior[groups[1]] <- propset.logprior[n.groups]
    use.loglike[groups[1]]  <- propset.loglike[n.groups]
  }
  
  if ( n.groups!=1 ) {										# make sure we are not done yet
    for(i in 1:(n.groups-1)) {		
      mh <- exp((propset.loglike[i] + propset.logprior[i]) - 
                        (currentset.loglike[i+1] + currentset.logprior[i+1]))
      if( !is.na(mh) && (runif(1) < mh) ) {
        use.theta[groups[i+1],] <- thetaset[i,]	
        use.logprior[groups[i+1]] <- propset.logprior[i]
        use.loglike[groups[i+1]] <- propset.loglike[i]
      }
    }
  }
  
  cbind (use.logprior, use.loglike, use.theta)
}

run.dmc <- function(samples,report=10,cores=1,p.migrate=0,
  gamma.mult=2.38,farjump=NA,force=FALSE,verbose=TRUE)
  # Run sampling. 
  # Force is a boolean for each interation forcing acceptance (not used in lessions).
{
  
  os <- get.os()
  if (os == "windows" & cores > 1) {
    require(snowfall,quietly=TRUE)
    require(rlecuyer,quietly=TRUE)
    sfInit(parallel=TRUE, cpus=cores, type="SOCK")
    sfClusterSetupRNG()
    sfLibrary(msm)
    sfLibrary(rtdists)
    sfLibrary(pracma)
    sfLibrary(statmod)
    sfExportAll() 
  } else if (cores > 1) {
    require(parallel, quietly=TRUE)
  } 
  
  if (!is.na(farjump) && !is.numeric(farjump) && (farjump<1))
    stop("farjump must be a number greater than or equal to 1") else 
    farjump <- round (farjump)
    
  n.trials <- dim(samples$data)[1]
  if ( is.null(samples$thin) || is.na(samples$thin) ) samples$thin <- 1
  nsamp <- 1+(samples$nmc-samples$start)*samples$thin

  if (is.numeric(force)) force <- c(rep(FALSE,force-1),TRUE)
  force <- c(TRUE,rep(force,length.out=nsamp-1))
  if ( !(all(is.logical(force))) )
      stop(paste("force argument must be a logical scalar or vector"))

  use.theta <- samples$theta[,,samples$start]
  use.logprior <- samples$summed_log_prior[samples$start,]
  use.loglike <- samples$log_likelihoods[samples$start,]

  store_i <- samples$start
  
  for (i in 2:nsamp) {

    gamma.mult.i <- gamma.mult
    if (!is.na(farjump) && (i %% farjump == 0)) gamma.multi.i <- .98 
    
    if ( runif(1)<p.migrate ) {         # Do migration
      temp <- migrate(use.theta          = use.theta,
                      use.logprior       = use.logprior,
                      use.loglike        = use.loglike,
                      p.prior            = samples$p.prior,
                      data               = samples$data,
                      rp                 = samples$rp) 
    } else {                            # Do crossover
      if ( cores==1 ) {                 # (Non-snowfall)
          temp=t(sapply(1:samples$n.chains,crossover,
            pars = 1:samples$n.pars,
            use.theta          = use.theta,
            use.logprior       = use.logprior,
            use.loglike        = use.loglike,
            p.prior            = samples$p.prior,
            data               = samples$data,
            rp                 = samples$rp,
            force              = force[i],
            gamma.mult         = gamma.mult.i))
        } else if (cores > 1 & os == "windows") { # via Snowfall
          temp=t(sfLapply(1:samples$n.chains,crossover,
            pars = 1:samples$n.pars,
            use.theta          = use.theta,
            use.logprior       = use.logprior,
            use.loglike        = use.loglike,
            p.prior            = samples$p.prior,
            data               = samples$data,
            rp                 = samples$rp,
            force              = force[i],
            gamma.mult         = gamma.mult.i))
          
          temp <- t(array(unlist(temp),dim=c(samples$n.pars+2, samples$n.chains)))
        } else { # Via parallel
          temp <- t(mclapply(1:samples$n.chains, crossover,
            pars = 1:samples$n.pars,
            use.theta          = use.theta,
            use.logprior       = use.logprior,
            use.loglike        = use.loglike,
            p.prior            = samples$p.prior,
            data               = samples$data,
            rp                 = samples$rp,
            force              = force[i],
            gamma.mult         = gamma.mult, mc.cores=cores))
          
          temp <- t(array(unlist(temp),dim=c(samples$n.pars+2,samples$n.chains)))
        }
      }
      
    use.logprior  <- temp[,1]                                  
    use.loglike   <- temp[,2]
    use.theta     <- temp[,3:(samples$n.pars+2)]

    if (i %% samples$thin == 0) { # store samples
      store_i <- store_i + 1
      if (verbose & (store_i %% report == 0)) cat(store_i," ")
      samples$summed_log_prior[store_i,]  <- use.logprior                                  
      samples$log_likelihoods[store_i,]   <- use.loglike
      samples$theta[,,store_i]            <- use.theta
    }
    
  }
  if (verbose) cat("\n")
  if (cores>1 & os == "windows") { sfStop() }
  samples
}

### Posterior predictives ----


get.dqp <- function(sim,facs,probs,n.post=NA,ns=NA,bw="nrd0") {

    quantile.names <- function(x,probs=seq(0, 1, 0.25),na.rm=FALSE,type=7, ...) {
      out <- quantile(x,probs=probs,na.rm=na.rm,type=type,names=FALSE,...)
      names(out) <- probs*100
      if ( all(is.na(out)) || length(x)==1) NULL else out
    }
    
    qs <- tapply(sim$RT,sim[,c(facs,"R")],quantile.names,probs=probs,na.rm=TRUE)

    # get probabilities 
    n <- tapply(sim$RT,sim[,c(facs,"R")],length)
    n[is.na(n)] <- 0 # In case some cells are empty
    nok <- tapply(sim$RT,sim[,c(facs,"R")],function(x){sum(!is.na(x))})
    nok[is.na(nok)] <- 0 # In case some cells are empty
    if ( is.null(facs) ) np <- sum(n) else
      np <- rep(apply(n,1:length(facs),sum),times=length(levels(sim$R)))
    p <- nok/np
    # p.all <- n/np

    # For a simulation get probability replicates
    if ( !is.na(n.post) && (n.post>1) ) { 
      repfac <- rep(1:n.post,each=sum(ns))
      # ns <- tapply(sim$RT,cbind(sim[,c(facs,"R"),drop=FALSE],rep=repfac),length)
      # ns[is.na(ns)] <- 0 # In case some cells are empty
      ps <- tapply(sim$RT,cbind(sim[,c(facs,"R"),drop=FALSE],rep=repfac),
        function(x){sum(!is.na(x))})
      ps[is.na(ps)] <- 0 # In case some cells are empty
      ps <- n.post*ps/np  
      # ps.all <- n.post*ns/np  
    } else ps=NULL

    # cell names
    cell.names <- dimnames(qs)[[1]]
    if (!is.null(facs)) {
      n.cell <- length(facs)+1
      for ( i in 2:n.cell )
        cell.names <- outer(cell.names,dimnames(qs)[[i]],"paste",sep=".")
    }
    cell.names <- as.vector(cell.names)
    # Get density and make defective
    dens <- tapply(sim$RT,sim[,c(facs,"R"),drop=FALSE],function(x){
      if (all(is.na(x))) NULL else {
        x <- x[x>=quantile(x,.01,na.rm=TRUE) & x<=quantile(x,.99,na.rm=TRUE)]
        if (length(x[!is.na(x)])<2) NULL else density(x[!is.na(x)],bw=bw)
      }
    })
    for (i in 1:length(p)) if ( is.finite(p[i]) && !(p[i]==0) ) 
    {
      if (!is.null(qs[i][[1]])) {
        names(qs[i][[1]]) <- as.numeric(names(qs[i][[1]]))*p[i]/100
        attr(qs[i][[1]],"cell.name") <- cell.names[i]
      }
      if (!is.null(dens[i][[1]]) ) {
        dens[i][[1]]$y <- dens[i][[1]]$y*p[i]
        attr(dens[i][[1]],"cell.name") <- cell.names[i]
      }
    }
    dnd <- dimnames(dens)
    dnq <- dimnames(qs)
    dens <- apply(dens,1:length(facs),function(x){x})
    qs <- apply(qs,1:length(facs),function(x){x})
    if ( is.null(dim(dens)) ) {
      dens <- array(dens,dim=c(length(dens)))
      dimnames(dens) <- dnd[-length(dnd)] 
      qs <- array(qs,dim=c(length(qs)))
      dimnames(qs) <- dnq[-length(dnq)] 
    }
    list(pdf=dens,cdf=qs,n=n,nok=nok,p=p,ps=ps) #,n=n,p.all,ps.all,ns=ns)
}



# n.post=100;probs=c(1:99)/100;random=TRUE
# bw="nrd0";report=10;save.simulation=FALSE; factors=NA
# save.simulation.as.attribute=FALSE;ignore.R2=TRUE
# gglist=TRUE; probs.gglist=c(0.1, 0.5, 0.9); CI.gglist=c(0.025, 0.975)
#  censor=c(NA,NA)

# samples=hsamples[[1]]; random=TRUE; n.post=2; gglist=TRUE; report=1

post.predict.dmc <-function(samples,n.post=100,probs=c(1:99)/100,random=TRUE,
  bw="nrd0",report=10,save.simulation=FALSE,factors=NA,
  save.simulation.as.attribute=FALSE,ignore.R2=FALSE,censor=c(NA,NA),
  gglist=TRUE, probs.gglist=c(0.1, 0.5, 0.9),CI.gglist=c(0.025, 0.975))
  # make list of posterior preditive density, quantiles and response p(robability)
  # NB: quantiles only calcualted for 2 or more RTs
{ 
  
  model <- attributes(samples$data)$model
  facs <- names(attr(model,"factors"))
  cvs <- samples$data[,attr(model,"cvs")]
  attr(cvs,"row.facs") <- apply(apply(
    samples$data[,facs,drop=FALSE],2,as.character),1,paste,collapse=".")
  if ( ignore.R2 & any(names(samples$data)=="R2") )
    samples$data <- samples$data[,names(samples$data)[names(samples$data)!="R2"]]
  if (!is.null(factors) ) {
    if (any(is.na(factors))) factors <- facs
    if (!all(factors %in% facs)) 
      stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
  }
  resp <- names(attr(model,"responses"))
  ns <- table(samples$data[,facs],dnn=facs)
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(3,1,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  if (is.na(n.post)) use <- c(1:dim(thetas)[1]) else {
    if (random) use <- sample(c(1:dim(thetas)[1]),n.post,replace=F) else
      use <- round(seq(1,dim(thetas)[1],length.out=n.post))
  }
  n.post <- length(use)
  posts <- thetas[use,]
  n.rep <- sum(ns)
  sim <- data.frame(matrix(nrow=n.post*n.rep,ncol=dim(samples$data)[2]))
  names(sim) <- names(samples$data)
  # Tweaks for Stop Signal
  if ( !any(names(samples$data)=="SSD") ) {
    SSD <- rep(Inf,sum(ns)) 
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data) %in% c("RT",names(cvs),"R2")]
  } else {
    # Assumes last two are SSD and RT! FIX ME. EG WONT WORK IF THERE ARE CVS
    if ( is.null(facs) ) SSD <- samples$data$SSD else
      SSD <- unlist(tapply(samples$data$SSD,samples$data[,facs],identity)) 
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data) %in% c("RT","SSD")]
    # leave.out <- -c((dim(samples$data)[2]-1):dim(samples$data)[2])
  }
  cat("\n")
  cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  for (i in names(samples$data)[leave.out])
    sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
  for (i in 1:n.post) {
    tmp <- simulate.dmc(posts[i,],model,n=ns,SSD=SSD,cvs=cvs)
    if (ignore.R2) tmp <- tmp[,names(tmp)[names(tmp)!="R2"]]
    sim[(1+(i-1)*n.rep):(i*n.rep),names(tmp)] <- tmp
    if ( (i %% report) == 0) cat(".")
  }
  if ( any(names(sim)=="R2") ) { # MTR model
    levs <- outer(levels(samples$data$R),sort(unique(samples$data$R2)),"paste",sep="")
    if (attributes(model)$type=="normDK") levs[,2] <- rep("DK",dim(levs)[1])
    levs <- sort(unique(as.vector(levs)))  
    sim$R <- paste(as.character(sim$R),as.character(sim$R2),sep="") 
    if (attributes(model)$type=="normDK") sim$R[sim$R2=="2"] <- "DK"
    sim$R <- factor(sim$R,levels=levs)
    samples$data$R <- paste(as.character(samples$data$R),as.character(samples$data$R2),sep="") 
    if (attributes(model)$type=="normDK") 
      samples$data$R[samples$data$R2=="2"] <- "DK"
    samples$data$R <- factor(samples$data$R,levels=levs)
  }
  reps <- rep(1:n.post,each=dim(samples$data)[1])
  
  if (!is.na(censor[1])) fast <- sim[,"RT"] < censor[1] else fast <- rep(FALSE,dim(sim)[1])
  if (!is.na(censor[2])) slow <- sim[,"RT"] > censor[2] else slow <- rep(FALSE,dim(sim)[1])
  ok <- !fast & !slow 
  sim <- sim[ok,]
  reps <- reps[ok]
  
  if ( save.simulation ) {
    sim <- cbind(reps,sim)
    attr(sim,"data") <- samples$data
    sim
  } else {
    sim.dqp <- get.dqp(sim,facs=factors,probs,n.post,ns=ns,bw=bw)
    dat.dqp <- get.dqp(sim=samples$data,factors,probs,bw=bw)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    out <- c(sim.dqp,dat.dqp)
    dpqs <- vector(mode="list",length=length(n.post))
    for (i in 1:n.post) {
      simi <- sim[reps==i,]
      dpqs[[i]] <- get.dqp(simi,factors,probs,1)
    }
    attr(out,"dpqs") <- dpqs
    if (save.simulation.as.attribute) 
      attr(out,"sim") <- cbind(reps,sim)
    if (gglist) attr(out, "gglist") <- 
      get.fitgglist.dmc(sim=cbind(reps,sim),data=samples$data,factors=factors, noR=FALSE, 
        quantiles.to.get= probs.gglist, CI = CI.gglist)
    out
  }
}           

### Automatic sampling ----

run.unstuck.dmc <- function(samples,nmc=NA,report=10,cores=1,
                            cut=10,nbad=0,max.try=100,p.migrate=0,
                            gamma.mult=2.38,verbose=FALSE,
                            end.no.migrate=FALSE)
  # Repeats sampling until <= nbad stuck chains as defined by cut or max.try
  # If samples has no content fills it in then repeatedly gets new sets of nmc
  # samples (nmc can be specified or taken from samples). If end.no.migrate
  # runs one final time with migration off. 
{
  if (is.null(samples$theta))
    stop("For multiple subjects use h.run.unstuck.dmc")
  n.chain <- dim(samples$theta)[1]
  if (any(is.na(samples$theta[,,2])))
    samples <- run.dmc(samples=samples,report=report,cores=cores,
                       gamma.mult=gamma.mult,p.migrate=p.migrate)
  if ( is.na(nmc) ) nmc <- samples$nmc
  try.num <- 1
  repeat {
    cat(paste("\nTry",try.num,"\n"))
    if ( length(pick.stuck.dmc(samples,verbose=verbose)) <= nbad ) break
    samples <- run.dmc(samples.dmc(nmc=nmc,samples=samples), 
                       report=report,cores=cores,p.migrate=p.migrate,
                      gamma.mult=gamma.mult)
    if (try.num >= max.try) break
    try.num <- try.num + 1
  }
  if (end.no.migrate) samples <- run.dmc(samples.dmc(nmc=nmc,samples=samples), 
                       report=report,cores=cores,gamma.mult=gamma.mult)
  samples
}


run.converge.dmc <- function(samples,nmc,report=10,cores=1,gamma.mult=2.38,
  cut=1.1,max.try=100,minN=NA,meanN=NA,
  transform=TRUE,autoburnin=FALSE,split=TRUE,verbose=FALSE)
  # Adds samples repeatedly, throws away intial samples if gelman.diag better
  # and repeats until gelman.diag Multivariate psrf < cut and once that is
  # fulfilled will continue if necessary to get either min or mean effectiveSize
{
  if (!is.na(minN) & !is.na(meanN)) {
    warning("Both minN and meanN specified, using minN")
    meanN <- NA 
  }
  if (!is.na(minN)) nfun <- "min"
  if (!is.na(meanN)) {
    nfun <- "mean"
    minN <- meanN
  }
  if ( is.null(samples$theta) )
    stop("For multiple subjects use h.run.converge.dmc")
  if ( any(is.na(samples$theta[,,2])) ) {
    samples <- run.dmc(samples=samples, 
                       report=report,cores=cores,gamma.mult=gamma.mult)
    gd <- gelman.diag.mpsrf(theta.as.mcmc.list(samples,split=split),
                      autoburnin=autoburnin,transform=transform)
    if (verbose) cat(paste("MPSRF: ",gd,"\n"))
  } else gd <- Inf
  if ( !is.na(minN) & (gd <= cut) ) 
    okN <- do.call(nfun,list(effectiveSize.dmc(samples))) > minN else
    okN <- TRUE
  if ( (gd > cut) | !okN ) 
    {  # Do more sampling
    if ( is.na(nmc) ) nmc <- samples$nmc
    try.num <- 0
    effectiveN <- NA
    repeat {
      samples <- run.dmc(samples.dmc(samples=samples,add=TRUE,nmc=nmc),
        report=report,cores=cores,gamma.mult=gamma.mult)
      gd <- gelman.diag.mpsrf(theta.as.mcmc.list(samples,split=split),
                      autoburnin=autoburnin,transform=transform)
      if ( try.num>0 ) {
        shorter <- samples.dmc(samples=samples,remove=1:nmc,nmc=0,add=TRUE) 
        gd.short <- gelman.diag.mpsrf(theta.as.mcmc.list(shorter,split=split),
                        autoburnin=autoburnin,transform=transform)
        if (gd.short < gd) {
          samples <- shorter
          gd <- gd.short
          if (verbose) cat(paste("Discarding initial",nmc,"samples.\n"))
        }
      }
      if ( try.num >= max.try ) break
      try.num <- try.num + 1
      if ( gd <= cut ) {
        if ( is.na(minN) ) break else {
          effectiveN <- do.call(nfun,list(effectiveSize.dmc(samples)))
          if (effectiveN > minN) break
        } 
      }
      if (verbose) print(paste("N =",dim(samples$theta)[3],
                               "Effective N =",effectiveN,
                               "Multivariate psrf achieved =",format(gd,digits=3)))
    }
  }
  if (verbose) {
    print(paste("Final multivariate psrf =",gd))
    cat("Effective sample size\n")
    print(effectiveSize.dmc(samples))
  }
  samples
}


# cores=1;p.migrate=.05;report=10;verbose=FALSE;split=TRUE; force=FALSE
# gamma.mult=2.38;minN=NA;meanN=NA;cut.unstuck=10;cut.converge=1.1;
# cut.flat.location=1/2;cut.flat.scale=1/2;max.try=100; n.add=NA

RUN.dmc <- function(samples,cores=1,report=10,p.migrate=.05,max.try=20,
  cut.unstuck=10,
  cut.flat.location=1/2,cut.flat.scale=1/2,
  cut.converge=1.1,split=TRUE,
  minN=NA,meanN=NA,use.effectiveSize = TRUE,
  n.add=NA,force=FALSE,
  verbose=FALSE,gamma.mult=2.38)
  
{

  StuckTests <- function(samples,verbose=FALSE,cut) {
    if (verbose) cat("Stuck chains check\n") 
    stucks <- pick.stuck.dmc(samples,cut=cut,verbose=verbose)
    fail <- length( stucks != 0)
    if (verbose) {
      if (!fail) cat(": OK\n") else
      cat(paste(":",length(stucks),"\n")) 
    }
    fail
  }

  FlatTests <- function(samples,p1=1/3,p2=1/3,cut.location=0.25,cut.scale=Inf,
    verbose=FALSE) {

   gmcmc <- function(samples) mcmc(matrix(aperm(samples$theta,c(1,3,2)),
    ncol=dim(samples$theta)[2],dimnames=list(NULL,dimnames(samples$theta)[[2]])))

    mat <- gmcmc(samples)
    xlen <- round(dim(mat)[1] * p1)
    ylen <- round(dim(mat)[1] * p2)
    # change in mean relative to robst SD
    m.zs <- apply(mat,2,function(x){
      m1 <- median(x[1:xlen])
      m2 <- median(x[(length(x)-ylen):length(x)])
      abs(m1-m2)/IQR(x)
    })
    names(m.zs) <- paste("m",names(m.zs),sep="_")
    fail <- any(m.zs>cut.location)
    if (!fail) out <- "" else
      out <- paste(names(m.zs)[m.zs==max(m.zs)],"=",round(max(m.zs),2))
    if ( is.finite(cut.scale) ) {
      # Change ini IQR reltiave to overall IQR 
      s.zs <- apply(mat,2,function(x){
        m1 <- IQR(x[1:xlen])
        m2 <- IQR(x[(length(x)-ylen):length(x)])
        abs(m1-m2)/IQR(x)
      })
      names(s.zs) <- paste("s",names(s.zs),sep="_")
      if (out != "") out <- paste(out,", ",sep="")
      if (any(s.zs>cut.scale)) out <- 
        paste(out,names(s.zs)[s.zs==max(s.zs)],"=",round(max(s.zs),2))
      fail <- fail | any(s.zs>cut.scale) 
    } 
    if (verbose) {
      cat("Flat check\n")
      print(round(m.zs,2))
      if ( is.finite(cut.scale) ) 
        print(round(s.zs,2)) else 
          if (!fail) cat(": OK\n") else 
            cat(paste(":",out,"\n"))
    }
    fail
  }

  gelman.diag.robust <- function(samples,split) 
  {
    gd <- try(gelman.diag(theta.as.mcmc.list(samples,start=1,end=samples$nmc,split=split),
          autoburnin=FALSE,transform=TRUE),silent=TRUE)
    if (class(gd)=="try-error") list(mpsrf=Inf,psrf=matrix(Inf)) else gd
  }


  MixTests <- function(samples,verbose=FALSE,cut,split=TRUE) {
    tmp <- gelman.diag.robust(samples,split=split)
    gds <- c(tmp$mpsrf,tmp$psrf[,1])
    fail <- max(gds) > cut
    if (verbose) {
      cat("Mixing check\n")
      print(round(gds,2))
      if (!fail) cat(": OK\n") else {
        nam <- names(gds)[gds==max(gds)]
        cat(paste(":",nam,"=",round(max(gds),2),"\n"))
      }
    }
    fail
  }
  
  LengthTests <- function(samples,minN,nfun,verbose=FALSE) {
    n <- do.call(nfun,list(get.size(samples)))
    fail <- n < minN
    if (verbose) { 
      cat("Length check")
      if (!fail) cat(": OK\n") else cat(paste(":",n,"\n"))  
    }
    fail
  }

  if (!verbose) report <- 1e8
  
  if (use.effectiveSize) get.size <- effectiveSize.dmc else
    get.size <- function(x){prod(dim(x$theta)[-2])}
  
  if ( !is.na(minN) & !is.na(meanN) ) {
    warning("Both minN and meanN specified, using minN")
    meanN <- NA 
  }
  if ( !is.na(minN) ) nfun <- "min"
  if ( !is.na(meanN) ) {
    nfun <- "mean"
    minN <- meanN
  }

  if ( is.na(n.add) ) n <-  ceiling(samples$nmc/3) else {
    if ( n.add >= floor(samples$nmc/2) )
      stop(paste("n.flat.test to large, must be less than",floor(samples$nmc/2)))
      n <- n.add
  }
  
  if ( any(!is.finite(samples$log_likelihoods[samples$nmc,])) ) { # New samples
    do.migrate=1
    if (verbose) cat("\nGetting initial set of samples\n")
    samples <- run.dmc(samples,report=report,verbose=verbose,
      cores=cores,p.migrate=p.migrate,gamma.mult=gamma.mult)
  } else do.migrate=0 # Samples already good, just want more.

  if ( !force && !is.null(attr(samples,"auto")) && 
       !is.na(attr(samples,"auto")) && attr(samples,"auto")!="GRID FAIL" )
    return(samples) # If already sucessfully run then dont run again unless force=TRUE
  n.try <- 0
  repeat {
    
    # DO TESTS
    
    # Stuck check
    if ( StuckTests(samples,cut=cut.unstuck,verbose=verbose) ) 
    { # Start again 
      do.migrate <- 1
      get.new <- test.again <- TRUE
    } else get.new <- test.again <- FALSE
    
    # Stationarity check
    if ( !get.new && any(FlatTests(samples,verbose=verbose,
      cut.location=cut.flat.location,cut.scale=cut.flat.scale)) ) 
    { # remove first 1/3, add new 1/3
      if (verbose) cat("Removing initial 1/3 of samples\n")
      nshift <- ceiling(samples$nmc/3)
      samples <- samples.dmc(samples=samples,remove=1:nshift,nmc=0,add=TRUE)  
      samples <- samples.dmc(samples=samples,add=TRUE,nmc=n)
      test.again <- TRUE  
      get.new <- FALSE
    } 
    
    # Mixing check
    if ( (!get.new & !test.again) && 
        MixTests(samples,verbose=verbose, cut=cut.converge, split=split) ) { 
    { # Make longer
        samples <- samples.dmc(samples=samples,nmc=n,add=TRUE)
        test.again <- TRUE
        get.new <- FALSE
      }
    }
    
    # Length check
    if ( (!get.new & !test.again) && ( !is.na(minN) && 
        LengthTests(samples,minN=minN,nfun=nfun,verbose=verbose)) ) { # Not long enough
        samples <- samples.dmc(samples=samples,nmc=n,add=TRUE)
        test.again <- TRUE
        get.new <- FALSE
    }
    
    # Sucess?
    if ( !get.new & !test.again ) {
      outcome <- "SUCESS"
      if (verbose) cat(paste("\nOUTCOME:",outcome,"AFTER TRY",n.try,"\n"))
      break
    }
    # GET MORE SAMPLES
    if ( get.new ) { # Stuck chains, start again
      if (verbose) cat("Getting new set of samples\n")
      samples <- run.dmc(samples.dmc(samples=samples,nmc=3*n,thin=samples$thin),
        report=report,cores=cores,p.migrate=p.migrate,gamma.mult=gamma.mult,verbose=verbose)
    } else {         # Test failed, update
      if (verbose) cat("Adding extra samples\n")
      samples <- run.dmc(samples,report=report,verbose=verbose,
        cores=cores,p.migrate=p.migrate*do.migrate,gamma.mult=gamma.mult)
    }
    
    # Give up?
    n.try <- n.try + 1
    if ( (n.try > max.try) ) {
      outcome <- "FAIL"
      if (verbose) cat(paste("\nOUTCOME:",outcome,"AFTER TRY",n.try,"\n"))
      break 
    } else if (verbose) cat(paste("COMPLETED TRY",n.try,"\n\n"))
  }
  if (outcome=="FAIL") attr(samples,"auto") <- NA else
    attr(samples,"auto") <- n.try
  samples
}



get.p.vector <- function(x,value=1){
  # Convenience function to get a named p.vector from a fit object
  # By default returns with all entries = value
  out <- x$theta[1,,1]
  out[1:length(out)] <- rep(value,length(out))
  out
}


### PDA ----

pda.dmc <- function(rts,srt,n=length(srt),adjust=.8) {
    minrt <- min(rts)
    maxrt <- max(rts)
    # Have decent number of simulated and range overlaps data.
    if ( (length(srt)<10) || (min(srt)>maxrt) || (max(srt)<minrt) )  
      return(numeric(length(rts)))
    bw <- bw.nrd0(srt)
    d <- density(srt,bw=bw,adjust=adjust,
      from=minrt-3*bw,to=maxrt+3*bw) # Do smooth over region of data
    d$y[d$y<0] <- 0
    d$y <- d$y*length(srt)/n
    out <- numeric(length(rts))
    ok <- (rts>d$x[1]) & (rts<d$x[length(d$x)])
    out[ok] <- approx(d$x,d$y,rts[ok])$y
    out[is.na(out) | !is.finite(out)] <- 0
    out
}
# 
# pda.dmc <- function(rts,srt,n=length(srt),adjust=.8) {
#     minrt <- min(rts)
#     maxrt <- max(rts)
#     # Have decent number of simulated and range overlaps data.
#     if ( (length(srt)<10) || (min(srt)>maxrt) || (max(srt)<minrt) )  
#       return(numeric(length(rts)))
#     srt <- log(srt); rts <- log(rts)
#     bw <- bw.nrd0(srt)
#     d <- density(srt,bw=bw,adjust=adjust)
#     d$y[d$y<0] <- 0
#     d$y <- d$y*length(srt)/n
#     out <- numeric(length(rts))
#     ok <- (rts>d$x[1]) & (rts<d$x[length(d$x)])
#     out[ok] <- approx(d$x,d$y,rts[ok])$y
#     out[is.na(out) | !is.finite(out)] <- 0
#     out
# }
# 
# pda.dmc <- function(rts,srt,n=length(srt),adjust=.8) {
#     minrt <- min(rts)
#     maxrt <- max(rts)
#     # Have decent number of simulated and range overlaps data.
#     if ( (length(srt)<10) || (min(srt)>maxrt) || (max(srt)<minrt) )  
#       return(numeric(length(rts)))
#     bw <- bw.nrd0(srt)
#     d <- density(srt,bw=bw,adjust=adjust) 
#     d$y[d$y<0] <- 0
#     d$y <- d$y*length(srt)/n
#     out <- numeric(length(rts))
#     ok <- (rts>d$x[1]) & (rts<d$x[length(d$x)])
#     out[ok] <- approx(d$x,d$y,rts[ok])$y
#     out[is.na(out) | !is.finite(out)] <- 0
#     out
# }


