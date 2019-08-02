# System functions for the DMC (Dynamic Models of Choice)
#    Functions for hierarchical models (versions of non-hierarchical functions 
#    and hierarchical model-specific functions)
#    Usually user does not need to edit

### Simulate data

# SSD=Inf; staircase=NA
h.simulate.dmc <- function(model,ps=NA,ns=1,n=2,
                           SSD=Inf,p.prior=NA,staircase=NA,
                           subject.cv=NULL)
  # Like simulate.dmc, but creates data for a set of ns subjects, either as fixed 
  # effects for from a hierarchical model, in which case p.prior must be supplied,
  # and ps (like p.vector in simulate.dmc) is ignored.
  # ps can be a vector exactly like p.vector, in which case each subject has 
  #   identical parameters. Can be a matrix with one row per subject, in which
  #   case must have ns rows. Saved (in expanded form) as "parameters" attribute. 
  # p.prior is a list of distributions from which subject parameters are
  #   samples. It should be created with prior.p.dmc. Saved as "p.prior" attribute 
  # ns is number of subjects to be simulated
  # n can be a single number for a balanced design or matrix for an unbalanced
  #   design, where rows are subjects and columns are design cells. If the matrix
#   has one row then all subjects have the same n in each cell, if it has one
#   column then all cells have the same n, otherwise each entry specifies the n
#   for a particular design subject x design cell combination
# SSD is for use only with stop-signal designs, specified like n except a vector
#   form is also allowed that is the same length as the data. It must have Inf 
#   in all go cells. staircase modifies SSD (see simulate.dmc)
# subject.cv: data frame, column names are in names(p.prior)
# model can be a list of models, one for each subject.

{
  # Check n
  ncell <- prod(unlist(lapply(attr(model,"factors"),length)))
  if (!is.matrix(n) & (is.vector(n) & length(n)!=1)) 
    stop("n must be a scalar or matrix")    
  if (is.matrix(n)) 
  {
    if (dim(n)[1]==1) # only cells differ 
      n <- matrix(rep(n,each=ns),nrow=ns) 
    if (dim(n)[2]==1) # only subjects differ 
      n <- matrix(rep(n,times=ncell),nrow=ns) 
  } else n <- matrix(rep(n,ns*ncell),nrow=ns)      
  if ( ns != dim(n)[1] )
    stop(paste("The n matrix must have",ns,"rows"))      
  if ( ncell != dim(n)[2] )
    stop(paste("The n matrix must have",ncell,"columns"))      
  
  # Check SSD
  ndata <- sum(n) 
  if ( is.matrix(SSD) ) 
  {
    if ( dim(SSD)[1]==1 ) # only cells differ 
      SSD <- matrix(rep(SSD,each=ns),nrow=ns) 
    if ( dim(SSD)[2]==1 ) # only subjects differ 
      SSD <- matrix(rep(SSD,times=ncell),nrow=ns) 
    if ( ns != dim(SSD)[1] )
      stop(paste("The SSD matrix must have",ns,"rows"))      
    if ( ncell != dim(SSD)[2] )
      stop(paste("The SSD matrix must have",ncell,"columns"))    
    SSD <- rep(t(SSD),as.vector(t(n)))    
  } else if ( length(SSD)==1 ) SSD <- rep(SSD,ndata)  
  if ( length(SSD) != ndata )
    stop(paste("SSD vector is not the same length as the data (",ndata,")"))
  
  # check ps/p.prior
  if ( any(is.na(p.prior)) ) 
  {
    if ( !is.matrix(ps) ) 
    {
      if (check.p.vector(ps,model)) stop() 
      ps <- matrix(rep(ps,each=ns),nrow=ns,dimnames=list(NULL,names(ps)))
    }
    if ((ns != dim(ps)[1]))
      stop("ps matrix must have ns rows")
    if (check.p.vector(ps[1,],model)) stop()
  } else { # Check p.prior
    if (!all(sort(names(attr(model,"p.vector")))==sort(names(p.prior))))
      stop("p.prior incompatible with model")   
  }
  
  # random effects
  if ( !any(is.na(p.prior)) ) 
    ps <- as.matrix(rprior.dmc(p.prior,ns))
  
  if (is.list(model)) { # Different model for each subject
    if (length(model) != ns) 
      stop("List of models must be same length as number of subjects")
    subject.models <- TRUE
  } else subject.models <- FALSE
  
  # Subject covariates
  if ( !is.null(subject.cv) ) {
    if ( !is.data.frame(subject.cv) || (dim(subject.cv)[1]!=ns) || 
      !all(names(subject.cv %in% dimnames(ps)[[2]])) )
      stop("subject.cv must be a data frame with names from p.prior and ns rows")
    for (i in names(subject.cv)) ps[,i] <- ps[,i] + subject.cv[[i]]
  }
  row.names(ps) <- 1:ns
  if (subject.models) modeli <- model[[1]] else modeli <- model
  ndatai <- cumsum(c(0,apply(n,1,sum)))
  datr <- (ndatai[1]+1):(ndatai[2])
  data <- cbind(s=rep(1,length(datr)),
    simulate.dmc(ps[1,],modeli,n[1,],SSD=SSD[datr],staircase = staircase))  
  if (ns>1) for (i in 2:ns)
  {
    if (subject.models) modeli <- model[[i]] else modeli <- model
    datr <- (ndatai[i]+1):(ndatai[i+1])
    data <- rbind(data,cbind(s=rep(i,length(datr)),
      simulate.dmc(ps[i,],modeli,n[i,],SSD=SSD[datr],staircase = staircase)))  
  }  
  data$s <- factor(data$s)
  attr(data,"parameters") <- ps
  if (!any(is.na(p.prior))) attr(data,"p.prior") <- p.prior
  data
}


### Hierarchical prior, likelihoods and posterior

h.summed.log.prior=function (pp, pp.prior)
{
  suppressWarnings( sum(log.prior.dmc(pp[[1]],pp.prior[[1]])) + 
                    sum(log.prior.dmc(pp[[2]],pp.prior[[2]])) )
}  


h.log.likelihood.dmc <- function(ps,pp,p.prior)
  # log-likelihood of subject parameters ps under population distribuiton p.prior  
{
  suppressWarnings( apply(apply(ps,1,log.prior.dmc,
                                p.prior=assign.pp(pp,p.prior)),2,sum) )  
}


h.log.posterior.dmc <- function(ps,p.prior,pp,pp.prior)
  # log-likelihood of subject parameters ps under population distribution 
  # p.prior,given population parameters pp (phi) with priors pp.prior
{
  # sum over subjects of likelihood of pp (pop pars) given ps (subjects pars)  
  sum(h.log.likelihood.dmc(ps,pp,p.prior)) +
  # prior probability of pp
  h.summed.log.prior(pp,pp.prior) 
}


### Setup samples


make.hstart <- function(fsamples,
  mu.sd=1,sigma.sd=1,verbose=TRUE,
  lower.mu=NULL,upper.mu=NULL,
  lower.sigma=NULL,upper.sigma=NULL,
  show.sd=FALSE,do.plot=FALSE,layout=c(2,5),digits=2)
  # Uses the results of fixed effects fitting to get a start prior for hierarchcial
  # By default sets prior sd to 1% of mean (tight), can specify as a vector
  # Prints prior parameters and can plot priors 
{

  mns <- lapply(fsamples,function(x){apply(x$theta,2,mean)})
  mns <- matrix(unlist(mns),nrow=length(mns[[1]]),dimnames=list(names(mns[[1]]),NULL))

  mn <- apply(mns,1,mean)
  if (length(mu.sd)==1) mu.sd <- abs(mn)*mu.sd/100
  if (length(mu.sd)!=length(mn))
    stop(paste("mu.sd must have length 1 or length of parameters:",length(mn)))
  if (is.null(lower.mu)) lower.mu <- rep(-Inf,length(mn))
  if (length(lower.mu)!=length(mn))
    stop(paste("lower.mu must have length of parameters:",length(mn)))
  if (is.null(upper.mu)) upper.mu <- rep(Inf,length(mn))
  if (length(upper.mu)!=length(mn))
    stop(paste("upper.mu must have length of parameters:",length(mn)))
  
  sd <- apply(mns,1,sd)
  if (length(sigma.sd)==1) sigma.sd <- sd*sigma.sd/100
  if (length(sigma.sd)!=length(sd))
    stop(paste("sigma.sd must have length 1 or length of parameters:",length(sd)))
  if (is.null(lower.sigma)) lower.sigma <- rep(0,length(mn))
  if (length(lower.sigma)!=length(mn))
    stop(paste("lower.sigma must have length of parameters:",length(mn)))
  if (is.null(upper.sigma)) upper.sigma <- rep(Inf,length(mn))
  if (length(upper.sigma)!=length(mn))
    stop(paste("upper.sigma must have length of parameters:",length(mn)))

  mu.prior <- prior.p.dmc(p1=mn,p2=mu.sd,lower=lower.mu,upper=upper.mu)
  sigma.prior <- prior.p.dmc(p1=sd,p2=sigma.sd,lower=lower.sigma,upper=upper.sigma)

  if (verbose) {
    cat("Mu prior\n")
    cat("Mean\n"); print(round(mn,digits))
    if (show.sd) {cat("SD\n"); print(round(mu.sd,digits))}
    cat("Sigma prior\n")
    cat("Mean\n"); print(round(sd,digits))
    if (show.sd) {cat("SD\n"); print(round(sigma.sd,digits))}
  }
  
  if (do.plot) {
    par(mfcol=layout)
    for (i in names(mu.prior)) plot.prior(i,mu.prior,main="MU")
    for (i in names(sigma.prior)) plot.prior(i,sigma.prior,main="SIGMA")
  }
  
  list(mu=mu.prior,sigma=sigma.prior)
}

make.theta1 <- function(fsamples) 
# Get last sampled thetas from each subject in fsamples, returning a list
# fsamples can also be a single subject, returns a chains x pars matrix
#
{
  
  get.theta1 <- function(samples)
    samples$theta[,,dim(samples$theta)[3]]
  
  if (is.list(fsamples)) {
    theta1 <- lapply(fsamples,get.theta1)
  } else {
    theta1 <- get.theta1(fsamples)    
  } 
  theta1
}


# cut=4;start=1;end=NA;verbose=FALSE;digits=2
pick.stuck.pars.dmc <- function(samples,cut=5,cut.v=5,start=1,end=NA) 
  # return index of stuck chains = average deviaton of value from median 
  #                                < or > cut*IQR from median  
{
  if ( is.na(end) ) end <- samples$nmc
  if ( end <= start ) stop("End must be greater than start")
  pars <- samples$theta[,,start:end]
  med <- apply(pars,2,median)
  dev <- apply(pars-rep(rep(med,each=dim(pars)[1]),times=dim(pars)[3]),1:2,mean)
  
  bad <- apply(dev,2,function(x){
    iqr <- IQR(x)
    m <- median(x)
    lo <- m - cut*iqr
    hi <- m + cut*iqr
    c(1:length(x))[x<lo | x > hi]
  })

  iqrs <- log(apply(pars,1:2,IQR))
  bad.v <- apply(iqrs,2,function(x){
    iqr <- IQR(x)
    m <- median(x)
    lo <- m - cut.v*iqr
    hi <- m + cut.v*iqr
    c(1:length(x))[x<lo | x > hi]
  })
  
  if (length(bad)!=0 & length(bad.v)!=0)
    for (i in 1:length(bad)) {
      tmp <- unique(c(bad[[i]],bad.v[[i]]))
      bad[[i]] <- tmp[!is.na(tmp)]
    }
  
  out <- bad[unlist(lapply(bad,function(x){length(x)>0}))] 
  if (length(out)==0) out <- NULL
  out
}

# start=1;end=NA;detailed=FALSE;remove.empty=TRUE
# remove.empty=FALSE
h.pick.stuck.pars.dmc <- function(hsamples,cut=5,start=1,end=NA,detailed=FALSE,remove.empty=TRUE) 
  # Do pick.stuck.pars to each subject and return non-zero cases (chain numbers by default)
{
  bad <- lapply(hsamples,pick.stuck.pars.dmc,cut=cut,start=start,end=end)
  if (remove.empty) out <- bad[unlist(lapply(bad,function(x){length(x)>0}))] else out <- bad
  if (detailed)  out else lapply(out,function(x){unique(unlist(x))})
}


# p.prior=NULL;data=NULL;pp.prior=NULL;
# samples=NULL;thin=NA;theta1=NULL;phi1=NULL
# start.prior=NULL;hstart.prior=NULL
# add=FALSE;remove=NA;rp=.001
# replace.bad.chains=FALSE;cut=5
# n.chains=ifelse(is.null(p.prior),length(attr(attr(data[[1]],
#   "model"),"p.vector"))*3,length(p.prior)*3)
# 
# nmc=100;samples=hsamples;replace.bad.chains=TRUE;cut=1

# New version with replace bad chains
h.samples.dmc <- function(nmc,p.prior=NULL,data=NULL,pp.prior=NULL,
                          samples=NULL,thin=NA,
                          theta1=NULL,phi1=NULL,
                          start.prior=NULL,hstart.prior=NULL,
                          add=FALSE,remove=NA,rp=.001,
                          replace.bad.chains=FALSE,cut=5,
                          n.chains=ifelse(is.null(p.prior),
                            length(attr(attr(data[[1]],"model"),"p.vector"))*3,
                            length(p.prior)*3)                       
)
  # Setup for hierarchical sampling
{
  
  do.samples.dmc <- function(s,nmc,p.prior,data,samples,theta1,start.prior=NULL,
                             add,rp,n.chains,thin,remove,replace.bad.chains)
  {
    if ( !add & is.null(samples) ) cat(".")
    samples.dmc(nmc=nmc,p.prior=p.prior,data=data[[s]],start.prior=start.prior,
      samples=samples[[s]],theta1=theta1[[s]][1:n.chains,names(p.prior)],add=add,
      rp=rp,verbose=FALSE,n.chains=n.chains,thin=thin,remove=remove,
      replace.bad.chains=replace.bad.chains[[s]])
  }
  
  if ( !is.null(data) && is.data.frame(data) )
    stop("For hierarchical case data cannot be a data frame")
  if ( !is.null(data) && !is.list(data) )
    stop("data arguement must be a list")
  
  if ( !is.null(samples) ) {
    if (!is.null(pp.prior) & add) stop("Cannot change pp.prior when adding!")
    if ( is.null(pp.prior) ) pp.prior <- attr(samples,"hyper")$pp.prior
    if ( is.null(p.prior)  )  p.prior <- samples[[1]]$p.prior
    if (is.na(thin)) thin <- samples[[1]]$thin
    if (replace.bad.chains) {
      replace.bad.chains <- h.pick.stuck.pars.dmc(samples,cut=cut,remove.empty = FALSE)
      if (!all(unlist(lapply(replace.bad.chains,is.null)))) {
        cat("Replacing bad chains\n")
        print(replace.bad.chains[unlist(lapply(replace.bad.chains,function(x){length(x)>0}))])
      }
    } else replace.bad.chains <- NULL
  }

  
  if ( is.null(data) ) 
  {
    n.subjects <- length(samples)
    s.names <- names(samples)
    if (!is.null(replace.bad.chains) && length(replace.bad.chains)!=n.subjects)
      stop("replace.bad.chains must be a list of the same length as samples")
  } else {
    n.subjects <- length(data)
    s.names <- names(data)
  }
  if ( is.list(theta1) && !all(sort(names(theta1))==sort(s.names)) )
    stop("Names of theta1 list do not match subject names in data")
    
   
  if ( is.null(pp.prior) ) { # non-hierarcical
    
    if ( !add & is.null(samples) ) 
      cat("Generating start points for each subject: ")
    out <- lapply(1:n.subjects,do.samples.dmc,nmc=nmc,p.prior=p.prior,add=add,
      rp=rp,start.prior=start.prior,data=data,samples=samples,theta1=theta1,
      n.chains=n.chains,thin=thin,remove=remove,
      replace.bad.chains=replace.bad.chains)
    names(out) <- s.names
    if(!add & is.null(samples)) cat("\n")    
  } else {                  # hierarchical model
    # check pp.prior
    if (!is.list(pp.prior)) stop("pp.prior must be a list")
    if ( length(pp.prior[[1]])<length(pp.prior[[2]]) )
      stop("Hyper-prior must have as many or more location than scale parameters")
    # Group parameters, !has.sigma, called fixed
    has.sigma <- names(pp.prior[[1]]) %in% names(pp.prior[[2]])
    fixed <- names(pp.prior[[1]])[ !has.sigma ]
    pp.names <- lapply(pp.prior,names)
    # Subject parameter has a hyper sigma
    has.hyper <- names(p.prior) %in% pp.names[[2]]
    if ( !all(names(p.prior)[has.hyper] %in% pp.names[[1]][has.sigma]) ) 
      stop("pp.prior location parameters not compatible with p.prior")
    if ( !all(names(p.prior)[has.hyper]==pp.names[[2]]) ) 
      stop("pp.prior scale parameters not compatible with p.prior")

    if ( !is.null(hstart.prior) ) { # check hstart.prior
      if (!is.list(hstart.prior)) stop("hstart.prior must be a list")
      if (!all(sort(names(pp.prior[[1]]))==sort(names(hstart.prior[[1]]))))
        stop("Incompatible location hstart.prior") 
      if (!all(sort(names(pp.prior[[2]]))==sort(names(hstart.prior[[2]])[has.sigma])))
        stop("Incompatible scale hstart.prior") 
    }
    if ( is.null(samples) ) { # setup new samples
      if ( is.null(data) )
        stop("Must specify data argument if samples argument not given.")
      bad <- fixed[!(fixed %in% names(attr(attr(data[[1]],"model"),"constants")))]
      if ( length(bad)>0 ) 
        stop(paste("Fixed hyper parameters not constants in data level:",bad))
      if (is.na(thin)) thin <- 1
      p.names <- names(pp.prior[[1]])
      n.pars <- length(p.names)
      phi <- array(NA,c(n.chains,n.pars,nmc))
      dimnames(phi)[[2]] <- p.names
      phi <- list(phi,phi)  # pairs of sampled hyper-parameters
      h_summed_log_prior <- array(-Inf,c(nmc,n.chains)) # hyper log-likelihoods
      h_log_likelihoods <- array(-Inf,c(nmc,n.chains))  # hyper log-likelihoods  
      if ( is.null(phi1) ) {
        cat("Generating hyper-start points for each chain: ")
        for( i in 1:n.chains ) { # ASSUMES PRIOR PRODUCES GOOD VALUES !!!
          cat(".")
          if ( is.null(hstart.prior) ) { # sample from prior
            phi[[1]][i,,1] <- rprior.dmc(pp.prior[[1]])[,p.names]  
            phi[[2]][i,has.sigma,1] <- 
              rprior.dmc(pp.prior[[2]])[,p.names[has.sigma]]  
          } else {                        # sample from hstart.prior
            phi[[1]][i,,1] <- rprior.dmc(hstart.prior[[1]])[,p.names] 
            phi[[2]][i,has.sigma,1] <- 
              rprior.dmc(hstart.prior[[2]])[,p.names[has.sigma]] 
          }
        }
        cat("\n")
      } else {
        if ( !is.list(theta1) )
          stop("If phi1 specified theta1 must be a list of thetas for each subject")
        phi[[1]][,p.names,1] <- phi1[[1]][,p.names]
        phi[[2]][,p.names,1] <- phi1[[2]][,p.names]
      }
      for (i in 1:n.chains) {
        h_summed_log_prior[1,i] <- 
            sum(log.prior.dmc(phi[[1]][i,,1],pp.prior[[1]])) + 
            sum(log.prior.dmc(phi[[2]][i,has.sigma,1],pp.prior[[2]]))
      }

      cat("Generating samples objects for each subject: ")
      out <- lapply(1:n.subjects,do.samples.dmc,nmc=nmc,p.prior=p.prior,
               data=data,samples=samples,add=add,rp=rp,n.chains=n.chains,
               theta1=theta1,thin=thin,remove=remove,
               replace.bad.chains=bad.chains)
      names(out) <- s.names
      
      # Fill in hyper-likelihoods
      cps <- lapply(out,function(x){x$theta[,,1]}) # subject list: chains x pars
      for( i in 1:n.chains ) {
        h_log_likelihoods[1,i] <- sum(h.log.likelihood.dmc(p.prior=p.prior[has.hyper],
          ps=t(data.frame(lapply(cps,function(x){x[i,]}))),                             
          pp=list(phi[[1]][i,has.sigma,1],phi[[2]][i,has.sigma,1]) ))
      }
      if ( any(!is.finite(h_log_likelihoods[1,])) )
        stop("Sampling start points was valid for data but not hyper level\n
              Try a tighter pp.prior or hstart.prior")
      cat("\n")
      
      # Correct data level prior given hypers
      consts <- phi[[1]][,,1][,!has.sigma,drop=FALSE]
      consts.names <- dimnames(consts)[[2]]
      pp.priors <- pp.prior[[1]][consts.names]
      for (k in 1:n.chains) {
        p.priork <- assign.pp(list(phi[[1]][k,has.sigma,1],
          phi[[2]][k,has.sigma,1]),p.prior=out[[1]]$p.prior[has.hyper]) 
        if ( any(!has.hyper) ) for (i in names(out[[1]]$p.prior)[!has.hyper])
          p.priork[[i]] <- out[[1]]$p.prior[[i]]
        for (j in 1:length(out))
          out[[j]]$summed_log_prior[1,k] <- summed.log.prior (
            p.vector=out[[j]]$theta[k,,1],p.prior=p.priork,
            p.names=names(out[[j]]$theta[k,,1])) + 
            summed.log.prior(consts[k,],pp.priors,consts.names)
      }

      attr(out,"hyper") <- list(phi=phi,h_log_likelihoods=h_log_likelihoods,
                                h_summed_log_prior=h_summed_log_prior,
                                pp.prior=pp.prior,start=1,n.pars=n.pars,
                                p.names=p.names,rp=rp,nmc=nmc,n.chains=n.chains,
                                thin=thin,has.sigma=has.sigma,has.hyper=has.hyper)
    } else { # restart from previous sampling
      
      # remove hyper samples 
      if ( !all(is.na(remove)) ) {
        attr(samples,"hyper")$phi[[1]] <- attr(samples,"hyper")$phi[[1]][,,-remove]
        attr(samples,"hyper")$phi[[2]] <- attr(samples,"hyper")$phi[[2]][,,-remove]
        attr(samples,"hyper")$nmc <- dim(attr(samples,"hyper")$phi[[1]])[3]
        attr(samples,"hyper")$h_summed_log_prior <- 
          attr(samples,"hyper")$h_summed_log_prior[-remove,]
        attr(samples,"hyper")$h_log_likelihoods  <- 
          attr(samples,"hyper")$h_log_likelihoods[-remove,]
      }

      out <- lapply(1:n.subjects,do.samples.dmc,nmc=nmc,data=data,
        samples=samples,theta1=theta1,p.prior=p.prior,add=add,rp=rp,
        n.chains=n.chains,thin=thin,remove=remove,replace.bad.chains=replace.bad.chains) 
      names(out) <- s.names
      
      hyper <- attr(samples,"hyper")
      if ( add ) { # add on
        phi <- array(NA,c(hyper$n.chains,hyper$n.pars,hyper$nmc+nmc))
        dimnames(phi)[[2]] <- hyper$p.names
        phi <- list(phi,phi)
        phi[[1]][,,1:hyper$nmc] <- hyper$phi[[1]]
        phi[[2]][,,1:hyper$nmc] <- hyper$phi[[2]]
        hyper$phi <- phi
        h_summed_log_prior <- array(-Inf,c(hyper$nmc+nmc,hyper$n.chains)) 
        h_log_likelihoods <- array(-Inf,c(hyper$nmc+nmc,hyper$n.chains)) 
        h_summed_log_prior[1:hyper$nmc,] <- hyper$h_summed_log_prior
        h_log_likelihoods[1:hyper$nmc,] <- hyper$h_log_likelihoods
        hyper$h_summed_log_prior <- h_summed_log_prior
        hyper$h_log_likelihoods <- h_log_likelihoods
        hyper$start <- hyper$nmc
        hyper$nmc <- hyper$nmc+nmc
      } else { # start afresh       
        old.nmc <- dim(hyper$h_log_likelihoods)[1]
        phi1 <- list(hyper$phi[[1]][,,old.nmc],hyper$phi[[2]][,,old.nmc])
        h_log_likelihoods1 <- hyper$h_log_likelihoods[old.nmc,]
        h_summed_log_prior1 <- hyper$h_summed_log_prior[old.nmc,]
        phi=array(NA,c(hyper$n.chains,hyper$n.pars,nmc))
        dimnames(phi)[[2]] <- hyper$p.names
        hyper$phi <- list(phi,phi)
        hyper$phi[[1]][,,1] <- phi1[[1]]
        hyper$phi[[2]][,,1] <- phi1[[2]]
        hyper$h_summed_log_prior <- array(-Inf,c(nmc,hyper$n.chains)) 
        hyper$h_log_likelihoods <- array(-Inf,c(nmc,hyper$n.chains)) 
        hyper$h_log_likelihoods[1,] <- h_log_likelihoods1
        # Recompute prior in case new pp.prior was provided.
        for (i in 1:n.chains) {
          hyper$h_summed_log_prior[1,i] <- 
            sum(log.prior.dmc(hyper$phi[[1]][i,,1],pp.prior[[1]])) + 
            sum(log.prior.dmc(hyper$phi[[2]][i,has.sigma,1],pp.prior[[2]]))
        }
        if (any(is.na(hyper$h_summed_log_prior[1,i])))
            stop("New pp.prior not compatible with old samples.")
        hyper$pp.prior <- pp.prior
        hyper$nmc <- nmc
        hyper$start <- 1
        hyper$thin <- thin
      }
      attr(out,"hyper") <- hyper
    }
  }
  out
}


# # p.prior=NULL;data=NULL;pp.prior=NULL; add=FALSE;remove=NA;rp=.001
# # samples=NULL;thin=NA;theta1=NULL;phi1=NULL;start.prior=NULL;hstart.prior=NULL
# # replace.bad.chains=NULL
# #                           
# # 
# # nmc=100;p.prior=p.priorV;pp.prior=pp.priorVVtight;samples=hVVl
# # 
# # n.chains=ifelse(is.null(p.prior),
# #     length(attr(attr(data[[1]],"model"),"p.vector"))*3,
# #     length(p.prior)*3)
# 
# h.samples.dmc <- function(nmc,p.prior=NULL,data=NULL,pp.prior=NULL,
#                           samples=NULL,thin=NA,
#                           theta1=NULL,phi1=NULL,
#                           start.prior=NULL,hstart.prior=NULL,
#                           add=FALSE,remove=NA,replace.bad.chains=NULL,rp=.001,
#                           n.chains=ifelse(is.null(p.prior),
#                             length(attr(attr(data[[1]],"model"),"p.vector"))*3,
#                             length(p.prior)*3)                       
# )
#   # Setup for hierarchical sampling
# {
#   
#   do.samples.dmc <- function(s,nmc,p.prior,data,samples,theta1,start.prior=NULL,
#                              add,rp,n.chains,thin,remove,replace.bad.chains)
#   {
#     if ( !add & is.null(samples) ) cat(".")
#     samples.dmc(nmc=nmc,p.prior=p.prior,data=data[[s]],start.prior=start.prior,
#       samples=samples[[s]],theta1=theta1[[s]][1:n.chains,names(p.prior)],add=add,
#       rp=rp,verbose=FALSE,n.chains=n.chains,thin=thin,remove=remove,
#       replace.bad.chains=replace.bad.chains[[s]])
#   }
#   
#   if ( !is.null(data) && is.data.frame(data) )
#     stop("For hierarchical case data cannot be a data frame")
#   if ( !is.null(data) && !is.list(data) )
#     stop("data arguement must be a list")
#   if ( is.null(data) ) 
#   {
#     n.subjects <- length(samples)
#     s.names <- names(samples)
#     if (!is.null(replace.bad.chains) && length(replace.bad.chains)!=n.subjects)
#       stop("replace.bad.chains must be a list of the same length as samples")
#   } else {
#     n.subjects <- length(data)
#     s.names <- names(data)
#   }
#   if ( is.list(theta1) && !all(sort(names(theta1))==sort(s.names)) )
#     stop("Names of theta1 list do not match subject names in data")
#     
#   if ( !is.null(samples) ) {
#     if (!is.null(pp.prior) & add) stop("Cannot change pp.prior when adding!")
#     if ( is.null(pp.prior) ) pp.prior <- attr(samples,"hyper")$pp.prior
#     if ( is.null(p.prior)  )  p.prior <- samples[[1]]$p.prior
#     if (is.na(thin)) thin <- samples[[1]]$thin
#   }
#   
#   if ( is.null(pp.prior) ) { # non-hierarcical
#     
#     if ( !add & is.null(samples) ) 
#       cat("Generating start points for each subject: ")
#     out <- lapply(1:n.subjects,do.samples.dmc,nmc=nmc,p.prior=p.prior,add=add,
#       rp=rp,start.prior=start.prior,data=data,samples=samples,theta1=theta1,
#       n.chains=n.chains,thin=thin,remove=remove,replace.bad.chains=replace.bad.chains)
#     names(out) <- s.names
#     if(!add & is.null(samples)) cat("\n")    
#   } else {                  # hierarchical model
#     # check pp.prior
#     if (!is.list(pp.prior)) stop("pp.prior must be a list")
#     if ( length(pp.prior[[1]])<length(pp.prior[[2]]) )
#       stop("Hyper-prior must have as many or more location than scale parameters")
#     # Group parameters, !has.sigma, called fixed
#     has.sigma <- names(pp.prior[[1]]) %in% names(pp.prior[[2]])
#     fixed <- names(pp.prior[[1]])[ !has.sigma ]
#     pp.names <- lapply(pp.prior,names)
#     # Subject parameter has a hyper sigma
#     has.hyper <- names(p.prior) %in% pp.names[[2]]
#     if ( !all(names(p.prior)[has.hyper] %in% pp.names[[1]][has.sigma]) ) 
#       stop("pp.prior location parameters not compatible with p.prior")
#     if ( !all(names(p.prior)[has.hyper]==pp.names[[2]]) ) 
#       stop("pp.prior scale parameters not compatible with p.prior")
# 
#     if ( !is.null(hstart.prior) ) { # check hstart.prior
#       if (!is.list(hstart.prior)) stop("hstart.prior must be a list")
#       if (!all(sort(names(pp.prior[[1]]))==sort(names(hstart.prior[[1]]))))
#         stop("Incompatible location hstart.prior") 
#       if (!all(sort(names(pp.prior[[2]]))==sort(names(hstart.prior[[2]])[has.sigma])))
#         stop("Incompatible scale hstart.prior") 
#     }
#     if ( is.null(samples) ) { # setup new samples
#       if ( is.null(data) )
#         stop("Must specify data argument if samples argument not given.")
#       bad <- fixed[!(fixed %in% names(attr(attr(data[[1]],"model"),"constants")))]
#       if ( length(bad)>0 ) 
#         stop(paste("Fixed hyper parameters not constants in data level:",bad))
#       if (is.na(thin)) thin <- 1
#       p.names <- names(pp.prior[[1]])
#       n.pars <- length(p.names)
#       phi <- array(NA,c(n.chains,n.pars,nmc))
#       dimnames(phi)[[2]] <- p.names
#       phi <- list(phi,phi)  # pairs of sampled hyper-parameters
#       h_summed_log_prior <- array(-Inf,c(nmc,n.chains)) # hyper log-likelihoods
#       h_log_likelihoods <- array(-Inf,c(nmc,n.chains))  # hyper log-likelihoods  
#       if ( is.null(phi1) ) {
#         cat("Generating hyper-start points for each chain: ")
#         for( i in 1:n.chains ) { # ASSUMES PRIOR PRODUCES GOOD VALUES !!!
#           cat(".")
#           if ( is.null(hstart.prior) ) { # sample from prior
#             phi[[1]][i,,1] <- rprior.dmc(pp.prior[[1]])[,p.names]  
#             phi[[2]][i,has.sigma,1] <- 
#               rprior.dmc(pp.prior[[2]])[,p.names[has.sigma]]  
#           } else {                        # sample from hstart.prior
#             phi[[1]][i,,1] <- rprior.dmc(hstart.prior[[1]])[,p.names] 
#             phi[[2]][i,has.sigma,1] <- 
#               rprior.dmc(hstart.prior[[2]])[,p.names[has.sigma]] 
#           }
#         }
#         cat("\n")
#       } else {
#         if ( !is.list(theta1) )
#           stop("If phi1 specified theta1 must be a list of thetas for each subject")
#         phi[[1]][,p.names,1] <- phi1[[1]][,p.names]
#         phi[[2]][,p.names,1] <- phi1[[2]][,p.names]
#       }
#       for (i in 1:n.chains) {
#         h_summed_log_prior[1,i] <- 
#             sum(log.prior.dmc(phi[[1]][i,,1],pp.prior[[1]])) + 
#             sum(log.prior.dmc(phi[[2]][i,has.sigma,1],pp.prior[[2]]))
#       }
# 
#       cat("Generating samples objects for each subject: ")
#       out <- lapply(1:n.subjects,do.samples.dmc,nmc=nmc,p.prior=p.prior,
#                data=data,samples=samples,add=add,rp=rp,n.chains=n.chains,
#                theta1=theta1,thin=thin,remove=remove,replace.bad.chains=replace.bad.chains)
#       names(out) <- s.names
#       
#       # Fill in hyper-likelihoods
#       cps <- lapply(out,function(x){x$theta[,,1]}) # subject list: chains x pars
#       for( i in 1:n.chains ) {
#         h_log_likelihoods[1,i] <- sum(h.log.likelihood.dmc(p.prior=p.prior[has.hyper],
#           ps=t(data.frame(lapply(cps,function(x){x[i,]}))),                             
#           pp=list(phi[[1]][i,has.sigma,1],phi[[2]][i,has.sigma,1]) ))
#       }
#       if ( any(!is.finite(h_log_likelihoods[1,])) )
#         stop("Sampling start points was valid for data but not hyper level\n
#               Try a tighter pp.prior or hstart.prior")
#       cat("\n")
#       
#       # Correct data level prior given hypers
#       consts <- phi[[1]][,,1][,!has.sigma,drop=FALSE]
#       consts.names <- dimnames(consts)[[2]]
#       pp.priors <- pp.prior[[1]][consts.names]
#       for (k in 1:n.chains) {
#         p.priork <- assign.pp(list(phi[[1]][k,has.sigma,1],
#           phi[[2]][k,has.sigma,1]),p.prior=out[[1]]$p.prior[has.hyper]) 
#         if ( any(!has.hyper) ) for (i in names(out[[1]]$p.prior)[!has.hyper])
#           p.priork[[i]] <- out[[1]]$p.prior[[i]]
#         for (j in 1:length(out))
#           out[[j]]$summed_log_prior[1,k] <- summed.log.prior (
#             p.vector=out[[j]]$theta[k,,1],p.prior=p.priork,
#             p.names=names(out[[j]]$theta[k,,1])) + 
#             summed.log.prior(consts[k,],pp.priors,consts.names)
#       }
# 
#       attr(out,"hyper") <- list(phi=phi,h_log_likelihoods=h_log_likelihoods,
#                                 h_summed_log_prior=h_summed_log_prior,
#                                 pp.prior=pp.prior,start=1,n.pars=n.pars,
#                                 p.names=p.names,rp=rp,nmc=nmc,n.chains=n.chains,
#                                 thin=thin,has.sigma=has.sigma,has.hyper=has.hyper)
#     } else { # restart from previous sampling
#       
#       # remove hyper samples 
#       if ( !all(is.na(remove)) ) {
#         attr(samples,"hyper")$phi[[1]] <- attr(samples,"hyper")$phi[[1]][,,-remove]
#         attr(samples,"hyper")$phi[[2]] <- attr(samples,"hyper")$phi[[2]][,,-remove]
#         attr(samples,"hyper")$nmc <- dim(attr(samples,"hyper")$phi[[1]])[3]
#         attr(samples,"hyper")$h_summed_log_prior <- 
#           attr(samples,"hyper")$h_summed_log_prior[-remove,]
#         attr(samples,"hyper")$h_log_likelihoods  <- 
#           attr(samples,"hyper")$h_log_likelihoods[-remove,]
#       }
# 
#       out <- lapply(1:n.subjects,do.samples.dmc,nmc=nmc,data=data,
#         samples=samples,theta1=theta1,p.prior=p.prior,add=add,rp=rp,
#         n.chains=n.chains,thin=thin,remove=remove,replace.bad.chains=replace.bad.chains)
#       names(out) <- s.names
#       hyper <- attr(samples,"hyper")
#       if ( add ) { # add on
#         phi <- array(NA,c(hyper$n.chains,hyper$n.pars,hyper$nmc+nmc))
#         dimnames(phi)[[2]] <- hyper$p.names
#         phi <- list(phi,phi)
#         phi[[1]][,,1:hyper$nmc] <- hyper$phi[[1]]
#         phi[[2]][,,1:hyper$nmc] <- hyper$phi[[2]]
#         hyper$phi <- phi
#         h_summed_log_prior <- array(-Inf,c(hyper$nmc+nmc,hyper$n.chains)) 
#         h_log_likelihoods <- array(-Inf,c(hyper$nmc+nmc,hyper$n.chains)) 
#         h_summed_log_prior[1:hyper$nmc,] <- hyper$h_summed_log_prior
#         h_log_likelihoods[1:hyper$nmc,] <- hyper$h_log_likelihoods
#         hyper$h_summed_log_prior <- h_summed_log_prior
#         hyper$h_log_likelihoods <- h_log_likelihoods
#         hyper$start <- hyper$nmc
#         hyper$nmc <- hyper$nmc+nmc
#       } else { # start afresh       
#         old.nmc <- dim(hyper$h_log_likelihoods)[1]
#         phi1 <- list(hyper$phi[[1]][,,old.nmc],hyper$phi[[2]][,,old.nmc])
#         h_log_likelihoods1 <- hyper$h_log_likelihoods[old.nmc,]
#         h_summed_log_prior1 <- hyper$h_summed_log_prior[old.nmc,]
#         phi=array(NA,c(hyper$n.chains,hyper$n.pars,nmc))
#         dimnames(phi)[[2]] <- hyper$p.names
#         hyper$phi <- list(phi,phi)
#         hyper$phi[[1]][,,1] <- phi1[[1]]
#         hyper$phi[[2]][,,1] <- phi1[[2]]
#         hyper$h_summed_log_prior <- array(-Inf,c(nmc,hyper$n.chains)) 
#         hyper$h_log_likelihoods <- array(-Inf,c(nmc,hyper$n.chains)) 
#         hyper$h_log_likelihoods[1,] <- h_log_likelihoods1
#         # Recompute prior in case new pp.prior was provided.
#         for (i in 1:n.chains) {
#           hyper$h_summed_log_prior[1,i] <- 
#             sum(log.prior.dmc(hyper$phi[[1]][i,,1],pp.prior[[1]])) + 
#             sum(log.prior.dmc(hyper$phi[[2]][i,has.sigma,1],pp.prior[[2]]))
#         }
#         if (any(is.na(hyper$h_summed_log_prior[1,i])))
#             stop("New pp.prior not compatible with old samples.")
#         hyper$pp.prior <- pp.prior
#         hyper$nmc <- nmc
#         hyper$start <- 1
#         hyper$thin <- thin
#       }
#       attr(out,"hyper") <- hyper
#     }
#   }
#   out
# }






### Run sampling

h.migrate <- function(use.phi,use.logprior,use.loglike,p.prior,ps,rp,pp.prior,
                      has.sigma,has.hyper,is.constant)
  # DEMCMC migrate set, all chains, hyper level
{
      
  # Which pars are not constants?
  de=list(!unlist(is.constant[[1]]),!unlist(is.constant[[2]]))
  
  n.chains <- dim(use.phi[[1]])[1]
  pars <- dim(use.phi[[1]])[2]
  lnum1 <- sample(c(1:n.chains),1)  # determine how many groups to work with
  lnum2 <- sort(sample(c(1:n.chains),lnum1,replace=F))  # get groups 
  phiset <- list( # initialize
    matrix(NA,lnum1,pars,dimnames=list(NULL,dimnames(use.phi[[1]])[[2]])),
    matrix(NA,lnum1,pars,dimnames=list(NULL,dimnames(use.phi[[1]])[[2]]))
  )	
  
  propset.logprior <- propset.loglike <- numeric(lnum1)
  currentset <- propset <- propw <- currw <- numeric(lnum1)
  
  index=numeric(lnum1)
  for (i in 1:lnum1) {
    index[i] <- sample(1:n.chains,1,replace=F)	
    # create a set of these particles to swap
    phiset[[1]][i,] <- use.phi[[1]][lnum2[i],] 
    phiset[[2]][i,] <- use.phi[[2]][lnum2[i],]  
    
    # perturb non-constant parameters
    phiset[[1]][i,de[[1]]] <- phiset[[1]][i,de[[1]]] + runif(1,-rp,rp)	
    phiset[[2]][i,de[[2]]] <- phiset[[2]][i,de[[2]]] + runif(1,-rp,rp)  
    
    propset.logprior[i] <- h.summed.log.prior (
      pp=list(phiset[[1]][i,has.sigma],phiset[[2]][i,has.sigma]),pp.prior=pp.prior)
    propset.loglike[i]  <- sum(h.log.likelihood.dmc(ps=ps[i,,], # has.hyper
      pp=list(phiset[[1]][i,has.sigma],phiset[[2]][i,has.sigma]),p.prior=p.prior[has.hyper]))
    
    propset[i] <- propset.logprior[i] + propset.loglike[i] 
    propw[i] <- propset[i]
    
    currentset[i] <- use.logprior[lnum2[i]] + use.loglike[lnum2[i]] 
    # yishin calcuatle new log likelihood!
    # # Update use.loglike for new ps
    # phi <- list(use.phi[[1]][k,],use.phi[[2]][k,])
    # use.loglike[k] <- sum(h.log.likelihood.dmc(ps=ps[k.theta,,], # has.hyper
    #   pp=list(phi[[1]][has.sigma],phi[[2]][has.sigma]),p.prior=p.prior[has.hyper]))

  

  }
  currw <- currentset

  mh <- exp(propw[lnum1] - currw[1])
  if ( !is.na(mh) && (runif(1) < mh) ) {
    use.phi[[1]][lnum2[1],] <- phiset[[1]][lnum1,]	# swap the 1st with last 
    use.phi[[2]][lnum2[1],] <- phiset[[2]][lnum1,]  #  (creating a circle)
    use.logprior[lnum2[1]] <- propset.logprior[lnum1] 
    use.loglike[lnum2[1]] <- propset.loglike[lnum1] 
  }
  if ( lnum1!=1 ) {										# make sure we are not done yet
    for(i in 1:(lnum1-1)){	
      mh <- exp(propw[i] - currw[i+1])
      if ( !is.na(mh) && (runif(1) < mh) ) {
        use.phi[[1]][lnum2[i+1],] <- phiset[[1]][i,]		
        use.phi[[2]][lnum2[i+1],] <- phiset[[2]][i,] 
        use.logprior[lnum2[i+1]] <- propset.logprior[i]
        use.loglike[lnum2[i+1]] <- propset.loglike[i] 
      }
    }
  }
  cbind(use.logprior,use.loglike,use.phi[[1]],use.phi[[2]])
}


blocked.h.crossover <- function(k,blocks,n.pars,use.phi,use.logprior,use.loglike,
  p.prior,ps,pp.prior,rp,has.sigma,has.hyper,is.constant,
  force=FALSE,gamma.mult=2.38,h.gamma.mult=2.38,random.theta=TRUE)
  # hyper level crossover for hypers with a sigma, as a series of blocks
{
  
  h.crossover <- function(k,pars,use.phi,use.logprior,use.loglike,p.prior,ps,
    is.constant,pp.prior,rp,has.sigma,has.hyper,h.gamma.mult=2.38,force=FALSE,
    random.theta=TRUE)
    # DEMCMC crossover update of one chain, phi level
  {
    
    
    if (random.theta) k.theta <- sample(1:dim(ps)[1],1) else k.theta <- k
      
    # Which pars are not constants?
    de=list(!unlist(is.constant[[1]][names(is.constant[[1]])[pars]]),
            !unlist(is.constant[[2]][names(is.constant[[1]])[pars]]))
    if ( all(!unlist(de)) ) # all are constants
      return(c(use.logprior[k],use.loglike[k],use.phi[[1]][k,],use.phi[[2]][k,]))
    
    # step size
    if ( is.na(h.gamma.mult) ) hgamma <- runif(1,0.5,1) else
      hgamma <-  h.gamma.mult/sqrt(2*length(pars)*2) # extra *2 as p1 and p2
    
    # Update use.loglike for new ps
    phi <- list(use.phi[[1]][k,],use.phi[[2]][k,])
    use.loglike[k] <- sum(h.log.likelihood.dmc(ps=ps[k.theta,,], # has.hyper
      pp=list(phi[[1]][has.sigma],phi[[2]][has.sigma]),p.prior=p.prior[has.hyper]))
    
    # Calcualte use.post
    use.post <-  use.logprior[k] + use.loglike[k]
    
    # DE step
    
    # pick two other chains
    index <- sample(c(1:dim(use.phi[[1]])[1])[-k],2,replace=F)
    
    # update mu
    phi[[1]][pars[de[[1]]]] <- use.phi[[1]][k,pars[de[[1]]]] + runif(1,-rp,rp) +
      hgamma*(use.phi[[1]][index[1],pars[de[[1]]]]-use.phi[[1]][index[2],pars[de[[1]]]]) 
    
    # update sigma
    phi[[2]][pars[de[[2]]]] <- use.phi[[2]][k,pars[de[[2]]]] + runif(1,-rp,rp) + 
      hgamma*(use.phi[[2]][index[1],pars[de[[2]]]]-use.phi[[2]][index[2],pars[de[[2]]]]) 

    # Get new post
    logprior <- h.summed.log.prior(pp.prior=pp.prior,
      pp=list(phi[[1]][has.sigma],phi[[2]][has.sigma])) # only for has.sigma 
    loglike <- sum(h.log.likelihood.dmc(ps=ps[k.theta,,], #has.hyper
      p.prior=p.prior[has.hyper],pp=list(phi[[1]][has.sigma],phi[[2]][has.sigma])))
    post <- logprior + loglike
    
    # Metropolis step
    epup <- exp(post-use.post)
    if ( force || (!is.na(epup)  && (runif(1) < epup)) ) {
      use.phi[[1]][k,pars] <- phi[[1]][pars]
      use.phi[[2]][k,pars] <- phi[[2]][pars]
      use.logprior[k] <- logprior
      use.loglike[k] <- loglike
    }
    
    c(use.logprior[k],use.loglike[k],use.phi[[1]][k,],use.phi[[2]][k,])
  }
  
  
  temp <- h.crossover(k,pars=blocks[[1]],use.phi=use.phi,force=force,
                      use.logprior=use.logprior,use.loglike=use.loglike,
                      p.prior=p.prior,ps=ps,pp.prior=pp.prior,rp=rp,
                      h.gamma.mult=h.gamma.mult,is.constant=is.constant,
                      has.sigma=has.sigma,has.hyper=has.hyper,
                      random.theta=random.theta)
  if ( length(blocks)>1 ) for ( b in 2:length(blocks) ) {
    use.logprior[k]  <- temp[1]
    use.loglike[k]   <- temp[2]
    use.phi[[1]][k,] <- temp[3:(n.pars+2)]  
    use.phi[[2]][k,] <- temp[(n.pars+3):(2*n.pars+2)]  
    temp <- h.crossover(k,pars=blocks[[b]],use.phi=use.phi,force=force,
                        use.logprior=use.logprior,use.loglike=use.loglike,
                        p.prior=p.prior,ps=ps,pp.prior=pp.prior,rp=rp,
                        h.gamma.mult=h.gamma.mult,is.constant=is.constant,
                        has.sigma=has.sigma,has.hyper=has.hyper,
                        random.theta=random.theta)
  } 
  temp
}


crossover.h <- function(k,pars,use.theta,use.logprior,use.loglike,p.priors,data,
                        rp,gamma.mult=2.38,consts=NULL,pp.priors=NULL,force=FALSE)
  # Data level crossover wity different priors for each chain, p.priors is a list
{
  # step size
  if (is.na(gamma.mult)) gamma <- runif(1,0.5,1) else
    gamma <-  gamma.mult/sqrt(2*length(pars))

  # DE step
  # pick two other chains
  index <- sample(c(1:dim(use.theta)[1])[-k],2,replace=F)
  
  # Update theta
  theta <- use.theta[k,]
  names(theta)=names(attr(attributes(data)$model,"p.vector")) 
  theta[pars] <- 
    use.theta[k,pars] + 
    gamma*(use.theta[index[1],pars]-use.theta[index[2],pars]) + 
    runif(1,-rp,rp)  
  
  # Combine use to get old post
  summed.use.post <-  use.loglike[k] + use.logprior[k]
  
  # Get prior for new point
  log.prior <- summed.log.prior (p.vector=theta,p.prior=p.priors[[k]],
    p.names=names(attr(attributes(data)$model,"p.vector")))
  
  if ( !is.null(consts) ) { # !has.sigma hyper to incorporte
    log.prior <- log.prior + summed.log.prior(
      consts[k,],pp.priors,p.names=dimnames(consts)[[2]])
    attr(attr(data,"model"),"all.par")[dimnames(consts)[[2]]] <- consts[k,]
  }
  loglike  <- sum(log.likelihood (p.vector=theta,data=data))
  
  # post for new point
  post <- log.prior + loglike
  
  # Metropolis step
  epup <- exp(post-summed.use.post)
  if (force || (!is.na(epup) && (runif(1) < epup)) ) {
    use.theta[k,]   <- theta
    use.logprior[k] <- log.prior
    use.loglike[k]  <- loglike
  }
 
  c(use.logprior[k], use.loglike[k], use.theta[k,])
}


migrate.h <- function(use.theta,use.logprior,use.loglike,
                      p.priors,data,rp,consts=NULL,pp.priors=NULL)
  # Data level migrate with different priors for each chain, p.priors is a list
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

    # Proposed new prior
    propset.logprior[i] <- summed.log.prior (p.vector=thetaset[i,], 
      p.prior=p.priors[[i]],p.names=names(attr(attributes(data)$model,"p.vector")))
    
    if ( !is.null(consts) ) { # !has.sigma hyper to incorporte
      propset.logprior[i] <- propset.logprior[i] + summed.log.prior(
        consts[i,],pp.priors,p.names=dimnames(consts)[[2]])
      attr(attr(data,"model"),"all.par")[dimnames(consts)[[2]]] <- consts[i,]
    }

    propset.loglike[i]  <- sum(log.likelihood (p.vector=thetaset[i,],data=data))
    
    propw.logprior[i]      <- propset.logprior[i]
    propw.loglike[i]       <- propset.loglike[i]
    
    currw.logprior[i]      <- currentset.logprior[i]
    currw.loglike[i]       <- currentset.loglike[i]
    
  }

  epup <- exp( (propset.loglike[n.groups] + propset.logprior[n.groups]) - 
                      (currentset.loglike[1] + currentset.logprior[1]))
  if (!is.na(epup) && (runif(1) < epup)) {
    use.theta[groups[1],]   <- thetaset[n.groups,]	# swap the 1st with last 
    use.logprior[groups[1]] <- propset.logprior[n.groups]
    use.loglike[groups[1]]  <- propset.loglike[n.groups]
  }
  
  if ( n.groups!=1 ) {										# make sure we are not done yet
    for(i in 1:(n.groups-1)) {		
      epup <- exp((propset.loglike[i] + propset.logprior[i]) - 
                  (currentset.loglike[i+1] + currentset.logprior[i+1]))
      if(!is.na(epup) && (runif(1) < epup)) {
        use.theta[groups[i+1],] <- thetaset[i,]	
        use.logprior[groups[i+1]] <- propset.logprior[i]
        use.loglike[groups[i+1]] <- propset.loglike[i]
      }
    }
  }
  
  cbind (use.logprior, use.loglike, use.theta)
}


# samples=Hsamples.ctl; cores=1;report=samples[[1]]$nmc;blocks=NA
# p.migrate=0;h.p.migrate=0;force.hyper=NA;force.data=NA
# gamma.mult=2.38;h.gamma.mult=NA;random.phi=TRUE;random.theta=TRUE


h.run.dmc <- function(samples,cores=1,report=samples[[1]]$nmc,blocks=NA,
                      p.migrate=0,h.p.migrate=0,force.hyper=NA,force.data=NA,
                      gamma.mult=2.38,h.gamma.mult=NA,random.phi=TRUE,random.theta=TRUE)
  # Fixed effects: run sampling spreading subjects over cores. For cores=1 report 
  #   is made for each subject, so default is just reprot when nmc done for each
  # Hierarchical: chains spread over cores at both subject and hyper level. 
  #   Hyper sampling is blocked as specified in blocks (hyper-parameter pairs done 
  #   together), by default one pair at a time. 
{
  
  # mci=i;thin=hyper$thin;phi=use.phi;has.sigma=hyper$has.sigma;has.hyper=hyper$has.hyper
  run.chains <- function(samplesi,mci,thin,force=FALSE,
                         phi,has.sigma,has.hyper,pp.prior,
                         p.migrate=p.migrate,gamma.mult=2.38,random.phi=TRUE)
    # Data level for hierarchical
  {  
    
    p.priors <- vector(mode="list",length=samplesi$n.chains)
    consts <- phi[[1]][,!has.sigma,drop=FALSE]
    p.names <- dimnames(consts)[[2]]
    pp.priors <- pp.prior[[1]][p.names]

    for (k in 1:samplesi$n.chains) {
      if (random.phi) k.phi <- sample(1:samplesi$n.chains,1) else k.phi <- k
      p.priors[[k]] <- assign.pp(list(phi[[1]][k.phi,has.sigma],
        phi[[2]][k.phi,has.sigma]),p.prior=samplesi$p.prior[has.hyper])
      if (any(!has.hyper)) for (i in names(samplesi$p.prior)[!has.hyper])
        p.priors[[k]][[i]] <- samplesi$p.prior[[i]]
      # Update use.logprior based on new priors
      attr(samplesi,"use")$logprior[k] <- summed.log.prior(
        p.vector=attr(samplesi,"use")$theta[k,],p.prior=p.priors[[k]],
        p.names=names(attr(attributes(samplesi$data)$model,"p.vector"))) +
        summed.log.prior(consts[k,],pp.priors,p.names)
    }

    if ( runif(1) < p.migrate ) {        # Do migration
        temp <- migrate.h(use.theta      = attr(samplesi,"use")$theta,
                      use.logprior       = attr(samplesi,"use")$logprior,
                      use.loglike        = attr(samplesi,"use")$loglike,
                      p.priors           = p.priors,
                      data               = samplesi$data,
                      rp                 = samplesi$rp,
                      consts             = consts,
                      pp.priors          = pp.priors) 
    } else {                             # Do crossover
      temp <- t(sapply(1:samplesi$n.chains,crossover.h,
                       pars=1:samplesi$n.pars,
                       use.theta    = attr(samplesi,"use")$theta,
                       use.logprior = attr(samplesi,"use")$logprior,
                       use.loglike  = attr(samplesi,"use")$loglike,   
                       p.priors     = p.priors,
                       data         = samplesi$data,
                       rp           = samplesi$rp,
                       gamma.mult   = gamma.mult,
                       consts       = consts,
                       force        = force,
                       pp.priors    = pp.priors))
    }
    
    # Harvest results                  
    attr(samplesi,"use")$logprior <- temp[,1]
    attr(samplesi,"use")$loglike  <- temp[,2]   
    attr(samplesi,"use")$theta   <- temp[,-c(1,2)]
    attr(samplesi,"use")$p.priors <- p.priors
      
    # store samples
    if ( mci %% thin == 0 ) { 
          attr(samplesi,"use")$store_i <- attr(samplesi,"use")$store_i + 1
          samplesi$summed_log_prior[attr(samplesi,"use")$store_i,] <- temp[,1]                                  
          samplesi$log_likelihoods[attr(samplesi,"use")$store_i,]  <- temp[,2]
          samplesi$theta[,,attr(samplesi,"use")$store_i]     <- temp[,-c(1,2)]
    }
    
    # gc()
    
    samplesi
  }

  no.sig.crossover <- function(k,samples,pars,rp,hgamma,pp.prior,force=FALSE,
                               use.phi,use.logpriors,use.loglikes) 
    # one chain hyper crossover for no sigma parameters  
  {
    # DE step
    # pick two other chains
    index <- sample(c(1:dim(use.phi[[1]])[1])[-k],2,replace=F)
    use.consts <- use.phi[[1]][k,pars]
    consts <- use.phi[[1]][k,pars] + runif(1,-rp,rp) +
      hgamma*(use.phi[[1]][index[1],pars]-use.phi[[1]][index[2],pars]) 
    p.names <- names(consts)
    pp.priord <- hyper$pp.prior[[1]][p.names]
    loglikes <- logpriors <- numeric(length(samples))
    for ( s in 1:length(samples) ) { 
      use.logpriors[k,s] <- summed.log.prior(
        p.vector=attr(samples[[s]],"use")$theta[k,],
        p.prior=attr(samples[[s]],"use")$p.priors[[k]],
        p.names=names(attr(attributes(samples[[s]]$data)$model,"p.vector"))
      ) + summed.log.prior(use.consts,pp.priord,p.names)
      logpriors[s] <- summed.log.prior(
        p.vector=attr(samples[[s]],"use")$theta[k,],
        p.prior=attr(samples[[s]],"use")$p.priors[[k]],
        p.names=names(attr(attributes(samples[[s]]$data)$model,"p.vector"))
      ) + summed.log.prior(consts,pp.priord,p.names)  
      datai <- samples[[s]]$data
      attr(attr(datai,"model"),"all.par")[p.names] <- use.consts
      use.loglikes[k,s] <- sum(log.likelihood(
        attr(samples[[s]],"use")$theta[k,],datai))
      attr(attr(datai,"model"),"all.par")[p.names] <- consts
      loglikes[s] <- sum(log.likelihood(
        attr(samples[[s]],"use")$theta[k,],datai))
    }
    post <- sum(logpriors + loglikes)
    use.post <- sum(use.logpriors[k,] + use.loglikes[k,])
    # Metropolis step
    epup <- exp(post-use.post)
    if ( force || (!is.na(epup)  && (runif(1) < epup)) ) {
      use.phi[[1]][k,pars] <- consts
      use.logpriors[k,] <- logpriors
      use.loglikes[k,] <- loglikes
    }
    c(use.logpriors[k,],use.loglikes[k,],use.phi[[1]][k,pars])
  }

  os <- get.os()
  hyper <- attr(samples,"hyper")
  if ( !is.null(hyper) ) # hierarchical sampling
  {
    
    # Setup
    if ( any(is.na(blocks)) ) blocks <- as.list(1:hyper$n.pars) else { # check
      if (any(unlist(lapply(blocks,function(x){
        length(x)==1 || all(hyper$has.sigma[x][1]==hyper$has.sigma[x][-1])}))))
        stop("Cant mix hyper-paramaters with and without sigma in a block")
    }
    sigma.block <- unlist(lapply(blocks,function(x){all(hyper$has.sigma[x])}))

    if (is.null(hyper$thin)) hyper$thin <- 1
    nsamp <- 1+(hyper$nmc-hyper$start)*hyper$thin
    
    if (any(is.na(force.data))) force.data <- rep(FALSE,nsamp-1)
    if (any(is.na(force.hyper))) force.hyper <- rep(FALSE,nsamp-1)
    if (!is.logical(force.hyper) || (length(force.hyper)!=(nsamp-1)))
      stop(paste("force.hyper argument must be a logical vector of length",nsamp-1))
    force.hyper <- c(TRUE,force.hyper)
    if (!is.logical(force.data) || (length(force.data)!=(nsamp-1)))
      stop(paste("force.data argument must be a logical vector of length",nsamp-1))
    force.data <- c(TRUE,force.data)

    # Setup "use" for hyper level
    ps <- aperm(array(dim=c(dim(samples[[1]]$theta)[-3],length(samples)),
      data=unlist(lapply(samples,function(x){x$theta[,,hyper$start]}),use.names=FALSE),
      dimnames=list(NULL,samples[[1]]$p.names,NULL)),c(1,3,2))
    use.phi <- list(hyper$phi[[1]][,,hyper$start],hyper$phi[[2]][,,hyper$start])
    use.logprior <- hyper$h_summed_log_prior[hyper$start,]
    use.loglike <- hyper$h_log_likelihoods[hyper$start,]
    # Setup hyper level
    store_i <- hyper$start
    pp.prior <- hyper$pp.prior
    is.constant <- lapply(pp.prior,function(x){
      lapply(x,function(y){attr(y,"dist")=="constant"})})
    
    # Setup for data level
    p.priors <- vector(mode="list",length=samples[[1]]$n.chains)
    for (k in 1:samples[[1]]$n.chains) {
      p.priors[[k]] <- assign.pp(
        list(use.phi[[1]][k,hyper$has.sigma],
             use.phi[[2]][k,hyper$has.sigma]),
             p.prior=samples[[1]]$p.prior[hyper$has.hyper])
    }
    samples <- lapply(samples,function(x){attr(x,"use") <- list(
        theta=x$theta[,,x$start],
        logprior=x$summed_log_prior[x$start,],
        loglike=x$log_likelihoods[x$start,],
        store_i = x$start,
        p.priors = p.priors
      ); x })
    p.prior=samples[[1]]$p.prior
   
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
    
    for (i in 2:nsamp) 
    {
      
      # Update data and hyper level for group parameters (hypers with no sigma)
      if ( any(!sigma.block) ) { 
        no.sig.blocks=blocks[!sigma.block]
        for ( j in 1:length(no.sig.blocks) ) {
          pars <- no.sig.blocks[[j]]
          # step size
          if ( is.na(h.gamma.mult) ) hgamma <- runif(1,0.5,1) else
            hgamma <-  h.gamma.mult/sqrt(2*length(pars)) 
          # subject list of logprior and loglike for each chain
          use.logpriors <- matrix(unlist(lapply(samples,function(x){
            attr(x,"use")$logprior})),nrow=samples[[1]]$n.chains)
          use.loglikes <-  matrix(unlist(lapply(samples,function(x){
            attr(x,"use")$loglike})),nrow=samples[[1]]$n.chains)

          if (cores > 1 & os == "windows") {
            temp <- sfLapply(1:hyper$n.chains,no.sig.crossover,samples=samples,
              pars=pars,hgamma=hgamma,pp.prior=pp.prior,use.phi=use.phi,
              force=force.hyper[i],rp=hyper$rp,
              use.logpriors=use.logpriors,use.loglikes=use.loglikes)
          } else if (cores > 1) {
            temp <- t(mclapply(1:hyper$n.chains, no.sig.crossover,
              samples=samples,pars=pars,hgamma=hgamma,pp.prior=pp.prior,
              use.phi=use.phi,force=force.hyper[i],rp=hyper$rp,
              use.logpriors=use.logpriors,
              use.loglikes=use.loglikes,mc.cores=cores))
          } else {
            temp <- lapply(1:hyper$n.chains,no.sig.crossover,samples=samples,
              pars=pars,force=force.hyper[i],hgamma=hgamma,pp.prior=pp.prior,
              use.phi=use.phi,rp=hyper$rp,
              use.logpriors=use.logpriors,use.loglikes=use.loglikes)
          }
          
          temp <- t(matrix(unlist(temp),ncol=hyper$n.chains))         
          for (s in 1:length(samples)) {
            attr(samples[[s]],"use")$logprior <- temp[,s]
            attr(samples[[s]],"use")$loglike <- temp[,s+length(samples)]
          }
          use.phi[[1]][,pars] <- temp[,-c(1:(2*length(samples))),drop=FALSE]
        }    
      }
      
      # Update hyper level with sigma
      if ( runif(1)<h.p.migrate ) {
        temp <- h.migrate(use.phi=use.phi,
                          use.logprior=use.logprior,
                          use.loglike=use.loglike,
                          p.prior=p.prior,is.constant=is.constant,
                          ps=ps,pp.prior=pp.prior,
                          rp=hyper$rp,has.hyper=hyper$has.hyper,
                          has.sigma=hyper$has.sigma) 
      } else if ( cores>1 & os == "windows") {
        temp <- t(sfLapply(1:hyper$n.chains,blocked.h.crossover,
                      blocks=blocks[sigma.block],n.pars=hyper$n.pars,
                      use.phi=use.phi,force=force.hyper[i],
                      use.logprior=use.logprior,
                      use.loglike=use.loglike,random.theta=random.theta,
                      p.prior=p.prior,ps=ps,rp=hyper$rp,
                      pp.prior=pp.prior,is.constant=is.constant,
                      has.sigma=hyper$has.sigma,has.hyper=hyper$has.hyper,
                      h.gamma.mult=h.gamma.mult))
        temp=matrix(unlist(temp,use.names=FALSE),nrow=hyper$n.chains,byrow=TRUE) 
      } else if ( cores>1 ) {
        temp <- t(mclapply(1:hyper$n.chains,blocked.h.crossover,
                  blocks=blocks[sigma.block],n.pars=hyper$n.pars,
                  use.phi=use.phi,force=force.hyper[i],
                  use.logprior=use.logprior,
                  use.loglike=use.loglike,random.theta=random.theta,
                  p.prior=p.prior,ps=ps,rp=hyper$rp,
                  pp.prior=pp.prior,is.constant=is.constant,
                  has.sigma=hyper$has.sigma,has.hyper=hyper$has.hyper,
                  h.gamma.mult=h.gamma.mult, mc.cores=cores))
        temp=matrix(unlist(temp,use.names=FALSE),nrow=hyper$n.chains,byrow=TRUE) 
      } else {
        temp <- t(sapply(1:hyper$n.chains,blocked.h.crossover,
                         blocks=blocks[sigma.block],n.pars=hyper$n.pars,
                         use.phi=use.phi,force=force.hyper[i],
                         use.logprior=use.logprior,
                         use.loglike=use.loglike,random.theta=random.theta,
                         p.prior=p.prior,ps=ps,rp=hyper$rp,
                         pp.prior=pp.prior,is.constant=is.constant,
                         has.sigma=hyper$has.sigma,has.hyper=hyper$has.hyper,
                         h.gamma.mult=h.gamma.mult))  
      }

      use.logprior <- temp[,1]
      use.loglike <- temp[,2]
      use.phi[[1]][,] <- temp[,-c(1,2)][,1:hyper$n.pars]
      use.phi[[2]][,] <- temp[,-c(1,2)][,(hyper$n.pars+1):(2*hyper$n.pars)]
      
      # Update data level for all but !has.sigma parameters
     if (cores>1 & os=="windows") {
        samples <- sfLapply(samples,run.chains,mci=i,p.migrate=p.migrate,
          thin=hyper$thin,gamma.mult=gamma.mult,phi=use.phi,pp.prior=pp.prior,
          has.sigma=hyper$has.sigma,has.hyper=hyper$has.hyper,
          force=force.data[i],random.phi=random.phi) 
      } else if (cores>1) {
        samples <- mclapply(samples,run.chains,mci=i,
          p.migrate=p.migrate,thin=hyper$thin,gamma.mult=gamma.mult,
          phi=use.phi,pp.prior=pp.prior,has.sigma=hyper$has.sigma,
          has.hyper=hyper$has.hyper,force=force.data[i],
          random.phi=random.phi,mc.cores=cores) 
      } else {
        samples <- lapply(samples,run.chains,mci=i,p.migrate=p.migrate,
          thin=hyper$thin,gamma.mult=gamma.mult,phi=use.phi,pp.prior=pp.prior,
          has.sigma=hyper$has.sigma,has.hyper=hyper$has.hyper,
          force=force.data[i],random.phi=random.phi)
      }
     
      ps <- aperm(array(dim=c(dim(samples[[1]]$theta)[-3],length(samples)),
        data=unlist(lapply(samples,function(x){attr(x,"use")$theta}),use.names=FALSE),
        dimnames=list(NULL,samples[[1]]$p.names,NULL)),c(1,3,2))
      
      if ( i %% hyper$thin == 0 ) { # store samples
        store_i <- store_i + 1
        if (store_i %% report == 0) cat(store_i," ")
        hyper$h_summed_log_prior[store_i,] <- use.logprior
        hyper$h_log_likelihoods[store_i,] <- use.loglike
        hyper$phi[[1]][,,store_i] <- use.phi[[1]]
        hyper$phi[[2]][,,store_i] <- use.phi[[2]]
      }

  
    }
    cat("\n")
    if (cores>1 & os=="windows") { sfStop() }
    attr(samples,"hyper") <- hyper
  } else { # fixed effect sampling
    s.names <- names(samples)
    if (any(is.na(force.data))) force.data <- FALSE
    if ( cores>1 & length(samples)>1 & os=="windows") { 
      require(snowfall,quietly=TRUE)
      require(rlecuyer,quietly=TRUE)
      sfInit(parallel=TRUE, cpus=cores, type="SOCK")
      sfClusterSetupRNG()
      sfLibrary(msm)
      sfLibrary(rtdists)
      sfLibrary(statmod)
      sfLibrary(pracma)
      sfExportAll()
      samples <- sfLapply(samples,run.dmc,p.migrate=p.migrate,report=report,
                          force=force.data)
      sfStop() 
    } else if (cores>1) {
      require(parallel, quietly=TRUE)
      samples <- mclapply(samples,run.dmc,p.migrate=p.migrate,
        report=report,force=force.data,mc.cores=cores)
    } else {
      samples <- lapply(samples,run.dmc,p.migrate=p.migrate,report=report,
                        force=force.data)
    }
    names(samples) <- s.names
  }
  samples
}


### Posterior predictive

# n.post=100;probs=c(1:99)/100;bw="nrd0";save.simulation=FALSE;factors=NA
# save.subject.posts=FALSE; ignore.R2=FALSE; censor=c(NA,NA); gglist= FALSE
# probs.gglist =c(0.1, 0.5, 0.9);CI.gglist =c(0.025, 0.975); cores=1
# samples=hAhi1; gglist=TRUE

h.post.predict.dmc <- function(samples,n.post=100,probs=c(1:99)/100,bw="nrd0",
  save.simulation=FALSE,factors=NA,save.subject.posts=FALSE,cores=1,ignore.R2=FALSE, 
  gglist= FALSE, probs.gglist =c(0.1, 0.5, 0.9), CI.gglist =c(0.025, 0.975),
  censor=c(NA,NA))
  # apply post.predict to each subject
{
  os <- get.os()
  if ( cores>1 & length(samples)>1 & os=="windows") { 
    cat("No progress indication in multi-core\n")
    require(snowfall,quietly=TRUE)
    require(rlecuyer,quietly=TRUE)
    sfInit(parallel=TRUE, cpus=cores, type="SOCK")
    sfClusterSetupRNG()
    sfLibrary(msm)
    sfLibrary(rtdists)
    sfExportAll()
    out <- sfLapply(samples,post.predict.dmc,n.post=n.post,probs=probs,bw=bw,
                    factors=factors,save.simulation=save.simulation,gglist=TRUE,
                    save.simulation.as.attribute=TRUE,ignore.R2=ignore.R2,censor=censor)
      sfStop() 
  } else if (cores>1) {
    cat("No progress indication in multi-core\n")
    require(parallel, quietly=TRUE)
    out <- mclapply(samples,post.predict.dmc,n.post=n.post,probs=probs,bw=bw,
                    factors=factors,save.simulation=save.simulation,gglist=TRUE,
                    save.simulation.as.attribute=TRUE,ignore.R2=ignore.R2,censor=censor)
  } else {
    out <- lapply(samples,post.predict.dmc,n.post=n.post,probs=probs,bw=bw,
                  factors=factors,save.simulation=save.simulation,gglist=TRUE,
                  save.simulation.as.attribute=TRUE,ignore.R2=ignore.R2,censor=censor)
  }

  if ( !save.simulation ) { # Get averages
    sim <- do.call(rbind,lapply(out,function(x){attr(x,"sim")}))
    if ( (any(names(samples[[1]]$data)=="R2")) && !ignore.R2 ) for (i in 1:length(samples)) {
      levs <- outer(levels(samples[[i]]$data$R),sort(unique(samples[[i]]$data$R2)),"paste",sep="")
      if (attributes(attributes(samples[[i]]$data)$model)$type=="normDK") 
        levs[,2] <- rep("DK",dim(levs)[1])
      levs <- sort(unique(as.vector(levs)))  
      
      samples[[i]]$data$R <- 
        paste(as.character(samples[[i]]$data$R),as.character(samples[[i]]$data$R2),sep="") 
      if (attributes(attributes(samples[[i]]$data)$model)$type=="normDK") 
        samples[[i]]$data$R[samples[[i]]$data$R2=="2"] <- "DK"
      samples[[i]]$data$R <- factor(samples[[i]]$data$R,levels=levs)
    }
    dat <- do.call(rbind,lapply(samples,function(x){x$data}))
    facs <- names(attr(attributes(samples[[1]]$data)$model,"factors"))
    if (!is.null(factors)) {
      if (any(is.na(factors))) factors <- facs
      if (!all(factors %in% facs)) 
        stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
    }
    sim.dqp <- get.dqp(sim[,-1],factors,probs,n.post=1,bw=bw)
    dat.dqp <- get.dqp(sim=dat,factors,probs,bw=bw)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    av <- c(sim.dqp,dat.dqp)
    dpqs <- vector(mode="list",length=length(n.post))
    for (i in 1:n.post) {
      simi <- sim[sim$reps==i,-1]
      dpqs[[i]] <- get.dqp(simi,factors,probs,1)
    }
    attr(av,"dpqs") <- dpqs
    # Strip global sim attribute and dpqs for each participant
    out <- lapply(out,function(x){
      attr(x,"sim") <- NULL
      if (!save.subject.posts) attr(x,"dpqs") <- NULL
      x
    })
   if (gglist) attr(av, "gglist") <- get.fitgglist.dmc(sim,dat,factors=factors, noR=FALSE, 
     quantiles.to.get = probs.gglist, CI= CI.gglist)
   attr(out,"av") <- av
  }
  out
}



### Automatic sampling

# p.migrate=0;h.p.migrate=0; cores=1
# nmc=NA;report=10;cut=10;nbad=0;max.try=20;verbose=TRUE
# end.no.migrate=FALSE;h.gamma.mult=NA;gamma.mult=2.38;slaveOutfile=NULL
# 
# samples=hPES1BV

h.run.unstuck.dmc <- function(samples,nmc=NA,report=10,cores=1,
                            cut=10,nbad=0,max.try=20,verbose=TRUE,
                            p.migrate=0,h.p.migrate=0,end.no.migrate=FALSE,
                            h.gamma.mult=NA,gamma.mult=2.38,slaveOutfile=NULL)
  # Like run.unstuck but applied to list of subjects. If hyper present then
  # unsticks that as well as subjects.
{
  if ( !is.null(samples$theta) )
    stop("For a single subject use run.unstuck.dmc")
  if ( is.na(nmc) ) nmc <- samples[[1]]$nmc
  if ( any(is.na(samples[[1]]$theta[,,2])) ) {
    cat("Getting initial set of samples\n")
    samples <- h.run.dmc(samples=samples,report=report,cores=cores,
        gamma.mult=gamma.mult,h.gamma.mult=h.gamma.mult,
        p.migrate=p.migrate,h.p.migrate=h.p.migrate)
  }
  if ( any(names(attributes(samples))=="hyper") ) {
    try.num <- 0
    repeat {
      stucks <- lapply(samples,pick.stuck.dmc,cut=cut)
      ok <- c(hyper=length(pick.stuck.dmc(samples,hyper=TRUE,cut=cut)),
        unlist(lapply(stucks,function(x){length(x)})))
      
      if (verbose) {
        cat(paste("\nAfter try",try.num,"number of stuck chains:\n"))
        print(sort(ok[ok>0],decreasing=TRUE))
      }
      
      if (try.num > max.try | all(ok<=nbad)) break
      samples <- h.run.dmc(samples=h.samples.dmc(samples=samples,nmc=nmc),
        cores=cores,report=report,p.migrate=p.migrate,h.p.migrate=h.p.migrate,
        gamma.mult=gamma.mult)
      try.num <- try.num + 1
    }
    if (end.no.migrate) h.run.dmc(samples=h.samples.dmc(samples=samples,nmc=nmc),
        cores=cores,report=report,gamma.mult=gamma.mult)
  } else {
    os <- get.os()
    if ( cores>1 & os=="windows") {
      require(snowfall,quietly=TRUE)
      require(rlecuyer,quietly=TRUE)
      sfInit(parallel=TRUE, cpus=cores, type="SOCK",slaveOutfile=slaveOutfile)
      sfClusterSetupRNG()
      sfLibrary(msm)
      sfLibrary(rtdists)
      sfLibrary(statmod)
      sfLibrary(pracma)
      sfExportAll()
      samples <- sfLapply(samples,run.unstuck.dmc,nmc=nmc,report=report,cores=1,
                          cut=cut,nbad=nbad,max.try=max.try,p.migrate=p.migrate,
                          gamma.mult=gamma.mult,verbose=verbose)
      if (cores>1) sfStop()
    } else if (cores>1) {
      require(parallel, quietly=TRUE)
      samples <- mclapply(samples,run.unstuck.dmc,nmc=nmc,
        report=report,cores=1,cut=cut,nbad=nbad,max.try=max.try,
        p.migrate=p.migrate,gamma.mult=gamma.mult,verbose=verbose,
        mc.cores=cores)
    } else {
      samples <- lapply(samples,run.unstuck.dmc,nmc=nmc,report=report,
                          cores=1,gamma.mult=gamma.mult,verbose=verbose,
                          cut=cut,nbad=nbad,max.try=max.try,p.migrate=p.migrate)
    }
  }
  samples
}

# report=1;cores=1;gamma.mult=2.38; random.theta=TRUE
# h.gamma.mult=NA;random.phi=TRUE;cut=1.1;max.try=100;minN=NA;meanN=NA
# slaveOutfile=NULL;digits=2;transform=TRUE;autoburnin=FALSE;split=TRUE;verbose=FALSE
# finalrun=TRUE;finalI=NA;finalminN=NA,finalmeanN=NA; addtofinal=FALSE; removefromfinal=NA
# 
# samples=hsamples0;nmc=4;cores=18;report=1;verbose=TRUE


h.run.converge.dmc <- function(samples,nmc=NA,report=10,cores=1,gamma.mult=2.38,
  h.gamma.mult=NA,random.phi=TRUE,random.theta=TRUE,cut=1.1,max.try=20,minN=NA,meanN=NA,
  slaveOutfile=NULL,digits=2,transform=TRUE,autoburnin=FALSE,split=TRUE,
  save="",verbose=TRUE,thorough=FALSE,
  finalrun=FALSE,finalI=NA,finalminN=NA,finalmeanN=NA,addtofinal=FALSE,removefromfinal=NA)
  # Like run.unstuck but applied to list of subjects. 
  # If !thorough and hyper then removes based on only hyper gelman.diag and 
  # effective sample size at hyper level only but criterion on gelman.diag 
  # applied to hyper AND subjects. 
  # If thorough and hyper will check every individual parameter Rhat for both 
  # hyper and subjects, exiting when maximum is less than cut. Effective sample
  # size same as not thorough.
  # If finalrun=TRUE gets final fresh set of samples fulfilling finalI/minN/meanN
  # If addtofinal adds on to previous, and can also removefromfinal=1:N before doing so
  # h.save will save each random effects (heirarchical) try to an .RData file
{
  
  get.gd <- function(samples,thorough=TRUE,autoburnin=FALSE,transform=TRUE,split=TRUE,verbose=FALSE) {
    if (thorough)
      c(gelman.diag.dmc(samples,hyper=TRUE)$psrf[,1],
        unlist(lapply(gelman.diag.dmc(samples),function(x){x$psrf[,1]}))) else
      h.gelman.diag.dmc(samples,autoburnin=autoburnin,
        transform=transform,split=split,verbose=FALSE)
  }

  if ( !is.null(samples$theta) )
    stop("For a single subject use run.converge.dmc")
  if ( is.na(nmc) ) nmc <- samples[[1]]$nmc
  if ( any(is.na(samples[[1]]$theta)) )
    samples <- h.run.dmc(samples=samples,report=report,cores=cores,
        gamma.mult=gamma.mult,h.gamma.mult=h.gamma.mult,
      random.phi=random.phi,random.theta=random.theta)
  if ( any(names(attributes(samples))=="hyper") ) {
    try.num <- 0
    if (!is.na(minN) & !is.na(meanN)) {
      warning("Both minN and meanN specified, using minN")
      meanN <- NA 
    }
    if ( !is.na(minN) ) nfun <- "min"
    if (!is.na(meanN)) {
      nfun <- "mean"
      minN <- meanN
    }

    if ( finalrun & is.na(finalI) ) {
      if ( is.na(finalminN) & is.na(finalmeanN) )
        stop("Must specify a finalI or a finalminN or a finalmeanN if finalrun=TRUE")
      if ( !is.na(finalminN) & !is.na(finalmeanN) ) {
        warning("Both finalminN and finalmeanN specified, using finalminN")
        finalmeanN <- NA 
      }
      if ( !is.na(finalminN) ) finalfun <- "min"
      if (!is.na(finalmeanN)) {
        finalfun <- "mean"
        finalminN <- finalmeanN
      }
    }    

    gd <- get.gd(samples,thorough=thorough,autoburnin=autoburnin,
      transform=transform,split=split,verbose=FALSE)
    if ( !all(gd < cut) | is.na(minN) ) effectiveN <- NA else
      effectiveN <- do.call(nfun,list(effectiveSize.dmc(samples,hyper=TRUE)))
    if (verbose) {
      if (thorough)
        cat(paste("Iterations = ",dim(attr(samples,"hyper")$phi[[1]])[3],
          ", Effective N = ",effectiveN,", Maximum psrf = ",
          round(max(gd),digits=digits),sep="")) else 
        {
          cat(paste("Iterations = ",dim(attr(samples,"hyper")$phi[[1]])[3],
            ", Effective N = ",effectiveN,", Hyper mpsrf = ",
            round(gd["hyper"],digits=digits),
            "\nSubject mpsrf achieved (sorted):\n",sep=""))
          print(round(sort(gd[names(gd)!="hyper"],decreasing=TRUE),digits=digits))
        }
    }
    if ( !(all(gd < cut)) | (!is.na(effectiveN) && effectiveN < minN) ) {
      repeat {
        try.num <- try.num + 1
        samples <- h.run.dmc(samples=h.samples.dmc(samples=samples,nmc=nmc,add=TRUE),
          cores=cores,report=report,gamma.mult=gamma.mult,
          h.gamma.mult=h.gamma.mult,random.phi=random.phi,random.theta=random.theta)
        gd <- get.gd(samples,thorough=thorough,autoburnin=autoburnin,
          transform=transform,split=split,verbose=FALSE)
        shorter <- h.samples.dmc(samples=samples,remove=1:nmc,nmc=0,add=TRUE) 
        gd.short <- get.gd(shorter,thorough=thorough,autoburnin=autoburnin,
          transform=transform,split=split,verbose=FALSE)
        if (thorough) shorten <- max(gd.short) < max(gd) else
          shorten <- gd.short["hyper"] < gd["hyper"]
        if ( shorten ) {
            samples <- shorter
            gd <- gd.short
            if (verbose) cat(paste("Discarding initial",nmc,"samples.\n"))
        }
        if (verbose) {
          if (thorough)
          cat(paste("Iterations = ",dim(attr(samples,"hyper")$phi[[1]])[3],
            ", Effective N = ",effectiveN,", Maximum psrf = ",
            round(max(gd),digits=digits),sep="")) else 
          {
            cat(paste("Iterations = ",dim(attr(samples,"hyper")$phi[[1]])[3],
              ", Effective N = ",effectiveN,", Hyper mpsrf = ",
              round(gd["hyper"],digits=digits),
              "\nSubject mpsrf achieved (sorted):\n",sep=""))
            print(round(sort(gd[names(gd)!="hyper"],decreasing=TRUE),digits=digits))
          }
        }
        if ( try.num >= max.try ) break
        if ( all(gd <= cut) ) {
          if ( is.na(minN) ) break else {
            effectiveN <- do.call(nfun,list(effectiveSize.dmc(samples,hyper=TRUE)))
            if ( effectiveN > minN ) break
          } 
        }
        if (save != "") save(samples,file=paste(save,"RData",sep="."))
      }
    }
    if ( addtofinal && (finalrun & !any(is.na(removefromfinal))) && 
         (max(removefromfinal)>dim(samples[[1]]$theta)[3]) ) {
      finalrun <- FALSE
      warning("Final run aborted as cutfromfinal to large")
    }
    if (finalrun) {
      try.num <- 0
      if ( addtofinal & !any(is.na(removefromfinal)) )
        samples <- h.samples.dmc(samples=samples,nmc=0,add=TRUE,remove=removefromfinal)
      if (verbose) cat("\n\nDoing final run\n")
      if ( !is.na(finalI) ) {
        if ( addtofinal ) {
         finalI <- finalI-dim(samples[[1]]$theta)[3]
         if (finalI > 0) {
           samples <- h.samples.dmc(samples=samples,nmc=finalI,add=TRUE)
           samples <- h.run.dmc(samples,cores=cores,report=report,gamma.mult=gamma.mult,
            h.gamma.mult=h.gamma.mult,random.phi=random.phi,random.theta=random.theta) 
         } else {
           warning("Already more than finalI available")
         }
        } else samples <- h.run.dmc(samples=h.samples.dmc(samples=samples,nmc=finalI),
          cores=cores,report=report,gamma.mult=gamma.mult,
          h.gamma.mult=h.gamma.mult,random.phi=random.phi,random.theta=random.theta)
      } else {
        samples <- h.run.dmc(samples=h.samples.dmc(samples=samples,nmc=nmc,add=addtofinal),
          cores=cores,report=report,gamma.mult=gamma.mult,
          h.gamma.mult=h.gamma.mult,random.phi=random.phi,random.theta=random.theta)
        effectiveN <- do.call(finalfun,list(effectiveSize.dmc(samples,hyper=TRUE)))
        if (verbose) cat(paste("Iterations = ",
              dim(attr(samples,"hyper")$phi[[1]])[3],", Effective N = ",effectiveN,"\n\n",sep="")) 
        if ( effectiveN < finalminN ) {
          repeat {
            try.num <- try.num + 1
            samples <- h.run.dmc(samples=h.samples.dmc(samples=samples,nmc=nmc,add=TRUE),
              cores=cores,report=report,gamma.mult=gamma.mult,
              h.gamma.mult=h.gamma.mult,random.phi=random.phi,random.theta=random.theta)
            if ( try.num >= max.try ) break
            effectiveN <- do.call(finalfun,list(effectiveSize.dmc(samples,hyper=TRUE)))
            if (verbose) cat(paste("Iterations = ",
              dim(attr(samples,"hyper")$phi[[1]])[3],", Effective N = ",effectiveN,"\n\n",sep="")) 
            if ( effectiveN > finalminN ) break
            if (save != "") save(samples,file=paste(save,"RData",sep="."))
          }
        }
      }
    }
  } else {
    os <- get.os()
    if ( cores>1 & os=="windows") {
      require(snowfall,quietly=TRUE)
      require(rlecuyer,quietly=TRUE)
      sfInit(parallel=TRUE, cpus=cores, type="SOCK",slaveOutfile=slaveOutfile)
      sfClusterSetupRNG()
      sfLibrary(msm)
      sfLibrary(rtdists)
      sfLibrary(statmod)
      sfLibrary(pracma)
      sfLibrary(coda)
      sfExportAll()
      samples <- sfLapply(samples,run.converge.dmc,report=report,cores=1,
        gamma.mult=gamma.mult,cut=cut,max.try=max.try,minN=minN,meanN=meanN,
        nmc=nmc,transform=transform,autoburnin=autoburnin,split=split,
        verbose=verbose)
      sfStop()
    } else if (cores>1) {
      require(parallel, quietly=TRUE)
      samples <- mclapply(samples,run.converge.dmc,report=report,
        cores=1,gamma.mult=gamma.mult,cut=cut,max.try=max.try,minN=minN,
        meanN=meanN,nmc=nmc,transform=transform,autoburnin=autoburnin,
        split=split,verbose=verbose,mc.cores=cores)
    } else {
      samples <- lapply(samples,run.converge.dmc,report=report,cores=cores,
        gamma.mult=gamma.mult,cut=cut,max.try=max.try,minN=minN,meanN=meanN,
        nmc=nmc,transform=transform,autoburnin=autoburnin,split=split,verbose=verbose)
    }
  }
  samples
}

  



# cores=1;report=10;p.migrate=.05;max.try=100
# cut.unstuck=10;cut.flat.location=1/4;cut.flat.scale=Inf
# cut.converge=1.1;split=TRUE;minN=NA;meanN=NA
# verbose=FALSE;gamma.mult=2.38;subjects.to.cores=TRUE;slaveOutfile=NULL
# 
# cores=3;report=1;max.try=1

# cores=1;report=10;p.migrate=.05;max.try=100
# cut.unstuck=10;cut.flat.location=1/4;cut.flat.scale=Inf
# cut.converge=1.1;split=TRUE;minN=NA;meanN=NA
# verbose=FALSE;gamma.mult=2.38;subjects.to.cores=TRUE;slaveOutfile=NULL

# cores=1;report=10;p.migrate=.05;max.try=100
# cut.unstuck=10;cut.flat.location=1/4;cut.flat.scale=Inf
# cut.converge=1.1;split=TRUE;minN=NA;meanN=NA
# verbose=FALSE;gamma.mult=2.38;subjects.to.cores=TRUE;slaveOutfile=NULL
# 
# cores=3;report=1;max.try=1

# cores=1;report=10;p.migrate=.05;max.try=100
# cut.unstuck=10;cut.flat.location=1/4;cut.flat.scale=Inf
# cut.converge=1.1;split=TRUE;minN=NA;meanN=NA
# verbose=FALSE;gamma.mult=2.38;subjects.to.cores=TRUE;slaveOutfile=NULL

h.RUN.dmc <- function(hsamples,
                      cores=1,
                      report=10,
                      p.migrate=.05,
                      h.p.migrate=.05,
                      max.try=100,
                      cut.unstuck=10,
                      cut.flat.location=.5,
                      cut.flat.scale=.5,
                      cut.converge=1.1,
                      split=TRUE,
                      minN=NA,
                      meanN=NA,
                      use.effectiveSize=TRUE,
                      n.add=NA,
                      force=FALSE,
                      thorough=TRUE,
                      verbose=FALSE,
                      gamma.mult=2.38,
                      h.gamma.mult=NA,
                      subjects.to.cores=TRUE,
                      slaveOutfile=NULL)
  
{
  
  if ( !is.null(hsamples$theta) ) {
    
    ### single participant ###
    
    stop("For a single subject use RUN.dmc")
    
  } else if (any(names(attributes(hsamples))=="hyper")) {
    
    ### multiple participants, truly hierarchical ###
    
    # functions for checking
    h.StuckTests <- function(samples,
                             thorough=TRUE,
                             verbose=FALSE,
                             cut) {
      
      if (verbose) cat("Stuck chains check\n") 
      
      if (thorough) {
        # check both hyper & participant-level chains
        stucks <- lapply(samples,pick.stuck.dmc,cut=cut)
        ok <- c(hyper=length(pick.stuck.dmc(samples,hyper=TRUE,cut=cut)),
                unlist(lapply(stucks,length)))
      } else {
        # check only hyper chains
        ok <- c(hyper=length(pick.stuck.dmc(samples,hyper=TRUE,cut=cut)))
      }
      
      fail <- any( ok > 0 )
      if (verbose) {
        if (!fail) cat(": OK\n") else
          cat(paste(":",sum(ok),"\n"))
      }
      fail
    }
    
    h.FlatTests <- function(samples,
                            thorough=TRUE,
                            p1=1/3,
                            p2=1/3,
                            cut.location=0.25,
                            cut.scale=Inf,
                            verbose=FALSE) {
      
      gmcmc <- function(samples, thorough) {
        
        h <- attr(samples, "hyper")
        phi1 <- matrix(aperm(h$phi[[1]],c(1,3,2)),
                       ncol=dim(h$phi[[1]])[2],
                       dimnames=list(NULL, paste0(dimnames(h$phi[[1]])[[2]],
                                                  ".h1")))
        phi2 <- matrix(aperm(h$phi[[2]],c(1,3,2)),
                       ncol=dim(h$phi[[2]])[2],
                       dimnames=list(NULL, paste0(dimnames(h$phi[[2]])[[2]],
                                                  ".h2")))
        m <- cbind(phi1, phi2)
        
        if (thorough) {
          
          # add participant-level chains
          
          mindlist <- lapply(seq_along(samples), function(i) {
            m_i <- matrix(aperm(samples[[i]]$theta,c(1,3,2)),
                          ncol=dim(samples[[i]]$theta)[2],
                          dimnames=list(NULL,
                          paste0(i, ".", dimnames(samples[[i]]$theta)[[2]])))
          })
          
          mind <- do.call(cbind, mindlist)
          m <- cbind(m, mind)
          
        }
        
        mcmc(m)
        
      }
      
      mat <- gmcmc(samples, thorough = thorough)
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
        # Change in IQR relative to overall IQR 
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
    
    h.MixTests <- function(samples,
                           thorough=TRUE,
                           verbose=FALSE,
                           cut,
                           split=TRUE) {
      
      tmp <- gelman.diag.dmc(samples, hyper=TRUE, split=split)
      names(tmp$mpsrf) <- "hyper_mpsrf"
      gds <- c(tmp$mpsrf,tmp$psrf[,1])
      
      if (thorough) {
        # add participant-level
        tmp2 <- gelman.diag.dmc(samples, split=split)
        tmp3 <- lapply(seq_along(tmp2), function(i) {
          psrf_i <- tmp2[[i]]$psrf[,1]
          names(psrf_i) <- paste0(i, ".", names(psrf_i))
          mpsrf_i <- tmp2[[i]]$mpsrf
          names(mpsrf_i) <- paste0(i, ".mpsrf")
          return(c(mpsrf_i, psrf_i))
        })
        gds <- c(gds, unlist(tmp3))
      }
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
    
    h.LengthTests <- function(samples,
                              minN,
                              nfun,
                              thorough=TRUE,
                              verbose=FALSE) {
      
      n <- do.call(nfun,list(h.get.size(samples, thorough = thorough)))
      
      fail <- n < minN
      if (verbose) { 
        cat("Length check")
        if (!fail) cat(": OK\n") else cat(paste(":",n,"\n"))  
      }
      fail
    }
    
    
    if (!verbose) report <- 1e8
    
    if (use.effectiveSize) {
      h.get.size <- function(hsamples, thorough) {
        neff <- effectiveSize.dmc(hsamples, hyper = TRUE)
        if (thorough) {
          neff <- c(neff, unlist(effectiveSize.dmc(hsamples)))
        }
        return(neff)
      }
    } else {
      h.get.size <- function(x, thorough){prod(dim(x[[1]]$theta)[-2])}
    }
    
    if ( !is.na(minN) & !is.na(meanN) ) {
      warning("Both minN and meanN specified, using minN")
      meanN <- NA 
    }
    if ( !is.na(minN) ) nfun <- "min"
    if ( !is.na(meanN) ) {
      nfun <- "mean"
      minN <- meanN
    }
    
    if ( is.na(n.add) ) n <-  ceiling(hsamples[[1]]$nmc/3) else {
      if ( n.add >= floor(hsamples[[1]]$nmc/2) )
        stop(paste("n.flat.test to large, must be less than",
                   floor(hsamples[[1]]$nmc/2)))
      n <- n.add
    }
    
    if ( any(!is.finite(hsamples[[1]]$log_likelihoods[hsamples[[1]]$nmc,])) ) 
    { # New samples
      do.migrate=1
      if (verbose) {
        cat("\nGetting initial set of samples\n")
      }
      hsamples <- h.run.dmc(samples=hsamples,report=report,cores=cores,
                            gamma.mult=gamma.mult,h.gamma.mult=h.gamma.mult,
                            p.migrate=p.migrate,h.p.migrate=h.p.migrate)
    } else do.migrate=0 # Samples already good, just want more.
    
    if ( !force && !is.null(attr(hsamples,"auto")) && 
         !is.na(attr(hsamples,"auto")) && attr(hsamples,"auto")!="GRID FAIL" )
      return(hsamples) # If already sucessfully run then dont run again unless force=TRUE
    n.try <- 0
    repeat {
      
      # DO TESTS
      
      # Stuck check
      if ( h.StuckTests(hsamples,
                        thorough=thorough,
                        cut=cut.unstuck,
                        verbose=verbose) ) 
      { # Start again 
        do.migrate <- 1
        get.new <- test.again <- TRUE
      } else get.new <- test.again <- FALSE
      
      # Stationarity check
      if ( !get.new && any(h.FlatTests(hsamples,
                                       verbose=verbose,
                                       thorough=thorough,
                                       cut.location=cut.flat.location,
                                       cut.scale=cut.flat.scale)) ) 
      { # remove first 1/3, add new 1/3
        if (verbose) cat("Removing initial 1/3 of samples\n")
        nshift <- ceiling(hsamples[[1]]$nmc/3)
        hsamples <- h.samples.dmc(samples=hsamples,remove=1:nshift,
                                  nmc=0,add=TRUE)  
        hsamples <- h.samples.dmc(samples=hsamples,add=TRUE,nmc=n)
        test.again <- TRUE  
        get.new <- FALSE
      } 
      
      # Mixing check
      if ( (!get.new & !test.again) && 
           h.MixTests(hsamples,
                      thorough=thorough,
                      verbose=verbose,
                      cut=cut.converge,
                      split=split) ) { 
        { # Make longer
          hsamples <- h.samples.dmc(samples=hsamples,nmc=n,add=TRUE)
          test.again <- TRUE
          get.new <- FALSE
        }
      }
      
      # Length check
      if ( (!get.new & !test.again) &&
           ( !is.na(minN) && h.LengthTests(hsamples,
                                           thorough=thorough,
                                           minN=minN,
                                           nfun=nfun,
                                           verbose=verbose)) ) { # Not long enough
        hsamples <- h.samples.dmc(samples=hsamples,nmc=n,add=TRUE)
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
        hsamples <- h.run.dmc(h.samples.dmc(samples=hsamples,nmc=3*n,
                                            thin=hsamples[[1]]$thin),
                              report=report,cores=cores,
                              gamma.mult=gamma.mult,
                              h.gamma.mult=h.gamma.mult,
                              p.migrate=p.migrate,
                              h.p.migrate=h.p.migrate)
      } else {         # Test failed, update
        if (verbose) cat("Adding extra samples\n")
        hsamples <- h.run.dmc(hsamples,
                              report=report,cores=cores,
                              gamma.mult=gamma.mult,
                              h.gamma.mult=h.gamma.mult,
                              p.migrate=p.migrate*do.migrate)
      }
      
      # Give up?
      n.try <- n.try + 1
      if ( (n.try > max.try) ) {
        outcome <- "FAIL"
        if (verbose) cat(paste("\nOUTCOME:",outcome,"AFTER TRY",n.try,"\n"))
        break 
      } else if (verbose) cat(paste("COMPLETED TRY",n.try,"\n\n"))
    }
    if (outcome=="FAIL") attr(hsamples,"auto") <- NA else
      attr(hsamples,"auto") <- n.try
    
  } else {
    
    ### multiple participants, not truly hierarchical ###
    
    if ( cores == 1 | !subjects.to.cores ) 
      hsamples <- lapply(hsamples,RUN.dmc,cores=cores,report=report,
                         p.migrate=p.migrate,max.try=max.try,cut.unstuck=cut.unstuck,
                         cut.flat.location=cut.flat.location,cut.flat.scale=cut.flat.scale,
                         cut.converge=cut.converge,split=split,
                         minN=minN,meanN=meanN,use.effectiveSize = use.effectiveSize,
                         verbose=verbose,gamma.mult=gamma.mult) else 
    {
      os <- get.os()
      if ( os=="windows" ) {
        require(snowfall,quietly=TRUE)
        require(rlecuyer,quietly=TRUE)
        sfInit(parallel=TRUE, cpus=cores, type="SOCK",slaveOutfile=slaveOutfile)
        sfClusterSetupRNG()
        sfLibrary(msm)
        sfLibrary(rtdists)
        sfLibrary(statmod)
        sfLibrary(pracma)
        sfLibrary(coda)
        sfExportAll()
        if (subjects.to.cores==TRUE) hsamples <- 
          sfLapply(hsamples,RUN.dmc,cores=1,report=report,
                   p.migrate=p.migrate,max.try=max.try,cut.unstuck=cut.unstuck,
                   cut.flat.location=cut.flat.location,cut.flat.scale=cut.flat.scale,
                   cut.converge=cut.converge,split=split,
                   minN=minN,meanN=meanN,use.effectiveSize = use.effectiveSize,
                   verbose=verbose,gamma.mult=gamma.mult) 
        sfStop()
      } else {
        require(parallel, quietly=TRUE)
        hsamples <- mclapply(hsamples,RUN.dmc,cores=1,report=report,
                             p.migrate=p.migrate,max.try=max.try,cut.unstuck=cut.unstuck,
                             cut.flat.location=cut.flat.location,cut.flat.scale=cut.flat.scale,
                             cut.converge=cut.converge,split=split,
                             minN=minN,meanN=meanN,use.effectiveSize = use.effectiveSize,
                             verbose=verbose,gamma.mult=gamma.mult,mc.cores=cores)
      } 
    }
    
    cat(paste("Number of trys (NA implies fails after max.try =",max.try,")\n"))
    print(unlist(lapply(hsamples,function(x){attr(x,"auto")})))
    
  }
  hsamples
}

# h.RUN.dmc <- function(hsamples,cores=1,report=10,p.migrate=.05,max.try=100,
#   cut.unstuck=10,
#   cut.flat.location=.5,cut.flat.scale=.5,
#   cut.converge=1.1,split=TRUE,
#   minN=NA,meanN=NA,use.effectiveSize = TRUE,
#   verbose=FALSE,gamma.mult=2.38,
#   subjects.to.cores=TRUE,slaveOutfile=NULL)
#   
# {
#   if ( cores == 1 | !subjects.to.cores ) 
#     hsamples <- lapply(hsamples,RUN.dmc,cores=cores,report=report,
#       p.migrate=p.migrate,max.try=max.try,cut.unstuck=cut.unstuck,
#       cut.flat.location=cut.flat.location,cut.flat.scale=cut.flat.scale,
#       cut.converge=cut.converge,split=split,
#       minN=minN,meanN=meanN,use.effectiveSize = use.effectiveSize,
#       verbose=verbose,gamma.mult=gamma.mult) else 
#   {
#     os <- get.os()
#     if ( os=="windows" ) {
#       require(snowfall,quietly=TRUE)
#       require(rlecuyer,quietly=TRUE)
#       sfInit(parallel=TRUE, cpus=cores, type="SOCK",slaveOutfile=slaveOutfile)
#       sfClusterSetupRNG()
#       sfLibrary(msm)
#       sfLibrary(rtdists)
#       sfLibrary(statmod)
#       sfLibrary(pracma)
#       sfLibrary(coda)
#       sfExportAll()
#       if (subjects.to.cores==TRUE) hsamples <- 
#         sfLapply(hsamples,RUN.dmc,cores=1,report=report,
#           p.migrate=p.migrate,max.try=max.try,cut.unstuck=cut.unstuck,
#           cut.flat.location=cut.flat.location,cut.flat.scale=cut.flat.scale,
#           cut.converge=cut.converge,split=split,
#           minN=minN,meanN=meanN,use.effectiveSize = use.effectiveSize,
#           verbose=verbose,gamma.mult=gamma.mult) 
#       sfStop()
#     } else {
#       require(parallel, quietly=TRUE)
#       hsamples <- mclapply(hsamples,RUN.dmc,cores=1,report=report,
#         p.migrate=p.migrate,max.try=max.try,cut.unstuck=cut.unstuck,
#         cut.flat.location=cut.flat.location,cut.flat.scale=cut.flat.scale,
#         cut.converge=cut.converge,split=split,
#         minN=minN,meanN=meanN,use.effectiveSize = use.effectiveSize,
#         verbose=verbose,gamma.mult=gamma.mult,mc.cores=cores)
#     } 
#   }
#   cat(paste("Number of trys (NA implies fails after max.try =",max.try,")\n"))
#   print(unlist(lapply(hsamples,function(x){attr(x,"auto")})))
#   hsamples
# }

####### Use grid to fit each subject seperately


setup.grid <- function(hsamples,dname,froot="s") {
  dir.create(dname)
  for (i in 1:length(hsamples)) {
    onam <- paste(froot,i,sep=".")
    fnam <- paste(dname,"/",onam,".RData",sep="")
    assign(onam,hsamples[[i]])
    save(list=onam,file=fnam)
  }
}


# collect.grid <- function(dname,froot="s",verbose=TRUE,hsamples=NA) {
#   if (any(is.na(hsamples))) load(paste(dname,"/","hsamples.RData",sep=""))
#   for (i in 1:length(hsamples)) {
#     onam <- paste(froot,i,sep=".")
#     fnam <- paste(dname,"/results.",onam,".RData",sep="")
#     tmp <- try(load(fnam),silent=TRUE)
#     if (class(tmp)=="try-error") bad <- TRUE else {
#       hsamples[[i]] <- get(tmp)
#       bad <- FALSE
#     }
#     if (verbose) {
#       if (bad) cat("Read fail\n") else 
#         cat(paste(onam,":",attr(hsamples[[i]],"auto"),"\n")) 
#     } else {
#       if (bad) cat("F") else cat(".")
#     }
#   }
#   cat("\n")
#   auto <- unlist(lapply(hsamples,function(x){
#     y <- attr(x,"auto"); if (is.null(y)) y<-NA; y
#   }))
#   attr(hsamples,"auto") <- auto
#   hsamples
# }

collect.grid <- function(dname,froot="s",verbose=TRUE,hsamples=NA) {
  if (any(is.na(hsamples))) load(paste(dname,"/","hsamples.RData",sep=""))
  for (i in 1:length(hsamples)) {
    onam <- paste(froot,i,sep=".")
    fnam <- paste(dname,"/results.",onam,".RData",sep="")
    tmp <- try(load(fnam),silent=TRUE)
    if (class(tmp)=="try-error") bad <- TRUE else {
      hsamples[[i]] <- get(tmp)
      bad <- FALSE
    }
    if (bad) attr(hsamples[[i]],"auto") <- "GRID FAIL"  
    if (verbose) {
      if (bad) cat("Read fail\n") else 
        cat(paste(onam,":",attr(hsamples[[i]],"auto"),"\n")) 
    } else {
      if (bad) cat("F") else cat(".")
    }
  }
  cat("\n")
  auto <- unlist(lapply(hsamples,function(x){
    y <- attr(x,"auto"); if (is.null(y)) y<-NA; y
  }))
  attr(hsamples,"auto") <- auto
  hsamples
}



cleanup.grid <- function(dname,froot="s") {
  load(paste(dname,"/","hsamples.RData",sep=""))
  file.remove(paste(dname,"/","hsamples.RData",sep=""))
  for (i in 1:length(hsamples)) {
    onam <- paste(froot,i,sep=".")
    fnam <- paste(dname,"/",onam,".RData",sep="")
    file.remove(fnam)
    fnam <- paste(dname,"/",onam,"results.RData",sep="")
    file.remove(fnam)
  }
}


# model.file="exgSSprobit.R"
# fname="pat5"; user="ajh296"
# froot="s";verbose=TRUE
# model.file="lnrSS.R"
# GB=2 # GB per job
# nservers=1; ncpus=1
# wall.hours=100 # Wall time in hours
# n.add=40
# max.try=20
# minN=2000;meanN=NA
# cut.unstuck=10
# cut.flat.location=.5;cut.flat.scale=.5
# report=10;p.migrate=.05
# cut.converge=1.1;split=TRUE;gamma.mult=2.38
# farjump=NA; force=FALSE; RUN=TRUE
# module_load="source /etc/profile.d/modules.sh\nmodule load R/3.3.2"



run.grid.dmc <- function(fname,model.dir,model.file,user,n.add, # usually 1/3 of initial size
                         nservers=1, ncpus=1, # Number os servers and cpus/server
                         GB=2, # GB per job 
                         wall.hours=100, # Wall time in hours
                         # required by grid to access R
                         module_load="source /etc/profile.d/modules.sh\nmodule load R/3.3.2",
                         froot="s",
                         verbose=TRUE,
                         RUN=TRUE,
                         # h.RUN.dmc parameters
                         force.RUN=FALSE, # dont run if file has an auto attribute 
                         max.try=20,
                         minN=NA,meanN=NA,
                         cut.unstuck=10,
                         cut.flat.location=.5,
                         cut.flat.scale=.5,
                         cut.converge=1.1,split=TRUE,
                         # run.dmc parameters
                         farjump=NA, 
                         force=FALSE,
                         # common parameters
                         p.migrate=.05,
                         report=10,
                         gamma.mult=2.38
)
{
  oname <- load(paste(fname,"RData",sep="."))
  hsamples <- get(oname)
  setup.grid(hsamples,fname,froot)
  # Write sh and R run files
  cat(paste("#!/bin/bash\n#PBS -l select=",nservers,":ncpus=",ncpus,":mem=",GB,"gb\n",
  "#PBS -l walltime=",wall.hours,":00:00\n#PBS -j oe\n#PBS -N ",fname,
  "\n",module_load,"\n","cd $PBS_O_WORKDIR\nRscript run-",fname,".R\nexit 0",sep=""),
  file=paste("run-",fname,".sh",sep=""))
  if (RUN) run.call <- paste("RUN.dmc(get(tmp), cores =",ncpus,", max.try =",max.try,
         ", minN =",minN,", meanN =",meanN,
         ", cut.flat.location =",cut.flat.location,
         ", cut.flat.scale =",cut.flat.scale,", n.add =",n.add,
         ", cut.unstuck =",cut.unstuck,", p.migrate =",p.migrate,
         ", cut.converge =",cut.converge,", split =",split,
         ", gamma.mult =",gamma.mult,", force =",force.RUN,
         ", verbose =",verbose,", report =",report,")") else
    run.call <- paste("run.dmc(get(tmp),report=",report,",cores=",ncpus,
                      ",p.migrate=",p.migrate,",gamma.mult=",gamma.mult,
                      ",farjump=",farjump,",force=",force,")")
  cat(paste("fname <- \"",fname,"\"; froot <- \"",froot,"\"\n",
    "source (\"dmc/dmc.R\")\nload_model(\"",model.dir,"\",\"",model.file,"\")\n",
    "num <- Sys.getenv(\"PBS_ARRAY_INDEX\")\n",
    "onam <- paste(froot,num,sep=\".\")\n",
    "fnam <- paste(fname,\"/\",onam,\".RData\",sep=\"\")\n",
    "tmp=load(fnam)\ntmp <-",run.call,"\n",
    "assign(onam,tmp)\nfnam <- paste(fname,\"/results.\",onam,\".RData\",sep=\"\")\n",
    "save(list=onam,file=fnam)",sep=""),file=paste("run-",fname,".R",sep=""))
  jnum <- system(paste("qsub -J 1-",length(hsamples)," run-",fname,".sh",sep=""),intern=TRUE)
  jnum <- strsplit(jnum,"[].rcgbcm",fixed=TRUE)[[1]]
  
  # Wait for run to finish 
  was.done <- ""
  repeat {
    done <- dir(pattern=jnum)
    if ( length(done)>length(was.done) ) {
      cat(paste(length(was.done)+c(1:sum(!(done %in% was.done))),
                done[!(done %in% was.done)],"\n"))
      was.done <- done
    }
    qstat <- system(paste("qstat -u",user),intern=TRUE)
    if ( length(qstat[grep(jnum,qstat)])==0 ) break
    Sys.sleep(1)
  }

  cat("Harvesting files: ")
  for (i in 1:length(hsamples)) {
    onam <- paste(froot,i,sep=".")
    fnam <- paste(fname,"/results.",onam,".RData",sep="")
    tmp <- try(load(fnam),silent=TRUE)
    if (class(tmp)=="try-error") bad <- TRUE else {
      hsamples[[i]] <- get(tmp)
      bad <- FALSE
    }
    if (bad) attr(hsamples[[i]],"auto") <- "GRID FAIL"  
    if (verbose) {
      if (bad) cat("Read fail\n") else 
        cat(paste(onam,":",attr(hsamples[[i]],"auto"),"\n")) 
    } else {
      if (bad) cat("F") else cat(".")
    }
  }
  cat("\n")
  auto <- unlist(lapply(hsamples,function(x){attr(x,"auto")}))
  attr(hsamples,"auto") <- auto
  assign(oname,hsamples)
  save(list=oname,file=paste(fname,"RData",sep="."))

  
  # Clean up
  unlink(fname,recursive=TRUE)
  file.remove(paste("run-",fname,".R",sep=""))
  file.remove(paste("run-",fname,".sh",sep=""))
  unlink(paste(fname,".o",jnum,".*",sep=""))
}





