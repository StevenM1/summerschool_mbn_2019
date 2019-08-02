# System functions for the DMC (Dynamic Models of Choice)
#    Functions specific to extract and analyze samples
#    Usually user does not need to edit

### Utilities ----

get.os <- function()
{
  sysinf <- Sys.info()
  ostype <- .Platform$OS.type

  ## Probe using Sys.info: Windows, Linux or Darwin
  if (!is.null(sysinf)) {
    os <- sysinf["sysname"]
    if (os == "Darwin") os <- "osx"
  } else {
    ## If something gets wrong with Sys.info, probe using .Platform
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))   os <- "osx"
    if (grepl("unix", R.version$os))      os <- "osx"
    if (grepl("linux-gnu", R.version$os)) os <- "linux"
    if (grepl("mingw32", R.version$os))   os <- "windows"
  }
  tolower(os)
}


theta.as.mcmc.list <- function(samples,start=1,end=NA,split=FALSE,thin=1)
  # extract theta and convert theta to mcmc.list object
  # split doubles number of chains, first and second half
{
  if ( is.na(end) ) end <- dim(samples$theta)[3]
  n.chains <- dim(samples$theta)[1]
  lst <- vector(mode="list",length=n.chains*ifelse(split,2,1))
  indx <- start:end
  if (split) is.in <- !as.logical(indx %% 2) else 
             is.in <- rep(TRUE,length(indx))
  if (split) {
    not.is.in <- !is.in
    if ( sum(not.is.in) > sum(is.in) ) {
      not.is.in[1:2][not.is.in[1:2]] <- FALSE  
    } else if ( sum(not.is.in) < sum(is.in) ) {
      is.in[1:2][is.in[1:2]] <- FALSE  
    }
  }
  for (i in 1:n.chains)
    lst[[i]] <- mcmc(t(samples$theta[i,,indx[is.in]]),thin=thin)
  if ( split ) {
    for (i in 1:n.chains)
      lst[[i+n.chains]] <- mcmc(t(samples$theta[i,,indx[not.is.in]]),thin=thin)
  }
  mcmc.list(lst)
}


phi.as.mcmc.list <- function(hyper,start=1,end=NA,split=FALSE,thin=1)
  # extract phi and convert theta to mcmc.list object
{
  # Parameters that are not constants
  ok1 <- lapply(hyper$pp.prior,function(x){
    lapply(x,function(y){attr(y,"dist")!="constant"})})
  ok2 <- paste(names(ok1[[2]])[unlist(ok1[[2]])],"h2",sep=".")
  ok1 <- paste(names(ok1[[1]])[unlist(ok1[[1]])],"h1",sep=".")
  if ( is.na(end) ) end <- dim(hyper$phi[[1]])[3]
  n.chains <- dim(hyper$h_log_likelihoods)[2] 
  lst <- vector(mode="list",length=n.chains)
  indx <- start:end
  if (split) is.in <- !as.logical(indx %% 2) else 
             is.in <- rep(TRUE,length(indx))
  if (split) {
    not.is.in <- !is.in
    if ( sum(not.is.in) > sum(is.in) ) {
      not.is.in[1:2][not.is.in[1:2]] <- FALSE  
    } else if ( sum(not.is.in) < sum(is.in) ) {
      is.in[1:2][is.in[1:2]] <- FALSE  
    }
  }
  for (i in 1:n.chains) {
    tmp1 <- t(hyper$phi[[1]][i,,indx[is.in]])
    dimnames(tmp1)[[2]] <- paste(dimnames(tmp1)[[2]],"h1",sep=".")
    tmp1 <- tmp1[,ok1,drop=FALSE]
    tmp2 <- t(hyper$phi[[2]][i,,indx[is.in]])
    dimnames(tmp2)[[2]] <- paste(dimnames(tmp2)[[2]],"h2",sep=".")
    tmp2 <- tmp2[,ok2,drop=FALSE]
    # Remove cases with !has.sigma
    tmp2 <- tmp2[,!apply(tmp2,2,function(x){all(is.na(x))})]
    lst[[i]] <- mcmc(cbind(tmp1,tmp2),thin=thin)
  }
  if (split) {
    for (i in 1:n.chains) {
      tmp1 <- t(hyper$phi[[1]][i,,indx[not.is.in]])
      dimnames(tmp1)[[2]] <- paste(dimnames(tmp1)[[2]],"h1",sep=".")
      tmp1 <- tmp1[,ok1,drop=FALSE]
      tmp2 <- t(hyper$phi[[2]][i,,indx[not.is.in]])
      dimnames(tmp2)[[2]] <- paste(dimnames(tmp2)[[2]],"h2",sep=".")
      tmp2 <- tmp2[,ok2,drop=FALSE]
      # Remove cases with !has.sigma
      tmp2 <- tmp2[,!apply(tmp2,2,function(x){all(is.na(x))})]
      lst[[i+n.chains]] <- mcmc(cbind(tmp1,tmp2),thin=thin)
    }
  }
  mcmc.list(lst)
}


p.fun.dmc <- function(samples,fun,hyper=FALSE,ptype=1,pnams=NA) 
  # applies fun to samples
  {
   if (!hyper) {
     if (any(is.na(pnams))) 
       pnams <- dimnames(samples$phi)[[2]] 
     as.vector(apply(samples$theta[,pnams,],c(1,3),fun)) 
   } else {
     if (any(is.na(pnams))) 
       pnams <- dimnames(attr(samples,"hyper")$phi[[ptype]])[[2]] 
     as.vector(apply(attr(samples,"hyper")$phi[[ptype]][,pnams,],c(1,3),fun))
   }
}


### Stuck chains ----

pick.stuck.dmc <- function(samples,hyper=FALSE,cut=10,start=1,end=NA,
                           verbose=FALSE,digits=2) 
  # return index of stuck chains (deviaton of ll from median < -cut)  
{
  if ( hyper ) { # MG: TODO: Haven't split up hyper$pll yet
    hyper <- attr(samples,"hyper")
    if ( is.na(end) ) end <- dim(hyper$phi[[1]])[3]      
    if ( end <= start )
      stop("End must be greater than start")
    mean.ll <- apply(hyper$h_log_likelihoods[start:end,] + 
                     hyper$h_summed_log_prior[start:end,],2,mean) 
    names(mean.ll) <- 1:length(mean.ll)
  } else {
    if ( is.na(end) ) end <- samples$nmc
    if ( end <= start ) stop("End must be greater than start")

    mean.ll <- apply(samples$log_likelihoods[start:end,] + 
                     samples$summed_log_prior[start:end,],2,mean) 
    names(mean.ll) <- 1:length(mean.ll)
  }
  dev <- -(sort(mean.ll)-median(mean.ll))
  bad <- as.numeric(names(dev)[dev>cut])
  if (verbose) {
    cat("Deviation of mean chain log-likelihood from median of means\n")
    print(round(dev,digits))
    cat("Bad chains: ")
    if (length(bad)==0) cat("None\n") else {
      cat(bad)
      cat("\n")
    }
  }
  bad    
}


h.pick.stuck.dmc <- function(samples,cut=10,start=1,end=NA,
                           verbose=FALSE,digits=2)
# List named by subject with chain numbers that are stuck.
{
  tmp <- lapply(lapply(samples,
    pick.stuck.dmc,start=start,cut=cut,digits=digits,end=end),
    function(x){if (length(x)>0) x})
  empty <- unlist(lapply(tmp,is.null))
  if (any(empty)) tmp <- tmp[-c(1:length(tmp))[empty]] 
  if (any(names(attributes(samples))=="hyper"))
  tmp <- c(hyper=pick.stuck.dmc(samples,hyper=TRUE,
                            start=start,cut=cut,digits=digits,end=end))
  tmp
}


### Diagnostics ----


gelman.diag.mpsrf <- function(mcmclist,autoburnin,transform) 
  # robust version ONLY USED IN sampling.R IN run.converge.dmc 
  # SHOULD BE ROLLED OUT OVER FOLLOWING FUNCITONS TO AVOID CRASHES OF 
  # AUTO PROCEDURES.
{
  gd <- try(gelman.diag(mcmclist,
          autoburnin=autoburnin,transform=transform),silent=TRUE)
  if (class(gd)=="try-error") Inf else gd$mpsrf
}


# hyper=FALSE;digits=2;start=1
# autoburnin=FALSE;transform=TRUE;end=NA
# split=TRUE;verbose=TRUE
gelman.diag.dmc <- function(samples,hyper=FALSE,digits=2,start=1,
                            autoburnin=FALSE,transform=TRUE,end=NA,
                            split=TRUE,verbose=TRUE)
  # R hat for one or list of subjects or hyper level. 
  # split doubles the number of chains by spliting into 1st and 2nd halves.
{
  if ( hyper ) {
    hyper <- attr(samples,"hyper")
    if (is.na(end)) end <- hyper$nmc
    gelman.diag(phi.as.mcmc.list(hyper,start=start,end=end,split=split),
                autoburnin=autoburnin,transform=transform)
  } else {
    if ( !is.null(samples$theta) ) {
      if (is.na(end)) end <- samples$nmc
      gelman.diag(theta.as.mcmc.list(samples,start=start,end=end,split=split),
                  autoburnin=autoburnin,transform=transform) 
    } else {
      out <- vector(mode="list",length=length(samples))
      start <- rep(start,length.out=length(samples))
      end <- rep(end,length.out=length(samples))
      names(out) <- names(samples)
      for (i in 1:length(samples)) {
        if ( is.na(end[i]) ) end[i] <- samples[[i]]$nmc
        out[[i]] <- gelman.diag(theta.as.mcmc.list(samples[[i]],start=start[i],
          end=end[i],split=split),autoburnin=autoburnin,transform=transform) 
      }
      tmp <- unlist(lapply(out,function(x){x$mpsrf}))
      if (verbose) {
        print(round(sort(tmp),digits))
        cat("Mean\n")
         print(round(mean(tmp),digits))
      }
      invisible(out)
    }
  }
}


h.gelman.diag.dmc <- function(samples,digits=2,start=1,verbose=TRUE,
                            autoburnin=FALSE,transform=TRUE,end=NA,
                            split=TRUE,...)

{
  gd <- sort(unlist(lapply(gelman.diag.dmc(samples,start=start,verbose=FALSE,
    autoburnin=autoburnin,transform=transform,split=split,end=end),function(x){x$mpsrf})))
  if ( any(names(attributes(samples))=="hyper") ) gd <- 
      c(hyper=gelman.diag.dmc(samples,start=start,end=end,verbose=FALSE,hyper=TRUE,
        autoburnin=autoburnin,transform=transform,split=split)$mpsrf,gd)
  if (verbose) print(round(gd,digits))
  invisible(gd)
}


effectiveSize.dmc <- function(samples,hyper=FALSE,digits=0,start=1,end=NA) 
  # Effective size for either single or multiple subjects
{
  if (hyper) {
    hyper <- attr(samples,"hyper")
    if (is.na(end)) end <- hyper$nmc
    round(effectiveSize(phi.as.mcmc.list(attr(samples,"hyper"),start=start,end=end)),digits)
  } else {
    if (!is.null(samples$theta)) {
      if (is.na(end)) end <- samples$nmc
      round(effectiveSize(theta.as.mcmc.list(samples,start=start,end=end)),digits) 
    } else {
      out <- vector(mode="list",length=length(samples))
      start <- rep(start,length.out=length(samples))
      end <- rep(end,length.out=length(samples))
      names(out) <- names(samples)
      for (i in 1:length(samples)) {
        if ( is.na(end[i]) ) end[i] <- samples[[i]]$nmc
        out[[i]] <- round(
          effectiveSize(theta.as.mcmc.list(samples[[i]],start=start[i],end=end[i])),
        digits) 
      }
      out
    }
  }
}

get.thin <- function(samples,hyper=FALSE) 
  # Creates summaries to help pick level of thinning.
  {
  es <- effectiveSize.dmc(samples,hyper=hyper)
  if (hyper) {
    print(es)  
    samples <- attr(samples,"hyper")
      n <- prod(dim(samples$phi[[1]])[-2])
      cat("Thin\n")
      print(round(c(mean=n/mean(es),min=n/min(es)),0))
  } else {
    if (any(names(samples)=="theta")) {
      print(es)
      n <- prod(dim(samples$theta)[-2])
      cat("Thin\n")
      print(round(c(mean=n/mean(es),min=n/min(es)),0))
    } else {
      cat("Minimum Effective Size\n")
      print(unlist(lapply(es,min)))
      cat("Mean Effective Size\n")
      print(round(unlist(lapply(es,mean))))
      n <- unlist(lapply(samples,function(x){prod(dim(x$theta)[-2])}))
      mn <- n/unlist(lapply(es,mean))
      mi <- n/unlist(lapply(es,min)) 
      cat("Thin\n")
      print(round(rbind(mean=mn,min=mi),0))
    }
  }
}

summary.dmc <- function(samples,hyper=FALSE,digits=2,start=1,end=NA,
  hyper.means=FALSE,hyper.ci=FALSE) 
  # CODA summary of dmc parameters
{
  if (hyper) {
    hyper <- attr(samples,"hyper")
    if (is.na(end)) end <- hyper$nmc
    if (hyper.means) {
      tmp <- summary(phi.as.mcmc.list(hyper,start=start,end=end))
      matrix(tmp$statistics[,"Mean"],ncol=2,
             dimnames=list(samples[[1]]$p.names,c("h1","h2")))
    } else if (hyper.ci) {
      tmp <- summary(phi.as.mcmc.list(hyper,start=start,end=end))$quantiles[,c(1,3,5)]
      tmp <- round(cbind(tmp[1:(dim(tmp)[1]/2),],
        tmp[((dim(tmp)[1]/2)+1):dim(tmp)[1],]),2)
      dimnames(tmp) <- list(unlist(strsplit(dimnames(tmp)[[1]],".h1")),
        paste(c("MEAN: ","",""," SD:  ","",""),dimnames(tmp)[[2]],sep="") )
      round(tmp,digits)
    } else summary(phi.as.mcmc.list(hyper,start=start,end=end))
  } else {
    if (!is.null(samples$theta)) {
      if (is.na(end)) end <- samples$nmc
      summary(theta.as.mcmc.list(samples,start=start,end=end)) 
    } else {
      out <- vector(mode="list",length=length(samples))
      start <- rep(start,length.out=length(samples))
      end <- rep(end,length.out=length(samples))
      names(out) <- names(samples)
      for (i in 1:length(samples)) {
        if ( is.na(end[i]) ) end[i] <- samples[[i]]$nmc
        out[[i]] <- 
          summary(theta.as.mcmc.list(samples[[i]],start=start[i],end=end[i]))
      }
      tmp <- t(data.frame(lapply(out,function(x){x[[1]][,1]})))
      tmp <- rbind(tmp,apply(tmp,2,mean))
      row.names(tmp) <- c(names(samples),"Mean")
      print(round(tmp,digits))
      invisible(out)
    }
  } 
}


check.recovery.dmc <- function(samples,p.vector=NULL,digits=2,verbose=TRUE) 
  # Tables parameter recovery stats for a single subject  
{
  qs <- summary.dmc(samples)$quantiles
  p.names <- row.names(qs)
  if (!is.null(p.vector) && (!all(dimnames(qs)[[1]] %in% names(p.vector))))
    stop("Names of p.vector do not match parameter names in samples")
  est <- qs[,"50%"]
  lo <- qs[,"2.5%"]
  hi <- qs[,"97.5%"]
  out <- rbind('2.5% Estimate'=lo,'50% Estimate'=est,
               '97.5% Estimate'=hi)
  if (!is.null(p.vector)) out <- rbind('True'=p.vector[p.names],out,
    'Median-True'= est- p.vector[p.names])
  attr(out,"ci50") <- qs[,c("25%","75%")]
  if (verbose) print(round(out,digits))
  invisible(out)
}


h.check.recovery.dmc <- function(samples,ps=NULL,hyper=FALSE,verbose=TRUE,
  digits=3,ptype=1,do.coverage=FALSE) 
{
  if (hyper) {
    samples <- list(theta=attr(samples,"hyper")$phi[[ptype]])
    samples$nmc <- dim(samples$theta)[3]
    out <- check.recovery.dmc(samples,ps,digits=digits,verbose=verbose)    
  } else {
    out <- samples
    for (i in 1:length(samples)) {
      if (is.matrix(ps)) psi <- ps[i,] else psi <- ps
      out[[i]] <- check.recovery.dmc(samples[[i]],psi,verbose=FALSE)
      if (i==1) {
        av <- out[[i]]
        attr(av,"ci50") <- NULL
      } else av <- out[[i]] + av
    }
    if (do.coverage) {
      coverage <- 100*apply(do.call(rbind,lapply(out,function(x){
        apply(x,2,function(y){y[2]<y[1] & y[4]>y[1]})})),2,sum)
      coverage50 <- 100*apply(do.call(rbind,lapply(out,function(x){
        x[2,] <- attr(x,"ci50")[,1]
        x[4,] <- attr(x,"ci50")[,2]
        apply(x,2,function(y){y[2]<y[1] & y[4]>y[1]})
      })),2,sum)
      out <- rbind(av,Coverage50=coverage50,Coverage95=coverage)/length(samples)
    } else out <- av/length(samples)
    if (verbose) print(round(out,digits=digits))
  }
  invisible(out)  
}


### Model selection ----

pll.dmc <- function(samples,hyper=FALSE,start=1,end=NA,prior=FALSE,like=FALSE,subject=NA)
  # extracts posterior log-likelihoods
{
  if ( hyper ) {
    hyper <- attr(samples,"hyper")
    if (is.null(hyper)) 
      stop("There are no hyper-parameters to plot.")
    if ( is.na(end) ) end <- hyper$nmc
    if ( end <= start )
      stop("End must be greater than start")
    if (prior) out <- hyper$h_summed_log_prior[start:end,] else
    if (like) out <-  hyper$h_log_likelihoods[start:end,] else
    out <- hyper$h_summed_log_prior[start:end,] + 
           hyper$h_log_likelihoods[start:end,]
    dimnames(out) <- list(start:end,1:dim(out)[2])
  } else {
    if ( is.null(samples$theta) & !is.na(subject) ) 
      samples <- samples[[subject]]
    if ( is.null(samples$theta) ) { # multiple subjects
      if ( is.na(end) ) end <- samples[[1]]$nmc
      if ( end <= start ) stop("End must be greater than start")
      out <- lapply(samples,function(x){
      if (prior) out <- x$summed_log_prior[start:end,] else
        if (like) out <- x$log_likelihoods[start:end,] else
          out <- x$summed_log_prior[start:end,] + x$log_likelihoods[start:end,]
        dimnames(out) <- list(start:end,1:dim(out)[2])
        out
     }) 
    } else { # single subejct
      if ( is.na(end) ) end <- samples$nmc
      if ( end <= start ) stop("End must be greater than start")
      if (prior) out <- samples$summed_log_prior[start:end,] else
      if (like) out <-  samples$log_likelihoods[start:end,] else
      out <- samples$summed_log_prior[start:end,] + 
             samples$log_likelihoods[start:end,]
      dimnames(out) <- list(start:end,1:dim(out)[2])
    }
  }
  out  
}



Dstats.dmc <- function(samples,save=FALSE,fast=TRUE)
  # Calculates and returns list with mean (meanD), variance (varD) and min (minD) 
  # of Deviance, deviance of mean theta (Dmean) and (optionally) matrix of 
  # deviances for each theta. If fast=FALSE frist recalcualtes.
{
  if (fast)
    D <- -2*samples$log_likelihoods else 
    D <- apply(samples$theta,c(3,1),function(x){ 
      -2*sum(log(likelihood.dmc(x,samples$data)))
    })
  
  mtheta <- apply(samples$theta,2,mean)
  Dmean <- -2*sum(log(likelihood.dmc(mtheta,samples$data)))
  minD <- min(D)
  if (save) 
    list(np=dim(samples$theta)[2],meanD=mean(D),varD=var(as.vector(D)),
         minD=min(D),Dmean=Dmean,D=D) else
    list(np=dim(samples$theta)[2],meanD=mean(D),varD=var(as.vector(D)),
         minD=min(D),Dmean=Dmean,D=NA)
}


pd.dmc <- function(ds) 
  # effective number of parameters calculated by mean, min and var methods  
{
  list(Pmean=ds$meanD-ds$Dmean,Pmin=ds$meanD-ds$minD,Pvar=ds$varD/2)
}


posterior.lr.dmc <- function(D1,D2,main="",
                             plot=FALSE,plot.density=TRUE)
  # Aitkins posterior deviance likelihood ratio test
{
  if (is.list(D1)) D1 <- D1$D
  if (is.list(D2)) D2 <- D2$D
  n <- pmin(length(D1),length(D2))
  dD <- D1[1:n]-D2[1:n] 
  if (plot) if (plot.density) 
    plot(density(dD),xlab="D1-D2",main=main) else
      hist(dD,breaks="fd",xlab="D1-D2",main=main)
  if (min(dD)<0 & max(dD)>0) abline(v=0)
  c(pD1=mean(dD<0))
}


IC.dmc <- function(ds=NULL,samples=NULL,DIC=FALSE,fast=TRUE,use.pd=NA)
  # Calcualte IC1 (i.e., BPIC, default) or DIC  
{
  if (is.null(ds)) if (is.null(samples)) 
    stop("Must supply samples or deviance statistics") else 
      ds <- Dstats.dmc(samples,save=TRUE,fast=fast)    
    pds <- pd.dmc(ds)
    if (is.na(use.pd)) {
      if (ds$minD < ds$Dmean) pd <- pds$Pmin else pd <- pds$Pmean   
    } else {
      if (use.pd %in% names(pds))  
        pd <- pds[[use.pd]] else
          stop(paste("use.pd must be one of:",paste(names(pds),collapse=",")))
    }      
    if (DIC) ds$meanD+pd else ds$meanD+2*pd
}


h.IC.dmc <- function(hsamples,DIC=FALSE,fast=TRUE,use.pd=NA) 
  # Applies IC.dmc each suubject, prints sum.
  # Invisibly returns results for each subject.    
{
  if ( any(names(hsamples)=="theta") )
    stop("For a single subject use Dstats.dmc")
  ds <- lapply(hsamples,Dstats.dmc,fast=fast,save=FALSE)    
  pds <- lapply(ds,pd.dmc)
  if ( is.na(use.pd) ) {
    pd <- numeric(length(pds))
    for (i in 1:length(pds)) {
      if (ds[[i]]$minD < ds[[i]]$Dmean) 
        pd[i] <- pds[[i]]$Pmin else 
        pd[i] <- pds[[i]]$Pmean  
    }
  } else {
    if ( use.pd %in% names(pds[[1]]) )  
      pd <- unlist(lapply(pds,function(x){x[[use.pd]]})) else
        stop(paste("use.pd must be one of:",paste(names(pds[[1]]),collapse=",")))
  } 
  ICs <- numeric(length(ds))
  for (i in 1:length(ds)) if (DIC)  
    ICs[i] <- ds[[i]]$meanD+pd[i] else 
    ICs[i] <- ds[[i]]$meanD+2*pd[i]
  Dmean <- unlist(lapply(ds,function(x){x$Dmean}))
  if (DIC) cat("Summed Minimum Deviance and DIC\n") else 
           cat("Summed Minimum Deviance and BPIC\n")
  print(c(sum(Dmean),sum(ICs)))
  invisible(cbind(MinD=Dmean,IC=ICs))
}


#### Posterior Likelihoods ----

trial_log_likes <- function(samples,thin_pointwise=1,get.size=FALSE,
                            chain_dot=TRUE,subject_dot=FALSE) 
  # Get pointwise log-likelihoods  
{
  
  n.trials <- dim(samples$data)[1]
  
  nmc_thin <- seq(thin_pointwise,samples$nmc,by=thin_pointwise)
  trial_log_likes  <- 
    array(-Inf,c(length(nmc_thin),samples$n.chains, dim(samples$data)[1]))
  dimnames(trial_log_likes) <- list(nmc_thin,NULL,NULL)
  if (get.size) return(trial_log_likes)
  if (chain_dot) cat("Processing chains: ")
  for (j in 1:samples$n.chains) {
    for (i in nmc_thin)
      trial_log_likes[as.character(i),j,]  <- 
          log.likelihood(samples$theta[j,,i],samples$data)
      if (chain_dot) cat(".")
  }
  if (chain_dot) cat("\n")
  if (subject_dot) cat(".")
  trial_log_likes 
}


# samples <- hsamples0.0; thin_pointwise=50; max_size=1e8

group_trial_log_likes <- function(samples,thin_pointwise=1,max_size=1e8,
  cores=1,slaveOutfile=NULL,get.size=FALSE) 
  # extracts trial_log_likes from a list of subject fits and concatanates
  {

  if (get.size) {
    cat("Calculating\n")
    chain_dot <- TRUE
  } else {
    cat("Processing subjects: ")
    chain_dot <- FALSE
  }
  tll <- trial_log_likes(samples[[1]],thin_pointwise=thin_pointwise,
                         chain_dot=chain_dot,subject_dot=TRUE,get.size=get.size)
  sdim <- dim(tll)
  size <- sum(length(samples)*prod(sdim))
  if (get.size) {
    nsamp <- matrix(unlist(lapply(samples,function(x){dim(x$log_likelihoods)})),nrow=2)
    dimnames(nsamp) <- list(c("Iterations","Chains"),names(samples))
    cat("Available log-likelihoods\n")
    print(nsamp)
    names(sdim) <- c("Iterations","Chains","Trials")
    cat("Subject 1\n")
    print(sdim)
    cat("Total log-likelihoods (millions)\n")
    print(round(size/1e6,2))
    return()
  }
  if (size>max_size) stop(paste("Output size",size,
                                "too large,adjust max_size (",max_size,")"))
  
  os <- get.os()
  if ( cores>1 & os=="windows") {
    cat("No progress reporting with multicore\n")  
    require(snowfall,quietly=TRUE)
      sfInit(parallel=TRUE, cpus=cores, type="SOCK",slaveOutfile=slaveOutfile)
      sfLibrary(msm)
      sfLibrary(rtdists)
      sfLibrary(statmod)
      sfLibrary(pracma)
      sfExportAll()
      tlls <- sfLapply(samples[-1],trial_log_likes,thin_pointwise=thin_pointwise,
                 chain_dot=TRUE,subject_dot=FALSE)
      sfStop()
  } else if (cores>1) {
    cat("No progress reporting with multicore\n")  
    require(parallel, quietly=TRUE)
     tlls <- mclapply(samples[-1],trial_log_likes,
      thin_pointwise=thin_pointwise, chain_dot=FALSE,subject_dot=FALSE,
      mc.cores=cores)
  } else {
      tlls <- lapply(samples[-1],trial_log_likes,thin_pointwise=thin_pointwise,
                 chain_dot=FALSE,subject_dot=TRUE)
  }

  tlls[[length(samples)]] <- tll
  sdims <- cbind(matrix(unlist(lapply(tlls,dim)),nrow=3))
  nmc <- min(sdims[1,])
  if ( !all(sdims[1,1]==sdims[1,-1]) )
    warning(paste("Subjects do not all have the same number of interations, using first",nmc,"for all."))
  if ( !all(sdims[2,1]==sdims[2,-1]) )
    stop("Subjects must have the same number of chains")
  out <- array(dim=c(nmc,dim(tlls[[1]])[2],sum(sdims[3,])))
  start <- 1; end <- sdims[3,1]
  for (i in 1:length(samples)) {
    out[,,start:end] <- tlls[[i]][1:nmc,,]
    if (i<length(samples)) {
      start <- end+1
      end <- start - 1 + sdims[3,i+1]
    }
  }
  out
}


group_subject_log_likes <- function(samples) 
  # extracts summed_log_likes from a list of subject fits and concatanates
  {

  tlls <- lapply(samples,function(x){x$log_likelihoods})
  sdims <- matrix(unlist(lapply(tlls,dim)),nrow=2)
  if ( !all(sdims[2,1]==sdims[2,-1]) )
    stop("Subjects must have the same number of chains")
  nmc <- min(sdims[1,])
  if ( !all(sdims[1,1]==sdims[1,-1]) )
    warning(paste("Subjects do not all have the same number of interations, using first",nmc,"for all."))
  out <- array(dim=c(nmc,dim(tlls[[1]])[2],length(samples)))
  for (i in 1:length(samples)) 
    out[,,i] <- tlls[[i]][1:nmc,]
  out
}



wIC.dmc <- function(ds=NULL,samples=NULL,
                    DIC=FALSE,fast=TRUE,use.pd=NA,...)
  # Calculate weights for a set of models  
{
  
  ICs <- function(samples,DIC=FALSE,fast=TRUE,use.pd=NA)
    IC.dmc(samples=samples,DIC=DIC,fast=fast,use.pd=use.pd)
  
  if (is.null(ds)) if (is.null(samples)) 
    stop("Must supply samples or deviance statistics") 
  if (is.null(samples))
    ics <- unlist(lapply(ds,IC.dmc,DIC=DIC,fast=fast,use.pd=use.pd)) else
      ics <- unlist(lapply(ds,ICs,DIC=DIC,fast=fast,use.pd=use.pd))
    d = ics-min(ics)
    w = exp(-d/2)/sum(exp(-d/2))
    if (!is.null(names(ics))) mnams=names(ics) else
      mnams <- paste("M",1:length(d),sep="")
    print(t(matrix(c(d,w),nrow=length(d),
             dimnames=list(mnams,c("IC-min","w")))),...)
}


waic.dmc <- function(trial_log_likes,mc_se=FALSE,save=FALSE,...) 
  # Calclaute WAIC  
{
  out <- waic(matrix(trial_log_likes,ncol=dim(trial_log_likes)[3]))
  if (mc_se) {
    n.chains <- dim(trial_log_likes)[2]
    waics <- vector(mode="list",length=n.chains)
    for (i in 1:n.chains) waics[[i]] <- waic(trial_log_likes[,i,])
    mc <- sd(unlist(lapply(waics,function(x){x$waic})))/sqrt(n.chains)
    print(c(p=out$p_waic,se_p=out$se_p,waic=out$waic,se_waic=out$se_waic,
      mc_se_waic=mc),...) 
  } else 
    print(c(p=out$p_waic,se_p=out$se_p,waic=out$waic,se_waic=out$se_waic),...)
  if (save) out
}


looic.dmc <- function(trial_log_likes,mc_se=FALSE,save=FALSE,...) 
  # Calcuate looic 
{
  out <- loo(matrix(trial_log_likes,ncol=dim(trial_log_likes)[3]))
  if (all(out$pareto_k<.5))
    cat("All Pareto k estimates OK (k < 0.5)\n") else {
    msg <- "See PSIS-LOO description (?'loo-package') for more information" 
    if (any(out$pareto_k>1))
      msg1 <- paste(sum(out$pareto_k>1)," (",
                    round(100*sum(out$pareto_k>1)/length(out$pareto_k),1),
                    "%) Pareto k estimates greater than 1\n",sep="") else
      msg1 <- paste(sum(out$pareto_k>0.5)," (",
                    round(100*sum(out$pareto_k>1)/length(out$pareto_k),1),
                    "%) Pareto k estimates estimates between 0.5 and 1\n",sep="")
     warning(msg1,msg) 
  }
  if (mc_se) {
    n.chains <- dim(trial_log_likes)[2]
    loos <- vector(mode="list",length=n.chains)
    for (i in 1:n.chains) loos[[i]] <- loo(trial_log_likes[,i,])
    mc <- sd(unlist(lapply(loos,function(x){x$looic})))/sqrt(n.chains)
    print(c(p=out$p_loo,se_p=out$se_p,looic=out$looic,se_looic=out$se_looic,
      mc_se_loo=mc),...) 
  } else 
    print(c(p=out$p_loo,se_p=out$se_p,looic=out$looic,se_looic=out$se_looic),...)
  if (save) out
}


loocompare.dmc <- function(loo1,loo2=NULL,...) 
  # Model comparision of objects produced by waic.dmc or looic.dmc  
{
  if ( !is.null(loo2) ) {
    tmp <- compare(loo1,loo2)
    out <- c(waic_diff=-2*tmp[[1]],se=2*tmp[[2]]) 
    print(out,...)
  } else {
    if ( !is.list(loo1) )
      stop("loo1 argument must be a list if loo2 not given")
    if ( any(names(loo1[[1]])=="looic") ) icnam <- "looic" else 
                                          icnam <- "waic"
    ics <- unlist(lapply(loo1,function(x){x[[icnam]]}))
    d = ics-min(ics)
    w = exp(-d/2)/sum(exp(-d/2))
    if (!is.null(names(ics))) mnams=names(ics) else
      mnams <- paste("M",1:length(d),sep="")
    print(t(matrix(c(d,w),nrow=length(d),
      dimnames=list(mnams,c("IC-min","w")))),...)
  }
}


#### Priors from posteriors ----

make.tnorm.prior <- function(p.prior,theta,scale.sd=1,hpar=1) 
  # Takes a samples array and creates fitting truncted normal priors.

{

  est.tnorm <- function(dat,lower=0,upper=Inf) 
    # Fit to truncated normal.  
  {
    if (any(dat<lower | dat>upper)) stop("Lower and/or upper do not match data.")
    obj <- function(par,dat,low,up) { 
      out <- suppressWarnings(-sum(dtnorm(
        dat,mean=par[1],sd=par[2],lower=low,upper=up,log=TRUE)))
      if (is.na(out)) out <- -Inf
      pmax(out,log(1e-10))  
    }
    par  <- c(mean=mean(dat),sd=sd(dat))
    if (lower==-Inf & upper==Inf) par else
      optim(par,obj,dat=dat,method="L-BFGS-B",low=lower,up=upper)$par
  }

  p.names <- names(p.prior)
  if ( is.null(dim(theta)) ) 
    if ( any(names(theta)=="theta") ) theta <- theta$theta else {
      theta <- attr(theta,"hyper")$phi[[hpar]]
    }
  if ( is.null(theta) || !all(sort(dimnames(theta)[[2]])==sort(p.names)) ) 
    stop("theta must be a parameter array or samples object matching p.prior.")
  scale <- rep(scale.sd,length.out=length(p.names))
  names(scale) <- p.names
  if (!all(sort(p.names)==sort(names(p.prior))))
    stop("p.prior does not match parameters")
  for (i in p.names) {
    lower <- p.prior[[i]]$lower; if (is.null(lower)) lower <- -Inf
    upper <- p.prior[[i]]$upper; if (is.null(upper)) upper <- Inf
    pars <- est.tnorm(as.vector(theta[,i,]),lower=lower,upper=upper)
    p.prior[[i]] <- 
      list(mean=pars[1],sd=pars[2]*scale[i],lower=lower,upper=upper,log=TRUE)
    attr(p.prior[[i]],"dist") <- "tnorm"
    attr(p.prior[[i]],"untrans") <- "identity"
  }
  p.prior
}

get.ppars <- function(p.prior) {
  out <- data.frame(do.call(rbind,lapply(p.prior,function(x){c(x$mean,x$sd)})))
  colnames(out) <- c("mean", "sd")
  out
}

assign.ppars <- function(p.prior,ppars) {
  for (i in names(p.prior)) {
    p.prior[[i]]$mean <- ppars[i,"mean"]
    p.prior[[i]]$sd <- ppars[i,"sd"]
  }
  p.prior
}


### Plausible values ----

posteriorRho <- function(r, n, npoints=100, kappa=1)
  # Code provided by Dora Matzke, March 2016, from Alexander Ly
  # Reformatted into a single funciton. kappa=1 implies uniform prior.
  # Picks smart grid of npoints points concentrating around the density peak. 
  # Returns approxfun for the unnormalized density. 
{

  .bf10Exact <- function(n, r, kappa=1) {
	  # Ly et al 2015
	  # This is the exact result with symmetric beta prior on rho
	  # with parameter alpha. If kappa = 1 then uniform prior on rho
	  #
	  if (n <= 2){
		  return(1)
	  } else if (any(is.na(r))){
		  return(NaN)
	  }
	  # TODO: use which
	  check.r <- abs(r) >= 1 # check whether |r| >= 1
	  if (kappa >= 1 && n > 2 && check.r) {
		  return(Inf)
	  }

	  log.hyper.term <- log(hypergeo::genhypergeo(U=c((n-1)/2, (n-1)/2), 
	                                              L=((n+2/kappa)/2), z=r^2))
	  log.result <- log(2^(1-2/kappa))+0.5*log(pi)-lbeta(1/kappa, 1/kappa)+
		lgamma((n+2/kappa-1)/2)-lgamma((n+2/kappa)/2)+log.hyper.term
	  real.result <- exp(Re(log.result))
	  return(real.result)
  }
  
  .jeffreysApproxH <- function(n, r, rho) {	
	  result <- ((1 - rho^(2))^(0.5*(n - 1)))/((1 - rho*r)^(n - 1 - 0.5))
	  return(result)
  }

  .bf10JeffreysIntegrate <- function(n, r, kappa=1) {
	# Jeffreys' test for whether a correlation is zero or not
	# Jeffreys (1961), pp. 289-292
	# This is the exact result, see EJ
	##
	if (n <= 2){
		return(1)
	} else if ( any(is.na(r)) ){
		return(NaN)
	}
	
	# TODO: use which
	if (n > 2 && abs(r)==1) {
		return(Inf)
	}
	hyper.term <- Re(hypergeo::genhypergeo(U=c((2*n-3)/4, (2*n-1)/4), L=(n+2/kappa)/2, z=r^2))
	log.term <- lgamma((n+2/kappa-1)/2)-lgamma((n+2/kappa)/2)-lbeta(1/kappa, 1/kappa)
	result <- sqrt(pi)*2^(1-2/kappa)*exp(log.term)*hyper.term
	return(result)
}


  # 1.0. Built-up for likelihood functions
  .aFunction <- function(n, r, rho) {
	  #hyper.term <- Re(hypergeo::hypergeo(((n-1)/2), ((n-1)/2), (1/2), (r*rho)^2))
	  hyper.term <- Re(hypergeo::genhypergeo(U=c((n-1)/2, (n-1)/2), L=(1/2), z=(r*rho)^2))
	  result <- (1-rho^2)^((n-1)/2)*hyper.term
	  return(result)
  }

  .bFunction <- function(n, r, rho) {
	  #hyper.term.1 <- Re(hypergeo::hypergeo((n/2), (n/2), (1/2), (r*rho)^2))
	  #hyper.term.2 <- Re(hypergeo::hypergeo((n/2), (n/2), (-1/2), (r*rho)^2))
	  #hyper.term.1 <- Re(hypergeo::genhypergeo(U=c(n/2, n/2), L=(1/2), z=(r*rho)^2))
	  #hyper.term.2 <- Re(hypergeo::genhypergeo(U=c(n/2, n/2), L=(-1/2), z=(r*rho)^2))
	  #result <- 2^(-1)*(1-rho^2)^((n-1)/2)*exp(log.term)*
	  #	((1-2*n*(r*rho)^2)/(r*rho)*hyper.term.1-(1-(r*rho)^2)/(r*rho)*hyper.term.2)
	  #
	  hyper.term <- Re(hypergeo::genhypergeo(U=c(n/2, n/2), L=(3/2), z=(r*rho)^2))
	  log.term <- 2*(lgamma(n/2)-lgamma((n-1)/2))+((n-1)/2)*log(1-rho^2)
	  result <- 2*r*rho*exp(log.term)*hyper.term
	  return(result)
  }

  .hFunction <- function(n, r, rho) {
	  result <- .aFunction(n, r, rho) + .bFunction(n, r, rho)
	  return(result)
  }

  .scaledBeta <- function(rho, alpha, beta){
	  result <- 1/2*dbeta((rho+1)/2, alpha, beta)
	  return(result)
  }

  .priorRho <- function(rho, kappa=1) {
	  .scaledBeta(rho, 1/kappa, 1/kappa)
  }

  fisherZ <- function(r) log((1+r)/(1-r))/2
  
  inv.fisherZ <- function(z) {K <- exp(2*z); (K-1)/(K+1)}
  
  
  # Main body
  
  # Values spaced around mode
  qs <- qlogis(seq(0,1,length.out=npoints+2)[-c(1,npoints+2)])
  rho <- c(-1,inv.fisherZ(fisherZ(r)+qs/sqrt(n)),1)
  # Get heights
  if (!is.na(r) && !r==0) {
    d <- .bf10Exact(n, r, kappa)*.hFunction(n, r, rho)*.priorRho(rho, kappa)
	} else if (!is.na(r) && r==0) {
		d <- .bf10JeffreysIntegrate(n, r, kappa)*
		  .jeffreysApproxH(n, r, rho)*.priorRho(rho, kappa)
	} else return(NA)
  # Unnormalized approximation funciton for density
  approxfun(rho,d)
}


postRav <- function(r, n, spacing=.01, kappa=1,npoints=100,save=FALSE)
  # r is a vector, returns average density. Can also save unnormalized pdfs
{
  funs <- sapply(r,posteriorRho,n=n,npoints=npoints,kappa=kappa)
  rho <- seq(-1,1,spacing)
  result <- apply(matrix(unlist(lapply(funs,function(x){
    out <- x(rho); out/sum(out)
  })),nrow=length(rho)),1,mean)
  names(result) <- seq(-1,1,spacing)
  attr(result,"n") <- n
  attr(result,"kappa") <- kappa
  if (save) attr(result,"updfs") <- funs
  result
}


postRav.Density <- function(result)
  # Produces density class object
{
  x.vals <- as.numeric(names(result))
  result <- result/(diff(range(x.vals))/length(x.vals))
  out <- list(x=x.vals,y=result,has.na=FALSE,
    data.name="postRav",call=call("postRav"),
    bw=mean(diff(x.vals)),n=attr(result,"n")) 
  class(out) <- "density"
  out
}

postRav.mean <- function(pra) {
  # Average value of object produced by posteriorRhoAverage
  sum(pra*as.numeric(names(pra)))
}

postRav.p <- function(pra,lower=-1,upper=1) {
  # probability in an (inclusive) range of posteriorRhoAverage object
  x.vals <- as.numeric(names(pra))
  sum(pra[x.vals <= upper & x.vals >= lower])
}

postRav.ci <- function(pra,interval=c(.025,.975))
{
  cs <- cumsum(pra)
  rs <- as.numeric(names(pra))
  tmp <- approx(cs,rs,interval)
  out <- tmp$y
  names(out) <- interval
  out
}

cor.plausible <- function(hsamples,p.name,cv,plot=FALSE,cor.plot=FALSE,pick=1000,
                          xlab="r",ylab="Density",main=p.name,
                          fun=NULL,...)
  # Correlations for each interation and chain with of p.name with cv$p.name  
{
  if ( !is.null(fun) ) subject.pars <- do.call(rbind,lapply(hsamples,function(x){
    pars <- aperm(x$theta,c(1,3,2));
    pars <- matrix(as.vector(pars),ncol=dim(pars)[3]);
    dimnames(pars)[[2]] <- dimnames(x$theta)[[2]];
    apply(pars,1,fun)
  })) else subject.pars <- do.call(rbind,lapply(hsamples,function(x){
    as.vector(x$theta[,p.name,])}))
  out <- apply(subject.pars,2,function(x){cor(x,cv[[p.name]])})
  if (plot) {
    if (cor.plot) {
      ylim <- c(min(subject.pars),max(subject.pars))
      ord <- order(cv[[p.name]])
      pick <- sample(1:dim(subject.pars)[2],pick)
      plot(cv[[p.name]][ord],subject.pars[,pick[1]][ord],col="grey",type="l",
        main=main,xlab=xlab,ylab=ylab,ylim=ylim)
      for (i in pick[-1]) lines(cv[[p.name]][ord],subject.pars[,i][ord],
        col="grey")
      mn <- apply(subject.pars,1,mean)
      lines(cv[[p.name]][ord],mn[ord],type="l",lwd=2)
    } else {
      dns <- density(out)
      # dns$y <- dns$y*diff(range(dns$x))/length(dns$x)
      plot(dns,xlab=xlab,ylab=ylab,main=main,...)
    }
  }
  invisible(out)
}


compare.r <- function(r1,r2,plot=FALSE,add=FALSE,n=1,...) {
  # Tests difference between r1 and r2. n is number of samples to take from each
  # population distribution. Returns probability r1>r2

  sample.pop.r <- function(fun,n=1,spacing=.01) {
    pts <- seq(-1,1,spacing)
    cdf <- fun(pts)
    cdf <- cumsum(cdf)/sum(cdf)
    approx(cdf,pts,runif(n))[[2]]
  }
  
  if ( any(names(attributes(r1))=="updfs") ) { # population case
    r12 <- unlist(lapply(attr(r1,"updfs"),sample.pop.r,n=n)) -
      unlist(lapply(attr(r2,"updfs"),sample.pop.r,n=n))
  } else {
    r12 <- r1-r2
  }
  if (plot) plot(density(r12),xlab="Correlation Difference",...)
  if (add) lines(density(r12),...)
  mean(r12>0)
}


# New version with subsampling.
# fun=NULL;pnames=NULL;pretty=NULL;hyper=FALSE;hpar=1;subsample=NA
# show.plot=FALSE;main="";xlab="";show.legend=TRUE;lpos="topleft";xlim=NULL
# show.table=TRUE;lo.p=.025;hi.p=.975;digits=3;line0=TRUE
# 
# hsamples <- sOldGo[[1]]
# pnames=c("v.block.easy.true","v.block.hard.true"); pretty=c("easy","hard")
# show.plot=TRUE;xlab="v.block TRUE"
# 
compare.p=function(hsamples,fun=NULL,pnames=NULL,pretty=NULL,hyper=FALSE,hpar=1,subsample=NA,
                   show.plot=FALSE,main="",xlab="",show.legend=TRUE,lpos="topleft",xlim=NULL,
                   show.table=TRUE,lo.p=.025,hi.p=.975,digits=3,line0=TRUE,...)
{
  
  get_iin <- function(subsample,niter) {
    if (any(is.na(subsample))) 1:niter else {
      if (max(subsample)>niter) stop("Asking for more than you have.")
      if (length(subsample)==1) sample.int(niter,subsample) else
        subsample
    }
  }

  if (is.null(fun) & is.null(pnames))
      stop("Must define one of fun and pnames arguments")
  if (is.null(fun)) {
    if (length(pnames)!=2) stop("pnames must have length 2")
    if ( !all(pnames %in% dimnames(hsamples[[1]]$theta)[[2]]) )
      stop("pnames must correspond to names of parameters in hsamples")
    fun <- function(x){c(x[pnames],-diff(x[c(pnames)]))}
  }
  if (!any(is.na(subsample))) 
    if (any(!is.numeric(subsample)) || !is.integer(subsample) || any(subsample<1))
      stop("subsample must contain only positive integers")
  if (hyper) {
    hsamples <- attr(hsamples,"hyper")
    hsamples$theta <- hsamples$phi[[hpar]]
    iin <- get_iin(subsample,dim(hsamples$theta)[3])
    av <- apply(matrix(as.vector(aperm(hsamples$theta[,,iin],c(1,3,2))),ncol=dim(hsamples$theta)[2],
                       dimnames=list(NULL,dimnames(hsamples$theta)[[2]])),1,fun,...)
    n <- dim(av)[1]
  } else {
    cmats <- lapply(hsamples,function(x){
      iin <- get_iin(subsample,dim(x$theta)[3])
      apply(matrix(as.vector(aperm(x$theta[,,iin],c(1,3,2))),ncol=dim(x$theta)[2],
                 dimnames=list(NULL,dimnames(x$theta)[[2]])),1,fun,...)
    })
    n <- dim(cmats[[1]])[1]
    ns <- unlist(lapply(cmats,function(x){dim(x)[2]}))
    if (!all(ns[1]==ns[-1]))
      warning(paste("Different number of samples per subject, using min =",min(ns)))
    ns <- min(ns)
    cmats <- lapply(cmats,function(x){x[,1:ns]})
    av <- apply(
      array(do.call(rbind,cmats),dim=c(n,length(hsamples),ns)),
      c(1,3),mean)
  }

  if (!is.null(pretty)) {
    if (length(pretty) != (n-1) )
      stop(paste("pretty must be length",dim(cmats[[1]])[1]-1))
    dn <- c(pretty,"contrast") 
  } else {
    dn <- dimnames(av)[[1]]
    dn[length(dn)] <- c("contrast")
  }
  dimnames(av) <- list(dn,NULL)
  tab <- apply(av,1,quantile,probs=c(lo.p,.5,hi.p))
  if (show.table)
    print(round(tab,digits=digits))
    p.gt.0 <- mean(av[n,]>0)
    print(round(c(p.gt.0=p.gt.0),digits=digits))
  if (show.plot) {
    plot(density(av[n,]),main=main,xlab=xlab,xlim=xlim)
    if (line0) abline(v=0,lty=2)
    if (show.legend) legend(lpos,
      paste(c(dimnames(tab)[[2]],"p>0: "),
      round(c(tab[2,],p.gt.0),digits),sep=" = "),bty="n")
  }
  attr(av,"table") <- tab
  attr(av,"p.gt.0") <- p.gt.0
  invisible(t(av))
}



# compare.p <- function(hsamples,fun=NULL,pnames=NULL,pretty=c("1","2"),hyper=FALSE,hpar=1,
#   show.plot=FALSE,main="",xlab="",show.legend=TRUE,lpos="topleft",xlim=NULL,
#   show.table=TRUE,lo.p=.025,hi.p=.975,digits=3,line0=TRUE,verbose=TRUE,...)
#   # Plots and tests difference distribution on parameters calcualted by fun then
#   # averaged over subject.
#   # fun must return each element of the contrast and then the difference (i.e,. a 3 vector).
#   # OR difference between first and second parameter in pnames. 
#   # hsamples is the result of fixed or random effects sampling
# {
#   if (is.null(fun) & is.null(pnames))
#     stop("Must define one of fun and pnames arguments")
#   if ( is.null(fun) ) {
#       if ( !all(pnames %in% dimnames(hsamples[[1]]$theta)[[2]]) )
#         stop("pnames must correspond to names of parameters in hsamples")
#       if (length(pnames)==1) fun <- function(x){x[pnames]} else
#       if (length(pnames)==2) 
#         fun <- function(x){c(x[pnames],-diff(x[c(pnames)]))} else
#       stop("pnames must have length 1 or 2")
#   }
#   if ( hyper ) {
#     hsamples <- attr(hsamples,"hyper")
#     hsamples$theta <- hsamples$phi[[hpar]]
#     av <- apply(matrix(as.vector(aperm(hsamples$theta,c(1,3,2))),ncol=dim(hsamples$theta)[2],
#         dimnames=list(NULL,dimnames(hsamples$theta)[[2]])),1,fun,...)
#     n <- dim(av)[1]
#   } else {
#     cmats <- lapply(hsamples,function(x){
#       apply(matrix(as.vector(aperm(x$theta,c(1,3,2))),ncol=dim(x$theta)[2],
#         dimnames=list(NULL,dimnames(x$theta)[[2]])),1,fun,...)
#     })
#     n <- dim(cmats[[1]])[1]
#     if ( is.null(n) )  {
#       ns <- unlist(lapply(cmats,length))
#       if (!all(ns[1]==ns[-1]))
#         warning(paste("Different number of samples per subject, using min =",min(ns)))
#       ns <- min(ns)
#       cmats <- lapply(cmats,function(x){x[1:ns]})
#       av <- apply(do.call(rbind,cmats),2,mean)  
#     } else {
#       ns <- unlist(lapply(cmats,function(x){dim(x)[2]}))
#       if (!all(ns[1]==ns[-1]))
#         warning(paste("Different number of samples per subject, using min =",min(ns)))
#       ns <- min(ns)
#       cmats <- lapply(cmats,function(x){x[,1:ns]})
#       av <- apply(
#         array(do.call(rbind,cmats),dim=c(n,length(hsamples),ns)),
#         c(1,3),mean)
#     }
#   }
#   if ( is.null(n) ) av else {
#     if (!is.null(pretty)) {
#       if (length(pretty) != (n-1) )
#         stop(paste("pretty must be length",dim(cmats[[1]])[1]-1))
#       dn <- c(pretty,"contrast") 
#     } else {
#       dn <- dimnames(av)[[1]]
#       dn[length(dn)] <- c("contrast")
#     }
#     dimnames(av) <- list(dn,NULL)
#     tab <- apply(av,1,quantile,probs=c(lo.p,.5,hi.p))
#     if (verbose & show.table)
#       print(round(tab,digits=digits))
#     p.gt.0 <- mean(av[n,]>0)
#     if (verbose) print(round(c(p.gt.0=p.gt.0),digits=digits))
#     if (show.plot) {
#       plot(density(av[n,]),main=main,xlab=xlab,xlim=xlim)
#       if (line0) abline(v=0,lty=2)
#       if (show.legend) legend(lpos,
#         paste(c(dimnames(tab)[[2]],"p>0: "),
#           round(c(tab[2,],p.gt.0),digits),sep=" = "),bty="n")
#     }
#     attr(av,"table") <- tab
#     attr(av,"p.gt.0") <- p.gt.0
#     invisible(t(av))
#   }
# }

 
# between.fun=NULL;between.names=NULL;within.pnames=NULL;within.fun=NULL
# hyper=TRUE;hpar=1;pretty=NULL; show.plot=TRUE;main=""
# 
# within.pnames="t0"
# 
# fits=list(A1=hA1,A2=hA2)
# between.names=c("A1","A2")
# within.fun=function(x){mean(x[c("B.new","B.old")])}

compare.ps <- function(fits,between.fun=NULL,within.pnames=NULL,within.fun=NULL,
  hyper=TRUE,hpar=1,pretty=NULL,show.plot=TRUE,main="",xlab="",show.legend=TRUE,lpos="topleft",xlim=NULL,
  show.table=TRUE,lo.p=.025,hi.p=.975,digits=3,line0=TRUE,...)
  # Tests differences in parameters among two or more fits  
{
  if (length(fits)!=2) stop("Fits list must be of length two.")
  if (!is.null(within.pnames) ){
    if (length(within.pnames) != 1)
      stop("If specified within.pnames can only be length 1")
    if (!(within.pnames %in% fits[[1]][[1]]$p.names) | 
        !(within.pnames %in% fits[[2]][[1]]$p.names) )
      stop("within.pnames argument not in fits object parameter names.")
  } 
  
  if ( is.null(between.fun) ) {
    if ( is.null(names(fits)) ) 
      between.names <- c("G1","G2") else
      between.names <- names(fits)
  } 
  avs <- do.call(rbind,lapply(fits,compare.p,fun=within.fun,pnames=within.pnames,
         hyper=hyper,hpar=hpar))
  if ( is.null(between.fun) ) {
    if (is.null(names(fits)))
      stop("When using between.names fits must be named")
    if ( !all(between.names %in% names(fits)) )
      stop("between.names must correspond to names of fits")
    if (length(between.names)!=2) 
      stop("between.names must have length 2")
    fun <- function(x){c(x[between.names],-diff(x[c(between.names)]))}
    dimnames(avs)[[1]] <- names(fits)
  } else fun <- between.fun
  
  av <- apply(avs,2,fun)  
  n <- dim(av)[1]
  if ( !is.null(pretty) ) {
    if ( length(pretty) != (n-1) )
      stop(paste("pretty must be length",dim(cmats[[1]])[1]-1))
    dn <- c(pretty,"contrast") 
  } else {
    dn <- dimnames(av)[[1]]
    dn[length(dn)] <- c("contrast")
  }
  dimnames(av) <- list(dn,NULL)
  tab <- apply(av,1,quantile,probs=c(lo.p,.5,hi.p))
  if (show.table)
    print(round(tab,digits=digits))
  p.gt.0 <- mean(av[n,]>0)
  print(round(c(p.gt.0=p.gt.0),digits=digits))
  if (show.plot) {
    plot(density(av[n,]),main=main,xlab=xlab,xlim=xlim)
    if (line0) abline(v=0,lty=2)
    if (show.legend) legend(lpos,
      paste(c(dimnames(tab)[[2]],"p>0: "),
        round(c(tab[2,],p.gt.0),digits),sep=" = "),bty="n")
  }
  attr(av,"table") <- tab
  attr(av,"p.gt.0") <- p.gt.0
  invisible(t(av))
}

subject.average.ci <- function(hsamples,fun=NULL,pnames=NULL,lo.p=.025,hi.p=.975,pretty=NULL,...)
{
  if (is.null(fun) & is.null(pnames))
    pnames <- dimnames(hsamples[[1]]$theta)[[2]]
  if (is.null(pretty)) pretty <- pnames
  if ( is.null(fun) ) {
      if ( !all(pnames %in% dimnames(hsamples[[1]]$theta)[[2]]) )
        stop("pnames must correspond to names of parameters in hsamples")
      fun <- function(x){x[pnames]} 
  }
  cmats <- lapply(hsamples,function(x){
    apply(matrix(as.vector(aperm(x$theta,c(1,3,2))),ncol=dim(x$theta)[2],
      dimnames=list(NULL,dimnames(x$theta)[[2]])),1,fun,...)
  })
  n <- dim(cmats[[1]])[1]
  av <- apply(
    array(do.call(rbind,cmats),dim=c(n,length(hsamples),dim(cmats[[1]])[2])),
    c(1,3),mean)
  dimnames(av)[[1]] <- pretty
  apply(av,1,quantile,probs=c(lo.p,.5,hi.p))
}



