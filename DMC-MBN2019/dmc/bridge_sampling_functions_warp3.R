
### Written by Quentin F. Gronau

##########################################
###          helper functions          ###
##########################################

extract.data <- function(samples, type) {
  
  if (type == "not-h") {
    data <- samples$data
  } else {
    data <- lapply(samples, function(x) x$data)
  }
  
  return(data)
  
}

extract.p.prior <- function(samples, type) {
  
  if (type == "not-h") {
    p.prior <- samples$p.prior
  } else {
    # assumes that p.prior is the same for all participants
    p.prior <- samples[[1]]$p.prior
  }
  
  return(p.prior)
  
}

extract.pp.prior <- function(samples, type) {
  
  if (type == "not-h") {
    pp.prior <- NULL
  } else if (type == "h-fixed") {
    pp.prior <- NULL
  } else if (type == "h-random") {
    pp.prior <- attr(samples, "hyper")$pp.prior
  }
  
  return(pp.prior)
  
}

phi.as.mcmc.list.bridge <- function(hyper, which = 1, ...) {
  
  ### extract one of the two list elements of phi and
  ### convert it to mcmc.list object
  n.chains <- dim(hyper$phi[[which]])[1] 
  lst <- vector(mode="list",length=n.chains)
  nsamples <- dim(hyper$phi[[which]])[3]
  for (i in 1:n.chains)
    lst[[i]] <- mcmc(t(hyper$phi[[which]][i,,]),...)
  
  mcmc.list(lst)
  
}

split.samples <- function(samples, type) {
  
  # split samples into first/second half
  
  samples1 <- samples2 <- samples
  
  if (type == "not-h") {
    
    nperchain <- dim(samples$theta)[3]
    index_1 <- seq_len(nperchain) %in% seq_len(round(nperchain/2))
    samples1$theta <- samples$theta[,,index_1]
    samples2$theta <- samples$theta[,,!index_1]
    
  } else {
    
    minnperchain <- min(sapply(seq_along(samples), function(i) dim(samples[[i]]$theta)[3]))
    index_1 <- seq_len(minnperchain) %in% seq_len(round(minnperchain/2))
    
    for (i in seq_along(samples)) {
      samples1[[i]]$theta <- samples[[i]]$theta[,,index_1]
      samples2[[i]]$theta <- samples[[i]]$theta[,,!index_1]
    }
    
    if (type == "h-random") {
      attr(samples1,"hyper")$phi[[1]] <- attr(samples,"hyper")$phi[[1]][,,index_1]
      attr(samples2,"hyper")$phi[[1]] <- attr(samples,"hyper")$phi[[1]][,,!index_1]
      attr(samples1,"hyper")$phi[[2]] <- attr(samples,"hyper")$phi[[2]][,,index_1]
      attr(samples2,"hyper")$phi[[2]] <- attr(samples,"hyper")$phi[[2]][,,!index_1]
    }
  }
  
  return(list(samples1 = samples1, samples2 = samples2))
    
}

create.post.mcmc.list <- function(samples, type) {
  
  ### create mcmc.list with posterior samples
  
  if (type == "not-h") {

    post_samples_tmp <- theta.as.mcmc.list(samples)
    n.p.prior <- names(samples$p.prior)
    post_samples <- post_samples_tmp[,n.p.prior]
    
  } else {

    # create individual-level parameter mcmc.lists
    participant_samples <- vector("list", length(samples))
    for (i in seq_along(samples)) {
      
      participant_samples[[i]] <- theta.as.mcmc.list(samples[[i]])
      
      # make sure that order matches p.prior order
      n.p.prior <- names(samples[[i]]$p.prior)
      participant_samples[[i]] <- participant_samples[[i]][,n.p.prior]
      
      # add participant identifier
      varnames(participant_samples[[i]]) <- 
        paste0("part_", i, "_", varnames(participant_samples[[i]]))
      
    }
    
    # if necessary (I think only for fixed effects case), force equal number of samples
    # for each participant (necessary because we want to eventually create a matrix)
    minSize <- min(sapply(participant_samples, niter))
    participant_samples2 <- lapply(participant_samples, function(x) x[seq_len(minSize),,drop = FALSE])
    # convert back to mcmc.lists
    participant_samples <- lapply(participant_samples2, function(x) mcmc.list(lapply(x, mcmc)))

    if (type == "h-random") {

      # create group-level mcmc.lists
      hyper <- attr(samples,"hyper")
      hyperone_samples <- phi.as.mcmc.list.bridge(hyper, which = 1) 
      varnames(hyperone_samples) <- paste0("hyperone_", varnames(hyperone_samples))
      
      hypertwo_samples <- phi.as.mcmc.list.bridge(hyper, which = 2)
      varnames(hypertwo_samples) <- paste0("hypertwo_", varnames(hypertwo_samples))
      
      # combine to obtain overall mcmc.list
      post_samples <- c(participant_samples, list(hyperone_samples, hypertwo_samples))
      nchains <- length(participant_samples[[1]])
      post_samples <- lapply(seq_len(nchains),
                             function(i) {
                               mcmc(do.call(cbind, lapply(post_samples, function(x) x[[i]])))
                               })
      post_samples <- mcmc.list(post_samples)
      
    } else if (type == "h-fixed") {
      
      nchains <- length(participant_samples[[1]])
      post_samples <- lapply(seq_len(nchains),
                             function(i) {
                               mcmc(do.call(cbind, lapply(participant_samples, function(x) x[[i]])))
                             })
      post_samples <- mcmc.list(post_samples)
      
    }
  }
  
  return(post_samples)
  
}

getLbounds <- function(post_samples, samples, p.prior, type) {
  
  ### get lower bounds of parameters
  
  lb <- numeric( ncol(post_samples) )
  cn <- colnames(post_samples)
  names(lb) <- cn
  
  if (type == "not-h") {
    
    for (i in seq_along(cn)) {
      p <- cn[i]
      if (is.null(p.prior[[p]]$lower)) {
        lb[p] <- -Inf
      } else {
        lb[p] <- p.prior[[p]]$lower
      }
    }
    
  } else {
    
    nParticipants <- length(samples)
    
    # individual-level parameters
    for (i in seq_len(nParticipants)) {
      
      p.prior.i <- samples[[i]]$p.prior
      part_index <- grep(paste0("part_", i, "_"), cn)
      cn_i <- cn[part_index]
      cn_i_clean <- sapply(strsplit(x = cn_i, paste0("part_", i, "_")), function(x) x[2])
      
      for (j in seq_along(cn_i_clean)) {
        p <- cn_i_clean[j]
        if (is.null(p.prior.i[[p]]$lower)) {
          lb[cn_i[j]] <- -Inf
        } else {
          lb[cn_i[j]] <- p.prior.i[[p]]$lower
        }
      }
    }
    
    if (type == "h-random") {
      
      # group-level parameters
      hyper <- attr(samples, "hyper")
      pp.prior <- hyper$pp.prior
      
      for (k in 1:2) {
        pp.prior.k <- pp.prior[[k]]
        n <- names(pp.prior.k)
        if (k == 1) {
          post_mat_n <- paste0("hyperone_", n)
        } else if (k == 2) {
          post_mat_n <- paste0("hypertwo_", n)
        }
        
        for (l in seq_along(n)) {
          p <- n[l]
          if (is.null(pp.prior.k[[p]]$lower)) {
            lb[post_mat_n[l]] <- -Inf
          } else {
            lb[post_mat_n[l]] <- pp.prior.k[[p]]$lower
          }
        }
      }
    }
    
  }
  
  return(lb)
  
}

getUbounds <- function(post_samples, samples, p.prior, type) {
  
  ### get upper bounds of parameters
  
  ub <- numeric( ncol(post_samples) )
  cn <- colnames(post_samples)
  names(ub) <- cn
  
  if (type == "not-h") {
    
    for (i in seq_along(cn)) {
      p <- cn[i]
      if (is.null(p.prior[[p]]$upper)) {
        ub[p] <- Inf
      } else {
        ub[p] <- p.prior[[p]]$upper
      }
    }
    
  } else {
    
    nParticipants <- length(samples)
    
    # individual-level parameters
    for (i in seq_len(nParticipants)) {
      
      p.prior.i <- samples[[i]]$p.prior
      part_index <- grep(paste0("part_", i, "_"), cn)
      cn_i <- cn[part_index]
      cn_i_clean <- sapply(strsplit(x = cn_i, paste0("part_", i, "_")), function(x) x[2])
      
      for (j in seq_along(cn_i_clean)) {
        p <- cn_i_clean[j]
        if (is.null(p.prior.i[[p]]$upper)) {
          ub[cn_i[j]] <- Inf
        } else {
          ub[cn_i[j]] <- p.prior.i[[p]]$upper
        }
      }
    }
    
    if (type == "h-random") {
      
      # group-level parameters
      hyper <- attr(samples, "hyper")
      pp.prior <- hyper$pp.prior
      
      for (k in 1:2) {
        pp.prior.k <- pp.prior[[k]]
        n <- names(pp.prior.k)
        if (k == 1) {
          post_mat_n <- paste0("hyperone_", n)
        } else if (k == 2) {
          post_mat_n <- paste0("hypertwo_", n)
        }
        
        for (l in seq_along(n)) {
          p <- n[l]
          if (is.null(pp.prior.k[[p]]$upper)) {
            ub[post_mat_n[l]] <- Inf
          } else {
            ub[post_mat_n[l]] <- pp.prior.k[[p]]$upper
          }
        }
      }
    }
    
  }
  
  return(ub)
  
}

insert.constants <- function(gen_samples_tmp, theta_t, constant_index) {
  
  ### insert constants into matrix with samples from proposal distribution
  
  cn <- colnames(theta_t)
  gen_samples <- matrix(nrow = nrow(gen_samples_tmp), ncol = ncol(theta_t))
  colnames(gen_samples) <- cn
  for (i in seq_along(cn)) {
    if (constant_index[[cn[i]]]) {
      gen_samples[,i] <- rep(theta_t[[1, cn[i]]], nrow(gen_samples))
    } else {
      gen_samples[,i] <- gen_samples_tmp[,cn[i]]
    }
  }
  
  return(gen_samples)
  
}

transform2Real <- function(theta, lb, ub) {
  
  ### transform samples to real line
  
  theta_t <- theta
  transTypes <- character()
  cn <- colnames(theta)
  
  for (i in seq_len(ncol(theta))) {
    
    p <- cn[i]
    
    if (lb[[p]] < ub[[p]] && is.infinite(lb[[p]]) && is.infinite(ub[[p]])) {
      transTypes[[p]] <- "unbounded"
      theta_t[,i] <- theta[,i]
    } else if (lb[[p]] < ub[[p]] && is.finite(lb[[p]]) && is.infinite(ub[[p]])) {
      transTypes[[p]] <- "lower-bounded"
      theta_t[,i] <- log(theta[,i] - lb[[p]])
    } else if (lb[[p]] < ub[[p]] && is.infinite(lb[[p]]) && is.finite(ub[[p]])) {
      transTypes[[p]] <- "upper-bounded"
      theta_t[,i] <- log(ub[[p]] - theta[,i])
    } else if (lb[[p]] < ub[[p]] && is.finite(lb[[p]]) && is.finite(ub[[p]])) {
      transTypes[[p]] <- "double-bounded"
      theta_t[,i] <- qnorm( (theta[,i] - lb[[p]])/(ub[[p]] - lb[[p]]) )
    } else {
      stop("Could not transform parameters, possibly due to invalid lower and/or upper
           prior bounds.")
    }
    
  }
  
  colnames(theta_t) <- paste0("trans_", colnames(theta))
  
  return(list(theta_t = theta_t, transTypes = transTypes))
  
}

invtransform2Real <- function(theta_t, lb, ub) {
  
  ### transform transformed samples back to original scales
  
  theta <- theta_t
  colnames(theta) <- stringr::str_sub(colnames(theta), 7)
  cn <- colnames(theta)
  
  for (i in seq_len(ncol(theta_t))) {
    
    p <- cn[i]
    
    if (lb[[p]] < ub[[p]] && is.infinite(lb[[p]]) && is.infinite(ub[[p]])) {
      theta[,i] <- theta_t[,i]
    } else if (lb[[p]] < ub[[p]] && is.finite(lb[[p]]) && is.infinite(ub[[p]])) {
      theta[,i] <- exp(theta_t[,i]) + lb[[p]]
    } else if (lb[[p]] < ub[[p]] && is.infinite(lb[[p]]) && is.finite(ub[[p]])) {
      theta[,i] <- ub[[p]] - exp(theta_t[,i])
    } else if (lb[[p]] < ub[[p]] && is.finite(lb[[p]]) && is.finite(ub[[p]])) {
      theta[,i] <- pnorm(theta_t[,i])*(ub[[p]] - lb[[p]]) + lb[[p]]
    } else {
      stop("Could not transform parameters, possibly due to invalid lower and/or upper
           prior bounds.")
    }
    
  }
  
  return(theta)
  
}

logJacobian <- function(theta_t, transTypes, lb, ub) {
  
  ### compute log of Jacobian
  
  logJ <- matrix(nrow = nrow(theta_t), ncol = ncol(theta_t))
  cn <- stringr::str_sub(colnames(theta_t), 7)
  
  for (i in seq_len( ncol(theta_t) )) {
    
    p <- cn[i]
    
    if (transTypes[[p]] == "unbounded") {
      logJ[,i] <- 0
    } else if (transTypes[[p]] == "lower-bounded") {
      logJ[,i] <- theta_t[,i]
    } else if (transTypes[[p]] == "upper-bounded") {
      logJ[,i] <- theta_t[,i]
    } else if (transTypes[[p]] == "double-bounded") {
      logJ[,i] <- log(ub[[p]] - lb[[p]]) + dnorm(theta_t[,i], log = TRUE)
    }
    
  }
  
  return(apply(logJ, 1, sum))
  
}

log.likelihood.bridge <- function(p.vector, data, min.like = 1e-10) {
  
  ### Get log likelihood (without renaming p.vector)
  
  suppressWarnings( log(likelihood.dmc(p.vector, data, min.like = min.like) ) )
}

log.posterior.bridge <- function(p.vector, p.prior, data, min.like = 1e-10,
                                 pnames = NULL) {
  
  ### Summed log posterior likelihood
  
  if (!is.null(pnames)) names(p.vector) <- pnames
  sum (log.likelihood.bridge(p.vector, data, min.like = min.like)) + 
    summed.log.prior(p.vector, p.prior, p.names = names(p.vector))
  
}

assign.pp.bridge <- function(pp, p.prior) {
  
  ### Slot pp values into p.prior (indexed by name)
  
  for (i in names(p.prior)) 
    p.prior[[i]][1:2] <- c(pp[[1]][i],pp[[2]][i])
  p.prior
}

h.log.likelihood.bridge <- function(ps, pp, p.prior) {
  
  ### log-likelihood of subject parameters ps under population distribution p.prior
  
  suppressWarnings( apply(apply(ps, 1, log.prior.dmc,
                                p.prior = assign.pp.bridge(pp, p.prior)), 2, sum) )  
}

h.log.posterior.bridge <- function(ps, p.prior, pp, pp.prior) {
  
  ### log-likelihood of subject parameters ps under population distribution 
  ### p.prior, given population parameters pp (phi) with priors pp.prior
  ### sum over subjects of likelihood of pp (pop pars) given ps (subjects pars) 
  
  sum(h.log.likelihood.bridge(ps,pp,p.prior)) +
    # prior probability of pp
    h.summed.log.prior(pp,pp.prior) 
}

h.unnormalized.posterior <- function(s.row, data, p.prior, pp.prior,
                                     pnames = NULL, type) {
  
  if (!is.null(pnames)) names(s.row) <- pnames
  
  ### compute unnormalized posterior for a samples row
  
  if (type == "h-fixed") {
    
    n <- names(p.prior)
    nParticipants <- length(data)
    ps <- t( sapply(seq_len(nParticipants),
                    function(i) s.row[ paste0("part_", i, "_", n) ]) )
    colnames(ps) <- stringr::str_sub(colnames(ps), 8)
    
    out <-  sum( sapply(seq_len( length(data) ),
                        function(i) log.posterior.bridge(p.vector = ps[i,],
                                                         p.prior = p.prior,
                                                         data = data[[i]])) )
    
  } else if (type == "h-random") {
    
    # ps : individual-level parameters
    # pp: group-level parameters
    # pp.prior: group-level priors
    
    n <- names(pp.prior[[1]])
    p.hyperone <- s.row[ paste0("hyperone_", n) ]
    names(p.hyperone) <- stringr::str_sub(names(p.hyperone), 10)
    p.hypertwo <- s.row[ paste0("hypertwo_", n) ]
    names(p.hypertwo) <- stringr::str_sub(names(p.hypertwo), 10)
    pp <- list(p.hyperone, p.hypertwo)
    p.prior.tmp <- assign.pp.bridge(pp, p.prior)
    
    nParticipants <- length(data)
    ps <- t( sapply(seq_len(nParticipants),
                    function(i) s.row[ paste0("part_", i, "_", n) ]))
    colnames(ps) <- stringr::str_sub(colnames(ps), 8)
    
    out <-  sum( sapply(seq_len( length(data) ),
                        function(i) sum(log.likelihood.bridge(p.vector = ps[i,],
                                                              data = data[[i]])))) +
      h.log.posterior.bridge(ps = ps, p.prior = p.prior.tmp, pp = pp, pp.prior = pp.prior)
    
  }
  
  return(out)
  
}

eval.unnormalized.posterior <- function(samples_4_iter, gen_samples, data, lb, ub,
                                        p.prior, pp.prior, m, L, transTypes, type,
                                        constant_index, cores, repetitions) {
  
  ### evaluate unnormalized posterior for posterior and generated samples
  
  e <- as.brob( exp(1) )
  n_post <- nrow(samples_4_iter)
  
  twom_min_s <- matrix(2*m, nrow = n_post, ncol = length(m), byrow = TRUE) -
    samples_4_iter[,!constant_index]
  twom_min_s <- insert.constants(gen_samples_tmp = twom_min_s,
                                 theta_t = samples_4_iter,
                                 constant_index = constant_index)
  
  q21 <- m_min_gen <- m_plus_gen <- vector("list", length = repetitions)
  
  for (i in seq_len(repetitions)) {
    m_min_gen[[i]] <- matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) -
      gen_samples[[i]][,!constant_index] %*% t(L)
    m_min_gen[[i]] <- insert.constants(gen_samples_tmp = m_min_gen[[i]],
                                       theta_t = gen_samples[[i]],
                                       constant_index = constant_index)
    
    m_plus_gen[[i]] <- matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) +
      gen_samples[[i]][,!constant_index] %*% t(L)
    m_plus_gen[[i]] <- insert.constants(gen_samples_tmp = m_plus_gen[[i]],
                                        theta_t = gen_samples[[i]],
                                        constant_index = constant_index)
  }
  
  if (type == "not-h") {
    
    if (cores == 1) {
      
      q11 <- log(e^(apply(invtransform2Real(samples_4_iter, lb, ub), 1,
                          log.posterior.bridge, data = data, p.prior = p.prior) +
                      logJacobian(samples_4_iter, transTypes, lb, ub)) +
                   e^(apply(invtransform2Real(twom_min_s, lb, ub), 1,
                            log.posterior.bridge, data = data, p.prior = p.prior) +
                        logJacobian(twom_min_s, transTypes, lb, ub)))
      
      for (i in seq_len(repetitions)) {
        q21[[i]] <- log(e^(apply(invtransform2Real(m_min_gen[[i]], lb, ub), 1,
                                 log.posterior.bridge, data = data, p.prior = p.prior) +
                             logJacobian(m_min_gen[[i]], transTypes, lb, ub)) + 
                          e^(apply(invtransform2Real(m_plus_gen[[i]], lb, ub), 1,
                                   log.posterior.bridge, data = data, p.prior = p.prior) +
                               logJacobian(m_plus_gen[[i]], transTypes, lb, ub)))
      }
      
    } else if (cores > 1) {
      
      os <- get.os()
      
      if (os == "windows") {   # use a cluster
        
        require(snowfall, quietly = TRUE)
        require(rlecuyer, quietly = TRUE)
        sfInit(parallel = TRUE, cpus = cores, type = "SOCK")
        sfClusterSetupRNG()
        sfLibrary(msm)
        sfLibrary(rtdists)
        sfLibrary(pracma)
        sfLibrary(statmod)
        sfExportAll() 
        
        q11 <- log(e^(sfApply(invtransform2Real(samples_4_iter, lb, ub), 1,
                              log.posterior.bridge, data = data,
                              p.prior = p.prior) +
                        logJacobian(samples_4_iter, transTypes, lb, ub)) +
                     e^(sfApply(invtransform2Real(twom_min_s, lb, ub), 1,
                                log.posterior.bridge, data = data,
                                p.prior = p.prior) +
                          logJacobian(twom_min_s, transTypes, lb, ub))) 
        
        for (i in seq_len(repetitions)) {
          q21[[i]] <- log(e^(sfApply(invtransform2Real(m_min_gen[[i]], lb, ub), 1,
                                     log.posterior.bridge, data = data,
                                     p.prior = p.prior) +
                               logJacobian(m_min_gen[[i]], transTypes, lb, ub)) + 
                            e^(sfApply(invtransform2Real(m_plus_gen[[i]], lb, ub), 1,
                                       log.posterior.bridge, data = data,
                                       p.prior = p.prior) +
                                 logJacobian(m_plus_gen[[i]], transTypes, lb, ub)))
        }
        
        sfStop() 
        
      } else {                   # use forking
        
        require(parallel, quietly = TRUE)
        
        # q11
        pvecs.a <- data.frame(t(invtransform2Real(samples_4_iter, lb, ub)))
        pvecs.b <- data.frame(t(invtransform2Real(twom_min_s, lb, ub)))
        q11.a <- mcmapply(log.posterior.bridge, pvecs.a,
                          MoreArgs = list(data = data, p.prior = p.prior,
                                          pnames = dimnames(pvecs.a)[[1]]),
                          mc.cores = cores) +
          logJacobian(samples_4_iter, transTypes, lb, ub)
        q11.b <- mcmapply(log.posterior.bridge, pvecs.b,
                          MoreArgs = list(data = data, p.prior = p.prior,
                                          pnames = dimnames(pvecs.b)[[1]]),
                          mc.cores = cores) +
          logJacobian(twom_min_s, transTypes, lb, ub)
        q11 <- log(e^q11.a + e^q11.b)
        
        # q21
        for (i in seq_len(repetitions)) {
          pvecs.a <- data.frame(t(invtransform2Real(m_min_gen[[i]], lb, ub)))
          pvecs.b <- data.frame(t(invtransform2Real(m_plus_gen[[i]], lb, ub)))
          q21.a <- mcmapply(log.posterior.bridge, pvecs.a,
                            MoreArgs = list(data = data, p.prior = p.prior,
                                            pnames = dimnames(pvecs.a)[[1]]),
                            mc.cores = cores) +
            logJacobian(m_min_gen[[i]], transTypes, lb, ub)
          q21.b <- mcmapply(log.posterior.bridge, pvecs.b,
                            MoreArgs = list(data = data, p.prior = p.prior,
                                            pnames = dimnames(pvecs.b)[[1]]),
                            mc.cores = cores) +
            logJacobian(m_plus_gen[[i]], transTypes, lb, ub)
          q21[[i]] <- log(e^q21.a + e^q21.b)
        }
        
      }
    }
    
  } else {
    
    if (cores == 1) {
      
      q11 <- log(e^(apply(invtransform2Real(samples_4_iter, lb, ub), 1,
                          h.unnormalized.posterior, data = data,
                          p.prior = p.prior, pp.prior = pp.prior, type = type) +
                      logJacobian(samples_4_iter, transTypes, lb, ub)) +
                   e^(apply(invtransform2Real(twom_min_s, lb, ub), 1,
                            h.unnormalized.posterior, data = data,
                            p.prior = p.prior, pp.prior = pp.prior, type = type) +
                        logJacobian(twom_min_s, transTypes, lb, ub))) 
      
      for (i in seq_len(repetitions)) {
        q21[[i]] <- log(e^(apply(invtransform2Real(m_min_gen[[i]], lb, ub), 1,
                                 h.unnormalized.posterior, data = data,
                                 p.prior = p.prior, pp.prior = pp.prior, type = type) +
                             logJacobian(m_min_gen[[i]], transTypes, lb, ub)) + 
                          e^(apply(invtransform2Real(m_plus_gen[[i]], lb, ub), 1,
                                   h.unnormalized.posterior, data = data,
                                   p.prior = p.prior, pp.prior = pp.prior, type = type) +
                               logJacobian(m_plus_gen[[i]], transTypes, lb, ub)))
      }
      
    } else if (cores > 1) {
      
      os <- get.os()
      
      if (os == "windows") {   # use a cluster
        
        require(snowfall, quietly = TRUE)
        require(rlecuyer, quietly = TRUE)
        sfInit(parallel = TRUE, cpus = cores, type = "SOCK")
        sfClusterSetupRNG()
        sfLibrary(msm)
        sfLibrary(rtdists)
        sfLibrary(pracma)
        sfLibrary(statmod)
        sfExportAll() 
        
        q11 <- log(e^(sfApply(invtransform2Real(samples_4_iter, lb, ub), 1,
                              h.unnormalized.posterior, data = data,
                              p.prior = p.prior, pp.prior = pp.prior, type = type) +
                        logJacobian(samples_4_iter, transTypes, lb, ub)) +
                     e^(sfApply(invtransform2Real(twom_min_s, lb, ub), 1,
                                h.unnormalized.posterior, data = data,
                                p.prior = p.prior, pp.prior = pp.prior, type = type) +
                          logJacobian(twom_min_s, transTypes, lb, ub))) 
        
        for (i in seq_len(repetitions)) {
          q21[[i]] <- log(e^(sfApply(invtransform2Real(m_min_gen[[i]], lb, ub), 1,
                                     h.unnormalized.posterior, data = data,
                                     p.prior = p.prior, pp.prior = pp.prior, type = type) +
                               logJacobian(m_min_gen[[i]], transTypes, lb, ub)) + 
                            e^(sfApply(invtransform2Real(m_plus_gen[[i]], lb, ub), 1,
                                       h.unnormalized.posterior, data = data,
                                       p.prior = p.prior, pp.prior = pp.prior, type = type) +
                                 logJacobian(m_plus_gen[[i]], transTypes, lb, ub)))
        }
        
        sfStop() 
        
      } else {                # use forking
        
        require(parallel, quietly = TRUE)
        
        # q11
        pvecs.a <- data.frame(t(invtransform2Real(samples_4_iter, lb, ub)))
        pvecs.b <- data.frame(t(invtransform2Real(twom_min_s, lb, ub)))
        q11.a <- mcmapply(h.unnormalized.posterior, pvecs.a,
                         MoreArgs = list(data = data, p.prior = p.prior,
                                         pp.prior = pp.prior,
                                         pnames = dimnames(pvecs.a)[[1]],
                                         type = type), mc.cores = cores) +
          logJacobian(samples_4_iter, transTypes, lb, ub)
        q11.b <- mcmapply(h.unnormalized.posterior, pvecs.b,
                         MoreArgs = list(data = data, p.prior = p.prior,
                                         pp.prior = pp.prior,
                                         pnames = dimnames(pvecs.b)[[1]],
                                         type = type), mc.cores = cores) +
          logJacobian(twom_min_s, transTypes, lb, ub)
        q11 <- log(e^q11.a + e^q11.b)
        
        # q21
        for (i in seq_len(repetitions)) {
          pvecs.a <- data.frame(t(invtransform2Real(m_min_gen[[i]], lb, ub)))
          pvecs.b <- data.frame(t(invtransform2Real(m_plus_gen[[i]], lb, ub)))
          q21.a <- mcmapply(h.unnormalized.posterior, pvecs.a,
                            MoreArgs = list(data = data, p.prior = p.prior,
                                            pp.prior = pp.prior,
                                            pnames = dimnames(pvecs.a)[[1]],
                                            type = type), mc.cores = cores) +
            logJacobian(m_min_gen[[i]], transTypes, lb, ub)
          q21.b <- mcmapply(h.unnormalized.posterior, pvecs.b,
                            MoreArgs = list(data = data, p.prior = p.prior,
                                            pp.prior = pp.prior,
                                            pnames = dimnames(pvecs.b)[[1]],
                                            type = type), mc.cores = cores) +
            logJacobian(m_plus_gen[[i]], transTypes, lb, ub)
          q21[[i]] <- log(e^q21.a + e^q21.b)
        }
      }
    }
  }
  
  return(list(q11 = q11, q21 = q21))
  
}

run.iterative.scheme <- function(q11, q12, q21, q22, r0, tol,
                                 criterion, L, silent, maxiter,
                                 neff) {
  
  ### run iterative updating scheme (using "optimal" bridge function,
  ### see Meng & Wong, 1996)
  
  l1 <- -log(2) + determinant(L)$modulus + q11 - q12 # log(l)
  l2 <- -log(2) + determinant(L)$modulus + q21 - q22 # log(ltilde)
  
  lstar <- median(l1)
  n.1 <- length(l1)
  n.2 <- length(l2)
  
  if (is.null(neff)) {
    neff <- n.1
  }
  
  s1 <- neff/(neff + n.2)
  s2 <- n.2/(neff + n.2)
  r <- r0
  logml <- log(r) + lstar
  i <- 1
  
  r_vals <- r
  logml_vals <-  logml
  e <- as.brob( exp(1) )
  
  criterion_val <- 1 + tol
  
  while (criterion_val > tol && i <= maxiter) {
    
    if (! silent) {
      cat(paste0("Iteration: ", i, "\n")) 
    }
    
    rold <- r
    logmlold <- logml
    numi <- as.numeric( e^(l2 - lstar)/(s1 * e^(l2 - lstar) + s2 *  r) )
    deni <- as.numeric( 1/(s1 * e^(l1 - lstar) + s2 * r) )
    r <- (n.1/n.2) * sum(numi)/sum(deni)
    logml <- log(r) + lstar
    r_vals <- c(r_vals, r)
    logml_vals <- c(logml_vals, logml)
    i <- i + 1
    
    criterion_val <- switch(criterion, "r" = abs((r - rold)/r),
                            "logml" = abs((logml - logmlold)/logml))
    
  }
  
  if (i > maxiter) {
    logml <- NA 
  }

  return(list(logml = logml, niter = i - 1, r_vals = r_vals,
              logml_vals = logml_vals, tol = tol, neff = neff,
              criterion = criterion, maxiter = maxiter))
  
}

core.bridge.sampler.dmc <- function(samples, type, cores, repetitions,
                                    use_neff, silent, maxiter, r0,
                                    tol1, tol2) {
  
  # extract data
  data <- extract.data(samples = samples, type = type)
  
  # extract p.prior (assumes that p.prior structure is the same for all participants)
  # and extract pp.prior (if there are no group-level priors, returns NULL for pp.prior)
  p.prior <- extract.p.prior(samples = samples, type = type)
  pp.prior <- extract.pp.prior(samples = samples, type = type)
  
  # split samples in first/second half
  tmp1 <- split.samples(samples = samples, type = type)
  samples_4_fit <- tmp1$samples1
  samples_4_iter <- tmp1$samples2
  
  # convert to mcmc.lists
  samples_4_fit <- create.post.mcmc.list(samples = samples_4_fit, type = type)
  samples_4_iter <- create.post.mcmc.list(samples = samples_4_iter, type = type)
  
  # get parameter bounds
  lb <- getLbounds(post_samples = samples_4_iter[[1]], samples = samples,
                   p.prior =  p.prior, type = type)
  ub <- getUbounds(post_samples = samples_4_iter[[1]], samples = samples,
                   p.prior =  p.prior, type = type)
  
  # transform parameters to real line
  samples_4_fit <- do.call("rbind", samples_4_fit)
  tmp <- transform2Real(samples_4_fit, lb, ub)
  samples_4_fit <- tmp$theta_t
  transTypes <- tmp$transTypes
  samples_4_iter <- lapply(samples_4_iter, function(x) transform2Real(x, lb = lb, ub = ub)$theta_t)
  
  # compute effective sample size
  if (use_neff) {
    samples_4_iter <- mcmc.list(lapply(samples_4_iter, mcmc))
    neff <- tryCatch(median(coda::effectiveSize(samples_4_iter)), error = function(e) {
      warning("effective sample size cannot be calculated,
              has been replaced by number of samples.", call. = FALSE)
      return(NULL)
    })
  } else {
    neff <- NULL
  }
  
  # convert to matrix
  samples_4_iter <- do.call("rbind", samples_4_iter)

  # (if necessary) remove constants before fitting proposal
  constant_index <- apply(samples_4_fit, 2, function(x) length(unique(x)) == 1)
  samples_4_fit <- samples_4_fit[,!constant_index]
  n_post <- nrow(samples_4_iter)
  
  # get mean vector & covariance matrix and generate samples from proposal
  m <- apply(samples_4_fit, 2, mean)
  V_tmp <- cov(samples_4_fit)
  V <- as.matrix(nearPD(V_tmp)$mat) # make sure that V is positive-definite
  L <- t(chol(V))
  gen_samples_tmp <- gen_samples <- q22 <- vector("list", length = repetitions)
  
  for (i in seq_len(repetitions)) {
    gen_samples_tmp[[i]] <- rmvnorm(n_post, sigma = diag(ncol(samples_4_fit)))
    colnames(gen_samples_tmp[[i]]) <- colnames(samples_4_fit)
  }
  
  # evaluate multivariate normal distribution for posterior samples and generated samples
  q12 <- dmvnorm((samples_4_iter[,!constant_index] - matrix(m, nrow = n_post,
                                                            ncol = length(m),
                                                            byrow = TRUE)) %*%
                   t(solve(L)), sigma = diag(ncol(samples_4_fit)), log = TRUE)
  
  for (i in seq_len(repetitions)) {
    q22[[i]] <- dmvnorm(gen_samples_tmp[[i]], sigma = diag(ncol(samples_4_fit)), log = TRUE)
    
    # insert (if necessary) constants back into gen_samples matrix
    gen_samples[[i]] <- insert.constants(gen_samples_tmp = gen_samples_tmp[[i]],
                                         theta_t = samples_4_iter,
                                         constant_index = constant_index)
  }
  
  # evaluate unnormalized posterior for posterior samples and generated samples
  qList <- eval.unnormalized.posterior(samples_4_iter = samples_4_iter, gen_samples = gen_samples,
                                       data = data, lb = lb, ub = ub, p.prior = p.prior, m = m,
                                       L = L, pp.prior = pp.prior, transTypes = transTypes,
                                       constant_index = constant_index, type = type, cores = cores,
                                       repetitions = repetitions)
  q11 <- qList$q11
  q21 <- qList$q21
  
  # run iterative updating scheme to compute log of marginal likelihood
  logml <- niter <- numeric(repetitions)
  nonconverging_first <- logical(repetitions)
  complete_out <- vector("list", repetitions)
  
  for (i in seq_len(repetitions)) {
    nonconverging_first[i] <- FALSE
    tmp <- run.iterative.scheme(q11 = q11, q12 = q12, q21 = q21[[i]],
                                q22 = q22[[i]], r0 = r0, tol = tol1,
                                L = L, silent = silent,
                                maxiter = maxiter, neff = neff,
                                criterion = "r")
    
    if (is.na(tmp$logml)) {
      nonconverging_first[i] <- TRUE
      lr <- length(tmp$r_vals)
      # use geometric mean as starting value
      r0_2 <- sqrt(tmp$r_vals[lr - 1]*tmp$r_vals[lr])
      tmp <- run.iterative.scheme(q11 = q11, q12 = q12, q21 = q21[[i]],
                                  q22 = q22[[i]], r0 = r0_2, tol = tol2,
                                  L = L, silent = silent,
                                  maxiter = maxiter, neff = neff,
                                  criterion = "logml")
    }
    
    logml[i] <- tmp$logml
    niter[i] <- tmp$niter
    complete_out[[i]] <- tmp

  }
  
  return(list(logml = logml, niter = niter, nonconverging_first = nonconverging_first,
              complete_out = complete_out))
  
}


###############################################################################
### DMC users should use one of the following two functions to estimate the ###
### (log) marginal likelihood (logml) (depending on whether they fit a      ###
### hierarchical or non-hierarchical model; the previous functions are only ###
### helper functions that are called within the following two functions)    ###
###############################################################################

bridge.sampler.dmc <- function(samples, cores = 1, repetitions = 1, use_neff = TRUE,
                               silent = FALSE, maxiter = 500, r0 = 1e-5,
                               tol1 = 1e-10, tol2 = 1e-6) {
  
  ### computes the log of the marginal likelihood for a non-hierarchical model
  ### via bridge sampling
  
  core.bridge.sampler.dmc(samples, type = "not-h", cores = cores,
                          repetitions = repetitions, use_neff = use_neff,
                          silent = silent, maxiter = maxiter,
                          r0 = r0, tol1 = tol1, tol2 = tol2)
  
}

h.bridge.sampler.dmc <- function(samples, cores = 1, repetitions = 1,
                                 use_neff = TRUE, silent = FALSE, maxiter = 500,
                                 r0 = 1e-5, tol1 = 1e-10, tol2 = 1e-6) {
  
  ### computes the log of the marginal likelihood for a hierarchical model
  ### via bridge sampling
  
  has.hyper <- ! is.null(attr(samples, "hyper"))
  
  if (has.hyper) {
    # random effects
    type <- "h-random"
  } else {
    # fixed effects
    type <- "h-fixed"
  }
  
  core.bridge.sampler.dmc(samples, type = type, cores = cores,
                          repetitions = repetitions, use_neff = use_neff,
                          silent = silent, maxiter = maxiter,
                          r0 = r0, tol1 = tol1, tol2 = tol2)
  
}

###########################################
###          further functions          ###
############################################

bf.dmc <- function(bridge1, bridge2, log = FALSE) {
  
  ### computes the (log) Bayes factor in favor of bridge1 over bridge2
  
  bf <- bridge1$logml - bridge2$logml
  if(!log) {
    bf <- exp(bf)
  }
  
  return(bf)
  
}
