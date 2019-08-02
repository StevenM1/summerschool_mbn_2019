# System functions for the DMC (Dynamic Models of Choice)
#    Usually user does not need to edit

make.level.array <- function(factors=list(S=c("s1","s2")))
  # Create array of all factorial combinations of factor levels
{  
  all.levels <- factors[[1]]
  
  if ( length(factors)>1 ) {
    for (i in 2:length(factors)) 
      all.levels <- outer(all.levels,factors[[i]],"paste",sep=".") 
    dimnames(all.levels) <- factors
  } else 
    all.levels <- array(all.levels,dim=length(factors[[1]]),dimnames=factors)
  
  all.levels
}

# constant.prior=NULL; cvs=NULL
# responses2=NULL;match.map=NULL;constants=numeric(0)
# type="norm";posdrift=TRUE;verbose=TRUE;check.model=TRUE
# 
#   factors=list(RS=c("Rr","Gr","Br"),IS=c("ri","gi","bi"))
#   responses=c("red","green","blue")
#   p.map=list(A="1",B="1",t0="1",mean_v="RI",sd_v="1",st0="1")
#   match.map=list(M=list(Rr="red",Gr="green",Br="blue"),RI=map)
#   constants=c(sd_v=1,st0=0)
  
model.dmc <- function(
  p.map,                        # list factors and constants for parameters
  responses,                    # Response (accumulator) names
  factors=list(dummy="1"),      # Factor names and levels
  cvs=NULL,                     # Names of trial covariates (in data)
  responses2=NULL,              # Second response name (multi-threshold models)
  match.map=NULL,               # Scores responses   
  constants=numeric(0),         # Parameters set to constant value
  constant.prior=NULL,          # Parameter sampled from a fixed prior
  type="norm",                  # model type
  posdrift=TRUE,                # only used by norm
  verbose=TRUE,                 # Print p.vector, constants and type
  check.model=TRUE
) 
# Creates a matrix used by get.par.mat to arrange elements of a parameter
# vector appropriately. Attributes of output used by get.par.mat to 
# add in constants and check transform.par creates the right parameters  
{
  
  grepl.dot <- function(pattern,x) {
    # Splits up pattern at ".", to get n parts. Matches each part to x and 
    # returens TRUE for each x where exactly n matches in any order.
    grepl.exact <- function(xdot,pattern) {
      xvec <- strsplit(xdot,".",fixed=TRUE)[[1]]
      any(pattern==xvec)
    } 
    
    ps <- strsplit(pattern,".",fixed=TRUE)[[1]]
    out <- sapply(x,grepl.exact,pattern=ps[1])
    if (length(ps)>1) {
      for (i in ps[-1])
        out <- rbind(out,sapply(x,grepl.exact,pattern=i))
      apply(out,2,sum)==length(ps)
    } else out
  }
  
  # Check requried inputs supplied
  if (is.null(p.map)) stop("Must supply p.map")
  if (is.null(factors)) stop("Must supply factors")
  if (is.null(responses)) stop("Must supply responses")
  
  # Check factors 
  if ( length(unlist(factors)) != length(unique(unlist(factors))) )
    stop("All factors levels must be unqiue")
  if ( any(names(factors) %in% c("1","R","R2","s")) )
    stop("Do not use s, R, R2, or 1 as a factor name")
  # Check no parameter names have a dot
  has.dot <- unlist(lapply(strsplit(names(p.map),".",fixed=TRUE),length))>1
  if ( any(has.dot) )
    stop(paste("Dots not allowed in p.map names, fix:",paste(names(p.map)[has.dot])))
  # Check R last if used
  if (any(unlist(lapply(p.map,function(x){any(x=="R") && x[length(x)]!="R"}))))
    stop("R factors must always be last")
  
  # Check responnses
  if ( type =="rd" ) {
    if (is.null(match.map))
      stop("Must specify supply a match.map for the DDM")
    if ( length(responses)!=2 ) 
      stop("DDM only applicable for two responses") 
  }
  if (!is.null(responses2)) if (!is.character(responses2) || length(responses2)<2)
    stop("responses2 must be a character vector of length 2 or greater")
  
  # Check match.map (if supplied)
  if (!is.null(match.map)) {
    # Check structure
    if ( length(match.map)<1 || class(match.map[[1]]) != "list" ) 
      stop("match.map must be a list of lists")
    # Check match.map contains at least name M
    if ( !any(names(match.map) %in% "M") ) 
      stop("match.map must have a list named M")
    map.names <- names(match.map)[names(match.map)!="M"]
    map.levels <- unlist(lapply(match.map[names(match.map)!="M"],levels))
    # convert match.map$M to responses and check
    if ( is.numeric(unlist(match.map$M)) ) 
      match.map$M <- lapply(match.map$M,function(x){responses[x]})
    if ( !all(unlist(match.map$M) %in% responses) )
      stop("match.map$M has index or name not in response names")
    if ( !(all(sort(responses)==sort(unique(unlist(match.map$M))))) )
      stop("Not all response names are scored by match.map$M")
    if ( length(match.map)>1 && 
         !all(lapply(match.map[names(match.map)!="M"],class)=="factor") )
      stop("Entries in match.map besides M must be factors")
    if ( length(unlist(map.levels)) != length(unique(unlist(map.levels))) )
      stop("All match.map levels must be unqiue")
    # Check factors
    if ( any(names(factors) == "M") )
      stop("Do not use M as a factor name")
    if ( any(names(factors) %in% names(match.map)) )
      stop(paste(map.names,"used in match.map, can not use as a factor name"))
    if ( any(unlist(factors) %in% c("true","false")) )
      stop("\"true\" and \"false\" cannot be used as factor levels")
    if ( any(map.levels %in% c("true","false")) )
      stop("\"true\" and \"false\" cannot be used as match.map levels")
    if ( length(unlist(c(factors,map.levels))) != 
         length(unique(unlist(c(factors,map.levels)))) )
      stop("Factor levels cannot overlap match.map levels")
    # Check M and R are last
    if (any(unlist(lapply(p.map,function(x){any(x=="M") && x[length(x)]!="M"}))))
      stop("M factors must always be last")
    
  }
  
  factors.short <- factors
  factors$R <- responses
  if (!is.null(match.map)) factors$M <- c("true","false")
  
  # protect againt grep problems
  for (i in unlist(factors)) if ( length(grep(i,unlist(factors)))!=1 )
    stop("Factor, response or map level is not unique or is substring of another 
         level or of \"true\" or \"false\"!" )
  
  # Add in extra match.map factors (if any)
  if ( !is.null(match.map) ) for (i in map.names)
    factors[[i]] <- levels(match.map[[i]])
  
  
  # Make parameter names
  names.par <- character(0)
  for ( i in names(p.map) ) 
  {
    if ( length(p.map[[i]])==1 && p.map[[i]] == "1" ) new.names <- i else 
    {
      new.names <- paste(i,factors[[p.map[[i]][1]]],sep=".")
      if ( length(p.map[[i]])>1 ) for ( j in 2:length(p.map[[i]]) )
        new.names <- as.vector(outer(
          new.names,factors[[p.map[[i]][j]]],"paste",sep="."
        ))
    }
    names.par <- c(names.par,new.names)
  } 
  
  # Make level array for manifest design and accumulators
  level.array <- make.level.array(factors[1:(length(factors.short)+1)])
  
  # Check match.map 
  if ( !is.null(match.map) ) for ( i in names(match.map$M) ) {
    match.position <- grep(i,level.array)
    if ( length(match.position)==0 )
      stop(paste(i,"in match.map is not in the design"))
  }
  
  # Does the par use the match factor?
  is.M <- unlist(lapply(p.map,function(x){
    any(unlist(strsplit(x,".",fixed=T)=="M"))
  }))
  
  # Does the par use a response factor 
  is.R <- unlist(lapply(p.map,function(x){
    any(unlist(strsplit(x,".",fixed=T)=="R"))
  }))
  
  if ( type =="rd"  & ( any(is.M) | any(is.R) ) ) 
    stop("Cannot use M or R in DDM p.map") 
  
  # Does the par use a map factor (i.e., in match map but not M)
  if ( !is.null(match.map) ) {
    is.map <- unlist(lapply(p.map,function(x){
      any(unlist(strsplit(x,".",fixed=T) %in% map.names))
    }))
  } else {
    is.map <- logical(length(p.map))
    names(is.map) <- names(p.map)
  }
  
  if ( any(is.map) ) {
    p.map.name <- lapply(p.map,function(x){
      unlist(strsplit(x,".",fixed=T))[
        unlist(strsplit(x,".",fixed=T)) %in% map.names]
    })
    nr <- length(responses)
    n <- length(level.array)
    map.shuffle <- matrix(aperm(array(1:n,dim=c(n/nr,nr,nr)),c(1,3,2)),ncol=nr)
  }
  
  if ( any(apply(cbind(is.M,is.R,is.map),1,function(x){sum(x)>1})) )
    stop("Parameters cannot have more than one of match.map and R factors")
  
  # use.par = boolean matrix for parameter use, cells x pars x resposnes
  use.par <- array(NA,
                   dim=c(length(level.array),length(names.par),length(responses)))
  dimnames(use.par) <- 
    list(as.vector(level.array),names.par,responses)
  
  # col.par = column parameter type (1st name)
  if ( is.null(match.map) ) 
    col.par.levels <- responses else
      col.par.levels <- c(responses,"true","false",map.levels)
  
  col.par <- strsplit(dimnames(use.par)[[2]],".",fixed=T)
  col.fac <- lapply(col.par,function(x){x[-1]})  
  col.par <- unlist(lapply(col.par,function(x){x[1]}))
  # split into fac and resp
  col.fac <- lapply(col.fac,function(x){
    if ( length(x)==0 ) out <- c(NA,NA)
    if ( length(x)==1 ) {
      if ( x %in% col.par.levels ) 
        out <- c(NA,x) else out <- c(x,NA)
    }
    if ( length(x)>1 ) 
      if ( x[length(x)] %in% col.par.levels ) 
        out <- c(paste(x[-length(x)],collapse="."),x[length(x)]) else 
          out <- paste(x,collapse=".")
        out
  })
  col.resp <- unlist(lapply(col.fac,function(x){x[2]}))
  col.fac <- unlist(lapply(col.fac,function(x){x[1]}))
  
  row.fac <- strsplit(dimnames(use.par)[[1]],".",fixed=T)  
  #  row.resp <- unlist(lapply(row.fac,function(x){x[length(x)]}))
  row.fac <- unlist(lapply(row.fac,function(x){
    paste(x[-length(x)],collapse=".")}))
  
  # Fill use.par array
  for ( p in unique(col.par) ) 
  { # parameters
    is.col <- p==col.par
    ncols <- sum(is.col)
    if ( ncols==1 ) use.par[,is.col,] <- TRUE else 
    { # there are parameter subtypes
      for ( i in 1:ncols ) 
      { # each parameter subtype
        # select rows based on factors
        tmp <- col.fac[is.col][i]
        is.fac.row <- rep(TRUE,dim(use.par)[1]) 
        if ( !is.na(tmp) ) is.fac.row[!grepl.dot(tmp,row.fac)] <- FALSE
        # set not applicable rows to false
        use.par[!is.fac.row,is.col,][,i,] <- FALSE
        if ( is.M[p] ) 
        { # has a match factor
          for ( j in names(match.map$M) ) 
          { # response cell
            correct.response <- match.map$M[[j]]
            is.rcell <- is.fac.row & grepl.dot(j,row.fac) 
            for ( k in responses )
            { # responses
              if ( k==correct.response ) 
              {
                if ( grepl("true",col.resp[is.col][i]) ) 
                  use.par[,is.col,][is.rcell,i,k] <- TRUE else
                    use.par[,is.col,][is.rcell,i,k] <- FALSE 
              } else {
                if ( grepl("false",col.resp[is.col][i]) ) 
                  use.par[,is.col,][is.rcell,i,k] <- TRUE else
                    use.par[,is.col,][is.rcell,i,k] <- FALSE                  
              }
            }
          }       
        } else if ( is.R[p] ) {
          for ( k in responses )
            use.par[is.fac.row,is.col,k][,i] <- k==col.resp[is.col][i]
        }  else if ( is.map[p] ) { 
          use.par[is.fac.row,is.col,][,i,] <- 
            match.map[[ p.map.name[[p]] ]] [map.shuffle[is.fac.row,]]==col.resp[is.col][i]
        } else use.par[is.fac.row,is.col,][,i,] <- TRUE
      }    
    }
  }
  
  if ( any(is.na(use.par)) )
    stop("Some cells of the map were not assigned!")
  
  # add in constants
  all.par <- use.par[1,,1]
  all.par[1:length(all.par)] <- NA
  if ( length(constants)>0 ) {
    if ( !all(names(constants) %in% names(all.par)) )
      stop("Name(s) in constants not in p.map")
    all.par[names(constants)] <- constants
  }
  if ( !is.null(constant.prior) ) {
    if (!all(names(constant.prior) %in% names(all.par))) 
      stop("Name(s) in constant prior not in p.map")  
    all.par[names(constant.prior)] <- rprior.dmc(constant.prior)  
  }
  attr(use.par,"all.par") <- all.par
  attr(use.par,"p.vector") <- all.par[is.na(all.par)]
  
  if (length(attr(use.par,"p.vector"))<2)
    stop("DMC cant handle models with less than two parameters")
  
  if (verbose) {
    cat("\nParameter vector names are: ( see attr(,\"p.vector\") )\n")
    print(names(all.par[is.na(all.par)]))
    cat("\nConstants are (see attr(,\"constants\") ):\n")
    print(constants)
    if (!is.null(constant.prior)) cat(paste("\nConstant prior for parameters:",
      paste(names(constant.prior),collapse=" "),"\n"))
    mod <- paste("\nModel type =",type)
    if (type=="norm") mod <- paste(mod,"(posdrift=",posdrift,")")
    cat(paste(mod,"\n\n"))
  }
  
  # parameters x cells x accumulators, used by p.list.dmc
  attr(use.par,"pca") <- aperm(use.par,c(2,1,3))
  
  # Names of parameter types (cannot have a period)
  attr(use.par,"par.names") <- unique(col.par)
  
  attr(use.par,"type") <- type
  attr(use.par,"factors") <- factors.short
  attr(use.par,"cvs") <- cvs
  attr(use.par,"responses") <- responses
  attr(use.par,"responses2") <- responses2
  attr(use.par,"constants") <- constants
  if (!is.null(constant.prior))
    attr(use.par,"constant.prior") <- constant.prior
  attr(use.par,"posdrift") <- posdrift
  
  par.df <- matrix(nrow=0,ncol=length(p.map))
  dimnames(par.df) <- list(NULL,names(p.map))
  par.df <- data.frame(par.df)
  if (!is.null(cvs))
    attr(par.df,"cvs") <- data.frame(matrix(rep(1,length(cvs)),
      nrow=1,dimnames=list(NULL,cvs)))
  par.df <- try(transform.dmc(par.df),silent=TRUE)
  if ( class(par.df)=="try-error" ) { # check list version
    par.df <- vector(mode="list",length=length(p.map))
    names(par.df) <- names(p.map)
    for (i in names(p.map)) par.df[[i]] <- matrix(nrow=0,ncol=length(responses))
    if (!is.null(cvs))
      attr(par.df,"cvs") <- data.frame(matrix(rep(1,length(cvs)),
      nrow=1,dimnames=list(NULL,cvs)))
    par.df <- try(transform.dmc(par.df),silent=TRUE)
  } 
  if ( class(par.df)=="try-error" ) { # check model matrix version
    par.df <- vector(length=length(p.map),mode="numeric")
    names(par.df) <- names(p.map)
    if (!is.null(cvs))
      attr(par.df,"cvs") <- data.frame(matrix(rep(1,length(cvs)),
      nrow=1,dimnames=list(NULL,cvs)))
    par.df <- try(transform.dmc(par.df),silent=TRUE)
  }  
  if ( class(par.df)=="try-error") 
    stop(paste("p.map must have names for each of the possible external parameter
                name (see top of model file for definitions)"))
  attr(use.par,"internal.par.names") <- names(par.df)
  
  # save "n1" orders
  resp <- unlist(lapply(strsplit(level.array,".",fixed=TRUE),function(x){
    x[length(x)]}))
  nr <- length(responses)
  n1.order <- matrix(nrow=length(resp),ncol=nr)
  for (i in 1:length(resp)) 
    n1.order[i,] <- c(c(1:nr)[resp[i]==responses],c(1:nr)[resp[i]!=responses])
  row.names(n1.order) <- row.names(use.par)
  
  # Boolean for matching cells
  match.cell <- logical(length(row.names(n1.order)))
  names(match.cell) <- row.names(n1.order)
  if ( !is.null(match.map) ) for (i in 1:length(match.map$M)) {
    match.num <- grep(match.map$M[i],names(match.cell)) 
    match.cell[match.num[match.num %in% 
                           grep(names(match.map$M[i]),names(match.cell))]] <- TRUE  
  }
  
  attr(use.par,"n1.order") <- n1.order
  attr(use.par,"match.cell") <- match.cell
  if ( !is.null(match.map) ) attr(use.par,"match.map") <- match.map
  
  if (type=="rd") # r1 cells have z flipped
  {
    is.r1 <- rep(FALSE,length(row.names(use.par)))
    names(is.r1) <- row.names(use.par) 
    is.r1[unlist(lapply(
      lapply(
        as.list(names(match.map$M)[match.map$M==responses[1]]),
        function(x)grepl(x,row.names(use.par))
      ),
      function(x) which(x==TRUE)))] <- TRUE
    attr(use.par,"is.r1") <- is.r1
    
    # add bound attributes
    bound <- rep("lower",dim(use.par)[1])
    bound[as.vector(sapply(
      paste("",names(attr(use.par,"match.map")$M),
      attr(use.par,"match.map")$M,sep="*"),
      function(x){grep(glob2rx(x),row.names(use.par))}))] <- "upper"
    names(bound) <- row.names(use.par)
    attr(use.par,"bound") <- bound
  }
  
# CHECK NOT COMPATIBLE WITH p.list.dmc ONLY p.df.dmc
#   if ( check.model ) { # Check model is good
#     p.vector <- attr(use.par,"p.vector")
#     p.vector[1:length(p.vector)] <- 1:length(p.vector)
#     if (class(try(print.cell.p(p.vector,use.par,verbose=FALSE),silent=TRUE))=="try-error")
#       stop("There is something wrong with the model that is not handled by the
#            checks implemented so far, please send a bug report to 
#            andrew.heathcote@utas.edu.au.")
#   }

  use.par
}



# Is p.vector compatible with model?
check.p.vector <- function(p.vector,model) 
{
  is.match <- names(attr(model,"p.vector")) %in% names(p.vector)
  bad <- any(!is.match)
  if ( any(!is.match) ) warning(paste("Parameter",
    names(attr(model,"p.vector"))[!is.match],"in model not present in p.vector\n"))
  is.match <- names(p.vector) %in% names(attr(model,"p.vector"))
  bad <- bad | any(!is.match)
  if ( any(!is.match) ) warning(paste("Parameter",
    names(p.vector)[!is.match],"in p.vector not present in model\n"))
  invisible(bad)
}


print.cell.p <- function(p.vector,model,verbose=TRUE) 
  # Print accumulator x internal parameter type matrix for each cell
{
  for (i in 1:dim(model)[1]) 
  {
    if (verbose) {
      print(row.names(model)[i])
      print(p.df.dmc(p.vector,i,model))
      cat("\n")
    } else {tmp <- p.df.dmc(p.vector,i,model)}
  }
}


p.df.dmc <- function(p.vector,cell,model,n1order=TRUE,cvs=NULL)
  # Gets parameter data frame (one row for each accumulator) for 
  # a design cell (specified by name or index) with model picking 
  # out the appropriate elements of par, and function transform.par
  # transforms them appropriately for model specified in model
  # Returns rows in natural (r1, r2 etc., used by simulate.dmc) or 
  # "n1" order (used by likelihood.dmc)  
{
  # Fill in constant prior (if any)
  if ( !is.null(attr(model,"constant.prior")) ) 
    attr(model,"all.par")[names(attr(model,"constant.prior"))] <- 
      rprior.dmc(attr(model,"constant.prior"))
  # Fill in non-constants
  attr(model,"all.par")[is.na(attr(model,"all.par"))] <- 
    p.vector[names(attr(model,"p.vector"))]
  # Make parameter matrix
  par.mat <- matrix(rep(attr(model,"all.par"),times=dim(model)[3])[model[cell,,]],
                    byrow=TRUE,nrow=dim(model)[3])
  dimnames(par.mat) <- list(dimnames(model)[[3]],attr(model,"par.names"))
  
  # For rd flip z
  if ( attributes(model)$type=="rd" && attributes(model)$is.r1[cell] )
    par.mat[,dimnames(par.mat)[[2]]=="z"] <- 
    1-par.mat[,dimnames(par.mat)[[2]]=="z"]
  
  if (!n1order) 
    par.df <- data.frame(par.mat) else
    par.df <- data.frame(par.mat)[attr(model,"n1.order")[cell,],]
  if (!is.null(cvs)) attr(par.df,"cvs") <- cvs
  transform.dmc(par.df)
}


p.list.dmc <- function(p.vector,model,cells,n1order=TRUE,cvs=NULL,n1.index=NULL)
  # Gets list of parameter matrices, with one row for each data point and
  # columns for each accumulator, and function transform.par
  # transforms them appropriately for model specified in model
  # Returns columns in natural (r1, r2 etc., used by simulate.dmc) or 
  # "n1" order (used by likelihood.dmc)  
{
  # Fill in constant prior (if any)
  if ( !is.null(attr(model,"constant.prior")) ) 
    attr(model,"all.par")[names(attr(model,"constant.prior"))] <- 
      rprior.dmc(attr(model,"constant.prior"))
  # Fill in non-constants
  attr(model,"all.par")[is.na(attr(model,"all.par"))] <- 
    p.vector[names(attr(model,"p.vector"))]
  
  # Make parameter array: data rows x parameters x accumulators
  parr <- array(
    rep(attr(model,"all.par"),times=dim(model)[1]*dim(model)[3])[attr(model,"pca")],
    dim=c(length(attr(model,"par.names")),dim(model)[1],dim(model)[3]),
    dimnames=list(attr(model,"par.names"),dimnames(model)[[1]],dimnames(model)[[3]])
  )[,cells,]

 
# NEEDS TO BE PORTED TO NOT REQUIRE cell
#   # For rd flip z
#   if ( attributes(model)$type=="rd" && attributes(model)$is.r1[cell] )
#     par.mat[,dimnames(par.mat)[[2]]=="z"] <- 
#     1-par.mat[,dimnames(par.mat)[[2]]=="z"]
  
  p.list <- suppressWarnings(apply(parr,1,data.frame))
  if ( !is.null(cvs) ) attr(p.list,"cvs") <- cvs
  p.list <- transform.dmc(p.list)
  
  if ( !n1order ) p.list else
    lapply(p.list,function(x){matrix(x[n1.index],ncol=dim(x)[2])})
  
  
}


# n=2; SSD=Inf;staircase=NA;TRIALS=NA; cvs=NULL
# SSD=c(Inf,Inf,.25,.25); n=6
simulate.dmc <- function(p.vector,model,n=2,cvs=NULL,SSD=Inf,staircase=NA,TRIALS=NA)
  # Create a data frame of simulated data using model
  # n can be a single number for a balanced design or a vector (one number
  #   per cell) to create an unbalanced design.
  # cvs is a data frame of covariates, must have as many rows as generated data
  # SSD is for use only with stop-signal designs, it must be a scalar, or a vector
  #   the same length as the number of cells or the same length as the data and 
  #   have Inf in all go cells
  # staircase is used only in stop-signal designs, overrides SSD, and simulates a
  #   tracking algorighm setting SSD=start initially then moving up step for each
  #   stop fail and down step for each stop-success. 
  # TRIALS is a trial number covariate used in stop-signal paradimgs to account
  #   for slowing or speeding over the course of the experiment.

# !!! TO DO !!!
# !!! Add an update ability where n=data and RT is updated,
# !!! also on first creation add index to speed update
{   
  if (check.p.vector(p.vector,model)) stop()
  # create factor data frame
  facs <- lapply(strsplit(dimnames(model)[[1]],".",fixed=TRUE),
                 function(x){x[-length(x)]})
  facs <- facs[1:(length(facs)/length(attr(model,"responses")))] 
  levs <- attr(model,"factors")
  fnams <- names(levs)
  facs <- data.frame(t(matrix(unlist(facs),nrow=length(fnams))))
  names(facs) <- fnams
  if ( !is.null(cvs) ) {
    cv.lookup <- apply(facs,1,paste,collapse=".")
    cv.out <- cvs  
  } else cv.lookup <- NULL
  
  # check n
  if ( length(n)==1 ) n <- rep(n,dim(facs)[1])
  if ( length(n)!= dim(facs)[1] )
    stop(paste("n must either be length 1 or",dim(facs)[1],"for this model."))
  
  if (!is.null(cvs) && (dim(cvs)[1]!=sum(n)) )
    stop(paste("cvs must have ",sum(n),"rows."))
  
  if ( !is.null(dimnames(n)) )
    n <- merge(facs, data.frame(n), sort=F)$Freq
  
  # create data data frame  
  data <- data.frame(lapply(facs,rep,times=n))
  for (i in names(levs)) data[[i]] <- factor(as.character(data[[i]]),levels=levs[[i]])
  data$R <- NA
  data$RT <- NA
  row1 <- 1
  
  if ( attr(model,"type") %in% c("normDK","norm2C","norm3C","norm3TC") ) # Multi-threshold race model
  {
    data$R2 <- NA
    R2.name <- "R2"
  } else R2.name <- NULL
  
  is.SSD <- FALSE  #Used to call random.dmc as a SS model or not
  # If a new model type is created it will need to be inserted here.
  
  if ( attr(model,"type") %in% c("lnrss","exgss","waldss","lbass") ) # stop signal model
  {
    is.SSD <- TRUE
    if ( any(is.na(SSD)) ) stop("SSD cannot contain NAs!")
    if ( !(length(SSD) %in% c(1,length(n),dim(data))) )
      stop(paste("SSD must have length =",dim(facs)[1],
                 "(number of cells) or length =",dim(data)[1],"(number of data points)"))
    if ( length(SSD)==1 ) {
      if ( any(is.na(staircase)) ) warning(
        paste("You have specified SSD as a scalar, unless it is Inf that\n",
              " means all cells have a stop signal, if it is Inf it means none do!")
      )
      SSD <- rep(SSD,dim(data)[1])
    }
    if ( length(SSD)==length(n) ) SSD <- rep(SSD,n)
    data$SSD <- NA
    SSD.name <- "SSD"
  } else SSD.name <- NULL
  if ( attr(model,"type") %in% c("lnrss","waldss","lbass") ) # stop signal with slowing/speeding
  {
    if ( length(TRIALS)==1 ) {
      if ( !is.na(TRIALS) ) stop("TRIALS must be NA if length = 1 (default)")
      TRIALS.name <- NULL  
    } else {
      if ( length(TRIALS) != dim(data)[1] | any(is.na(TRIALS)) ) 
        stop(paste("TRIALS cannot contain NAs and must have length =",
                   dim(data)[1],"(number of data points)"))
      TRIALS.name <- "TRIALS"
      data$TRIALS <- NA
    }
  } else TRIALS.name <- NULL
  
  
  if ( is.SSD )   # Stop signal model 
  {
    for ( i in 1:dim(facs)[1] ) if ( n[i]>0 ) 
    {
      p.df <- p.df.dmc(p.vector,i,model,n1order=FALSE)
      rown <- row1+n[i]-1
    
      data[row1:rown,c("RT","R",R2.name,SSD.name,TRIALS.name)] <- 
        random.dmc(n=n[i],p.df,model,SSD=SSD[row1:rown],staircase=staircase,
                   TRIALS=TRIALS[row1:rown])
      #Call the random function in the model file.
      #TRIALS defaults to null and is not used for exgss
      row1 <- rown+1 
    
    }
  } else {  #If not an SSD model, only n[i],p.df and model needs to be called
    for ( i in 1:dim(facs)[1] ) if ( n[i]>0 ) 
    {
      p.df <- p.df.dmc(p.vector,i,model,n1order=FALSE,
        cvs=cvs[attr(cvs,"row.facs")==cv.lookup[i],])
      rown <- row1+n[i]-1   
      data[row1:rown,c("RT","R",R2.name,SSD.name,TRIALS.name)] <- 
        random.dmc(n[i],p.df,model) #random is set in the model files.
      if ( !is.null(cvs) ) cv.out[row1:rown,] <- 
        cvs[attr(cvs,"row.facs")==cv.lookup[i],]
      row1 <- rown+1  
    }
  }
  
  data$R <- factor(data$R,levels=1:length(attr(model,"responses")),
                   labels=attr(model,"responses"))
  if ( attr(model,"type") %in% c("normDK","norm2C") ) # Multi-threshold race model
  {
    data$R2 <- factor(data$R2,levels=1:length(attr(model,"responses2")),
                      labels=attr(model,"responses2"))
  } 
  
  if ( attr(model,"type") == "rd" ) 
    # Flip responses for cells mapped to response 1
  {
    cell.names <- apply(data[,1:length(names(facs)),drop=F],1,paste,collapse=".")
    M <- attr(model,"match.map")$M
    R <- attr(model,"responses")
    for ( i in names(M) )
      if ( M[[i]] == R[1] ) 
        data[grep(i,cell.names),"R"] <- as.character(
          factor(as.character(data[grep(i,cell.names),"R"]),
                 levels=R,labels=R[2:1]))  
  }
  if ( is.null(cvs) ) data else
    cbind.data.frame(data,cv.out)
}



data.model.dmc <- function(data,model,n.pda=1e4,report=Inf)
  # Combines data and model, checking they are compatible and adding
  # a cell.index attribute to data speeding likelihood computation
{
  
  make.n1.index <- function(model,cells) 
    # make index to n1 re-order theta array   
  {
    has.na <- unlist(lapply(
      strsplit(cells,".",fixed=TRUE),function(x){x[length(x)]=="NA"}))
    nacc <- dim(attr(model,"n1.order"))[2]
    n1.order <- matrix(rep(1:nacc,each=length(cells)),
      nrow=length(cells),ncol=nacc,dimnames=list(cells,NULL))
    n1.order[!has.na,] <- attr(model,"n1.order")[cells[!has.na],]
    npar <- length(attr(model,"internal.par.names"))
    nacc <- dim(n1.order)[2]
    ndat <- dim(n1.order)[1]
    n1.index <- matrix(nrow=0,ncol=2)
    for (i in 1:nacc) 
      n1.index <- rbind(n1.index,
        cbind(1:ndat,n1.order[,i])) 
    n1.index
  }
  
  if ( is.list(model) ) { # Different model for each subject
    if (!any(names(data)=="s")) 
      stop("Only use a list of models with multiple subjects")
    if (length(model) != length(levels(data$s)))
      stop("List of models must be same length as number of subjects")
    subject.models <- TRUE
    modeli <- model[[1]]
  } else {
    subject.models <- FALSE
    modeli <- model
  }

  # check data
  if ( !is.data.frame(data) )
    stop("data must be a data frame")
  fnams <- names(attr(modeli,"factors"))
  if ( !all(c(fnams,"R","RT") %in% names(data)) )
    stop(paste("data must have columns named:",
               paste(fnams,collapse=" "),"R","RT"))
  for ( i in names(attr(modeli,"factors")) )
    if ( !all(attr(modeli,"factors")[[i]] == levels(data[,i])) ) 
      stop(paste("Factor",i,"must have these levels in this order:",
                 paste(attr(model,"factors")[[i]],collapse=" ")))
  for ( i in names(attr(modeli,"factors")) )
    if ( any(is.na(data[,i])) ) stop(paste("Factor",i,"has NAs"))
  if ( !all(sort(attr(modeli,"responses")) == sort(levels(data[,"R"]))) ) 
    stop(paste("data$R must have levels:",
               paste(attr(modeli,"responses"),collapse=" ")))
  if ( !is.numeric(data$RT) )
    stop("data$RT must be of type numeric")
  if (!is.null(attr(model,"cvs")) && 
      !all(attr(model,"cvs") %in% names(data)))
        stop("Some cvs do not have coresponding names in the data")
  

  if ( any(names(data)=="s") ) # more than one subject
  {
    dat <- data
    data <- vector(mode="list",length=length(levels(dat$s)))
    names(data) <- levels(dat$s)
    
    if (subject.models) names(model) <- levels(dat$s)
    is.sim <- !is.null(attr(dat,"parameters"))
    if (is.finite(report)) cat("Processing subjects: ")
    ss <- 0
    
    for ( s in names(data) )
    {
      if (subject.models) modeli <- model[[s]] else modeli <- model
      data[[s]] <- dat[dat$s==s,names(dat)!="s"]
      # add model and index attribute to data
      cells <- apply(data[[s]][,c(fnams,"R")],1,paste,collapse=".")
      cell.index <- vector(mode="list",length=dim(modeli)[1])
      names(cell.index) <- row.names(modeli)
      for ( i in names(cell.index) )
        cell.index[[i]] <- cells %in% i  
      attr(data[[s]],"cell.index") <- cell.index
      attr(data[[s]],"cell.empty") <- 
        unlist(lapply(cell.index,function(x){sum(x)}))==0
      attr(data[[s]],"model") <- modeli 
      if (is.sim) attr(data[[s]],"parameters") <- attr(dat,"parameters")[s,]
      attr(data[[s]],"cells") <- apply(data[[s]][,c(fnams,"R")],1,paste,collapse=".") 
      attr(data[[s]],"n1.index") <- make.n1.index(modeli,cells=attr(data[[s]],"cells"))
      
      if ( attr(modeli,"type") == "rd" ) {
        bounds <- character(dim(data[[s]])[1])
        for ( i in row.names(modeli) ) if ( !attr(data[[s]],"cell.empty")[i] )
          bounds[attr(data[[s]],"cell.index")[[i]]] <- attr(modeli,"bound")[i]
        attr(data[[s]],"bounds") <- bounds 
        attr(data[[s]],"flip") <- attr(data[[s]],"cells")  %in%
          names(attributes(modeli)$is.r1[attributes(modeli)$is.r1])
      }
      attr(data[[s]],"n.pda") <- n.pda
      ss <- ss+1
      if (is.finite(report) && (ss %% report == 0)) cat(".")
    }
    if (is.finite(report)) cat("\n")
  } else { # one subject
    # add model and index attribute to data
    cells <- apply(data[,c(fnams,"R")],1,paste,collapse=".")
    cell.index <- vector(mode="list",length=dim(model)[1])
    names(cell.index) <- row.names(model)
    for ( i in names(cell.index) )
      cell.index[[i]] <- cells %in% i  
    attr(data,"cell.index") <- cell.index
    attr(data,"cell.empty") <- 
      unlist(lapply(cell.index,function(x){sum(x)}))==0
    attr(data,"model") <- model
    attr(data,"cells") <- 
      apply(data[,c(names(attr(model,"factors")),"R")],1,paste,collapse=".")
    attr(data,"n1.index") <- make.n1.index(model,cells)
    if ( attr(model,"type") == "rd" ) {
      bounds <- character(dim(data)[1])
      for ( i in row.names(model) ) if ( !attr(data,"cell.empty")[i] )
        bounds[attr(data,"cell.index")[[i]]] <- attr(model,"bound")[i]
      attr(data,"bounds") <- bounds
      attr(data,"flip") <- attr(data,"cells")  %in%
        names(attributes(model)$is.r1[attributes(model)$is.r1])
    }
    attr(data,"n.pda") <- n.pda
  }
  data 
}


data.model.cvs.dmc <- function(p.vector,model,data,factor.names,cv.names,n.pda=1e4) 
  # Simulate using p.vector and model, and the factors and covariates in data
  # and return a data.model object.
  {
 
  if (!all(c(factor.names,cv.names) %in% names(data)))
    stop("Data does not contain some of factor and/or cv names.")
  cvs <- data[,cv.names]
  attr(cvs,"row.facs") <- apply(apply(
    data[,factor.names,drop=FALSE],2,as.character),1,paste,collapse=".")
  ns <- table(data[,factor.names],dnn=factor.names)
  data.model.dmc(simulate.dmc(p.vector,model=model,n=ns,cvs=cvs),
    model,n.pda = n.pda)
}



######### CUSTOM MAP MAKING ----

empty.map <- function(FR,levels)
{
  if (length(grep(".",levels,fixed=TRUE))!=0) 
    stop("Cannot have a dot in level names.")
  level.array <- make.level.array(FR)
  map <- rep("",length=length(level.array))
  names(map) <- level.array
  factor(map,levels=levels)
}

# include=NA;match.values=NA
# value="rm";eq.list=list(c(1,3));funs=funs[-2]
assign.map <- function(map,value="",eq.list=list(),
                       funs=NULL,include=NA,match.values=NA)
{
  
  if ( any(is.na(include)) ) ok <- rep(TRUE,length(map)) else 
  {
    ok <- grepl(include[1],names(map))
    if (length(include)>1) for (i in 2:length(include))
      ok <- ok | grepl(include[i],names(map)) 
  }
  if ( length(eq.list)==0 ) # Fill in empty elements if included
    map[is.na(map) & ok] <- value else if ( all(unlist(lapply(eq.list,length))==1) )
    {
      if (all(is.na(match.values)) || length(match.values)!=length(eq.list))
        stop("If eq.list only has length one entries match.value must contain target for each entry")
      match.mat <- matrix(unlist(lapply(eq.list,function(x){
        (unlist(lapply(strsplit(names(map),".",fixed=T),function(y){
          y[x]})))})),ncol=length(eq.list))
      map[apply(match.mat,1,function(x){all(x==match.values)})] <- value
    } else {
      if ( is.null(funs) ) funs <- lapply(eq.list,function(x){
        list("identity","identity")
      })
      if (length(funs)!=length(eq.list))
        stop("Must specify one function pair in funs for each entry in eq.list")
      if ( class(funs)!="list" || !all(unlist(lapply(funs,class))=="list") )
        stop("funs must be  list of lists")
      if ( !all(unlist(lapply(funs,length))==2) )
        stop("Each entry in funs must be length 2")
      funs <- lapply(funs,function(x){lapply(x,function(y){
        if ( is.character(y) && y=="" ) "identity" else y
      })})
      pair.list <- lapply(eq.list,function(x){
        matrix(unlist(lapply(strsplit(names(map),".",fixed=T),function(y){
          y[x]})),nrow=2)})
      map.mat <- matrix(nrow=length(map),ncol=length(eq.list))
      for ( i in 1:length(eq.list) )
        map.mat[,i] <- apply(pair.list[[i]],2,function(x){
          do.call(funs[[i]][[1]],list(x[1])) == 
            do.call(funs[[i]][[2]],list(x[2]))
        })
      map[apply(map.mat,1,all) & ok] <- value
    }
  map
}

# map <- empty.map()
# map <- assign.map( map,value="true",eq.list=list(c(1,2)) )
# map <- assign.map(map,value="false")
