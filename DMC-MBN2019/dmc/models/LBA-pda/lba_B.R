# Template setup for n-choice LBA, B=b-A parameterization WITH PDA
#   External parameters types: A, B, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all applications

transform.dmc <- function(par.df) 
{
  # User supplied tranforms go here
  par.df$b <- par.df$B+par.df$A
  
#   # COMMENT OUT this check for speed after debugging
#   if ( !all(type.par.names %in% names(par.df)) )
#     stop("Trasform has not created parameter(s) required by the model.")
  
  par.df[,c("A","b","t0","mean_v","sd_v","st0")]
}


random.dmc <- function(n,p.df,model)
{
  rlba.norm(n,A=p.df$A,b=p.df$b,t0=p.df$t0, 
    mean_v=p.df$mean_v,sd_v=p.df$sd_v,st0=p.df$st0[1],
    posdrift = attr(model,"posdrift"))
}




likelihood.dmc <- function(p.vector,data,min.like=1e-10)   
{

  common.cell <- function(cell.response) {
    # Gets indexes of model row.names that share cell
    cells <- unlist(lapply(
      strsplit(cell.response,".",fixed=TRUE),function(x){paste(x[-length(x)],collapse=".")}))
    list(use=!duplicated(cells),index=as.numeric(factor(cells))) 
  }
  
  likelihood <- numeric(dim(data)[1])
  ui <- common.cell(row.names(attr(data,"model")))
  conditions <- row.names(attr(data,"model"))[ui$use]
  
  for ( i in 1:length(conditions) ) { 
    
    responses <- row.names(attr(data,"model"))[ui$index==i]
    rts <- sapply(responses,function(x){
      data$RT[attr(data,"cell.index")[[x]]]  
    })
    p.df <- p.df.dmc(p.vector,conditions[i],attributes(data)$model,n1order=FALSE)

    # Get PDA sample
    samp <- try(random.dmc(attr(data,"n.pda"),p.df,attr(data,"model")),silent=TRUE)
    
    for (j in 1:length(responses))  if ( length(rts[[j]]) > 0 ) {
      if ( (class(samp)!="try-error") ) {
        is.samp <- samp$response==j
        if ( any(is.samp) ) {
         likelihood[ attr(data,"cell.index")[[ responses[j] ]] ] <- 
           pda.dmc(rts=rts[[j]], srt=samp[is.samp,"rt"], n=attr(data,"n.pda")) 
        } else likelihood[ attr(data,"cell.index")[[ responses[j]]] ] <- 0
      } else likelihood[ attr(data,"cell.index")[[ responses[j]]] ] <- 0
    }
      
  }
  pmax(likelihood,min.like)
}


