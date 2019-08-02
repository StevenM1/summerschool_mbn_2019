# n-choice LNR, t0, mu, sigma parameterization WITH PDA
#    External parameters types: meanlog, sdlog, t0, st0 (optional)
#    Internal parameters types: meanlog, sdlog, t0, st0 (optional)

# User edited funcitons for the DMC (Dynamic Models of Choice)
#    Source in all appliciatons

# source("rtdists_extras.R")

transform.dmc <- function(par.df) 
# Transfroms parameters to a form suitbale for the model being used. 
{
  par.df[,c("meanlog","sdlog","t0","st0")]

}

random.dmc<- function(n,p.df,model)
# Retruns a data frame with columns RT and R 
{
  rlnr(n,meanlog=p.df$meanlog,sdlog=p.df$sdlog,
                      t0=p.df$t0,st0=p.df$st0[1])
}



likelihood.dmc <- function(p.vector,data,min.like=1e-10)   
{

  common.cell <- function(cell.response) {
    # Gets indexes of model row.names that share cell
    cells <- unlist(lapply(
      strsplit(cell.response,".",fixed=TRUE),function(x){x[[-length(x)]]}))
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


