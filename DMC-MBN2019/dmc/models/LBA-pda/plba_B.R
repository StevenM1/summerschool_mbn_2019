# Template setup for n-choice LBA, B=b-A parameterization WITH PDA
#   External parameters types: A, B, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all applications

transform.dmc <- function(par.df) 
{
  try(par.df$w <- par.df$w[2:1],silent=TRUE)
  try(par.df$sw <- par.df$sw[2:1],silent=TRUE)
  par.df[,c("A","B","C","v","sv","w","sw","rD","tD","t0","swt")]
}


random.dmc <- function(n,p.df,model)
{
  rplba(n,A=p.df$A,B=p.df$B,C=p.df$C,
          v=p.df$v, w=p.df$w, sv=p.df$sv,sw=p.df$sw,
          rD=p.df$rD[1], tD=p.df$tD[1], t0=p.df$t0[1], swt=p.df$swt[1])
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
    samp <- try(random.dmc(attr(data,"n.pda"),p.df),silent=TRUE)

    for (j in 1:length(responses))  if ( length(rts[[j]]) > 0 ) {
      if ( (class(samp)!="try-error") ) {
        is.samp <- samp[,"choice"]==j
        if ( any(is.samp) ) {
         likelihood[ attr(data,"cell.index")[[ responses[j] ]] ] <-
           pda.dmc(rts=rts[[j]], srt=samp[is.samp,"rt"], n=attr(data,"n.pda"))
        } else likelihood[ attr(data,"cell.index")[[ responses[j]]] ] <- 0
      } else likelihood[ attr(data,"cell.index")[[ responses[j]]] ] <- 0
    }

  }
  pmax(likelihood,min.like)
}



