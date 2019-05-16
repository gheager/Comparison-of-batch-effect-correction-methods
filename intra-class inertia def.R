library(magrittr)
library(purrr)
library(matrixStats)
batch_inertia <- function(data, outcome){
  outcome%<>%factor
  vars<-NULL
  for(x in outcome %>% levels){
    if(sum(outcome==x)>1){
      vars%<>%c(rowVars(data[,which(outcome==x)],na.rm=TRUE) %>% mean %>% sqrt)
    }
  }
  return(vars %>% mean)
}

batch.inertia<-function(data,group,doPCA=TRUE,transpose=doPCA,center=TRUE,scale=FALSE){
  group%<>%factor
  group %>% levels -> groups
  if(transpose) data %<>% t
  if(doPCA) data %<>% scale(center,scale) %<>% na.omit %<>% svd %<>% use_series(u)
  cumvars<-matrix(nrow=length(groups),ncol=ncol(data)) %>% set_rownames(groups)
  for(g in groups){
    if(sum(group==g)>1){
      cumvars[g,]<-data %>% apply(2,as_mapper(
        ~var(.[group==g])
      ))
    }
  }
  cumvars %<>% na.omit
  return(list(
    cumvars=cumvars,
    total=cumvars %>% colMeans %>% cumsum,
    volume.index=cumvars %>% sqrt %>% apply(1,cumprod) %>% apply(1,mean)
  ))
}
