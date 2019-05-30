library(magrittr)
library(pamr)
library(sva)
library(RUVSeq)
library(harmony)

correct.batch.effect<-function(data,batch,
                               method=c('none','bmc','combat','harmony',
                                        'ruvg','ruvs','ruvr'),
                               model,log=TRUE,
                               gPCA=TRUE,genemap=TRUE){
  if(log) data%<>%log1p
  if(method=='none') return(data)
  else if(method=='bmc'){
    if(missing(model)) return(pamr.batchadjust(list(x=data,batchlabels=batch))$x)
  }
  else if(method=='combat'){
    if(missing(model)) return(ComBat(data,batch))
    else return(ComBat(data,batch,mod=model.matrix(model)))
  }
  else if(method=='ruvg'){
    if(missing(model)) 
  }
}