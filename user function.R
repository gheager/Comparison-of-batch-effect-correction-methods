library(magrittr)
library(pamr)
library(sva)
library(RUVSeq)
library(batchelor)
#library(harmony)

library(SummarizedExperiment)

correct.batch.effect<-function(data,batch,
                               method=c('none','combat','ruv','mnn','bmc'),
                               model,log=TRUE,model.data,k=1){
  for(u in model.data %>% seq_along){
    eval(parse(text=paste0(
      names(model.data)[u],'<-model.data[[u]]'
    )))
  }
  if(log) data%<>%log1p
  if(method=='none') return(data)
  else if(method=='bmc'){
    if(missing(model)) return(pamr.batchadjust(list(x=data,batchlabels=batch))$x)
  }
  else if(method=='combat'){
    if(missing(model)) return(ComBat(data,batch))
    else return(ComBat(data,batch,mod=model.matrix(model,data=model.data)))
  }
  else if(method=='ruv'){
    if(missing(model)) stop('Biological model is needed to run RUVs')
    else return(RUVs(data,cIdx=seq_len(nrow(data)),k=k,
                     scIdx=makeGroups(expand.grid(model.data) %>% apply(1,paste)),isLog=log)$normalizedCounts)
  }
  else if(method=='mnn'){
    return(mnnCorrect(data,batch=batch,k=k) %>% assays %$% corrected)
  }
}

remove.batch.effect<-function(...,list=NULL,model=NULL,method=c('none','combat','ruv','bmc','mnn'),log=TRUE,k=1){
  if(is.null(list)) experiments<-list(...) else experiments<-list; if(is.null(names(experiments))) names(experiments)<-paste0('batch',experiments %>% seq_along)
  genes<-experiments %>% map(rownames)
  common.genes<-genes %>% purrr::reduce(intersect)
  filtered.genes<-genes %>% map(setdiff %>% partial(y=common.genes))
  data<-NULL;batch<-NULL;vars<-list()
  for(i in experiments %>% seq_along){
    experiments[[i]]->exp
    data%<>%cbind(exp %>% assays %$% counts %>% extract(rownames(exp)%in%common.genes,))
    batch%<>%c(names(experiments)[[i]] %>% rep(dim(exp)[2]))
    for(v in all.vars(model)) vars[[v]]%<>%c(exp[[v]])
  }; batch%<>%factor
  filter <- data %>% t %>% data.frame %>% split(batch) %>% map(~colSums(.)!=0) %>% purrr::reduce(`&`)
  data%<>%extract(filter,)
  data%<>%correct.batch.effect(batch,method,model,log,model.data=model.frame(model,vars),k=k)
  return(SummarizedExperiment(
    assays=list(counts=data),
    colData=gtools::smartbind(list=experiments %>% map(colData)) %>% set_rownames(experiments %>% map(colnames) %>% purrr::reduce(c)),
    metadata=list(batch=batch)
  ))
}

eigenangles.summaryexperiment<-function(experiment,scale=FALSE){
  eigenangles(
    experiment %>% assays %$% counts,
    experiment %>% metadata %$% batch,
    experiment$organism_part,
    scale=scale
  )
}

gPCA.integrate<-function(experiment,group){
  experiment %>% assays %$% counts %>% gPCA(experiment[[group]],scaleY=TRUE)
}