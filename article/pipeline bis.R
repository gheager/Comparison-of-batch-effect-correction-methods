library(magrittr)
library(stringr)
library(purrr)
library(ggplot2)
library(igraph)
library(gtools)

#batch effect correction packages
library(sva) #for ComBat
library(RUVSeq) #for RUVs
library(batchelor) #for mnnCorrect
library(pamr) #for batch mean centering

load_experiments<-function(directory){
  directory %>% dir %>% map(~get(load(paste0(directory,'/',.x)))$rnaseq) %>% 
    set_names(directory %>% dir %>% str_split('-') %>% map(~.x[2:3] %>% paste(collapse='')) %>% unlist)
}

download_experiments_from_ExpressionAtlas<-function(..., destdir=getwd() %>% paste('experiments',sep='/')){
  if(destdir %>% dir.exists){
    stop('Attempted to create directory `',desdir,'`` but this address already exists. Rename it or modify `destdir` argument in `download_experiments`.')
  }else{
    destdir %>% dir.creates
  }
  for(experiment in list(...)){
    paste0('https://www.ebi.ac.uk/gxa/experiments-content/',experiment,'/static/',experiment,'-atlasExperimentSummary.Rdata') %>% download.file(destdir)
  }
  desdir %>% load_experiments
}

remove_isolated_experiments<-function(experiments,factor){
  batch<-experiments %>% imap(~.y %>% rep(dim(.x)[2])) %>% unlist(use.names=FALSE)
  group<-experiments %>% map(~.[[factor]]) %>% unlist(use.names=FALSE)
  groups <- group %>% split(batch)
  intersections<-NULL
  for(i in groups){
    for(j in groups){
      intersections%<>%c(length(intersect(i,j)))
    }
  }
  intersections%<>%matrix(length(groups))%<>%set_colnames(names(groups))%<>%set_rownames(names(groups))
  intersections %>% graph_from_adjacency_matrix %>% plot
  experiments[rownames(intersections)[rowSums(intersections!=0)<=1]]<-NULL
  batch<-experiments %>% imap(~.y %>% rep(dim(.x)[2])) %>% unlist(use.names=FALSE)
  group<-experiments %>% map(~.[[factor]]) %>% unlist(use.names=FALSE)
  groups <- group %>% split(batch)
  intersections<-NULL
  for(i in groups){
    for(j in groups){
      intersections%<>%c(length(intersect(i,j)))
    }
  }
  intersections%<>%matrix(length(groups))%<>%set_colnames(names(groups))%<>%set_rownames(names(groups))
  intersections %>% graph_from_adjacency_matrix %>% plot
  return(experiments)
}

merge_experiments<-function(experiments, log=TRUE, filter=TRUE){
  genes<-experiments %>% map(rownames)
  common.genes<-genes %>% purrr::reduce(intersect)
  filtered.genes<-genes %>% map(setdiff %>% partial(y=common.genes))
  data<-experiments %>% map(~.x@assays$data$counts %>% extract(rownames(.x)%in%common.genes,)) %>% purrr::reduce(cbind)
  batch<-experiments %>% imap(~.y %>% rep(ncol(.x))) %>% unlist(use.names=FALSE) %>% factor
  if(filter){
    filter <- data %>% t %>% data.frame %>% split(batch) %>% map(~colSums(.)!=0) %>% purrr::reduce(`&`)
    data%<>%extract(filter,)
  }
  if(log) data%<>%log1p
  return(SummarizedExperiment(
    assays=if(log) list(log_counts=data) else list(counts=data),
    colData=experiments %>% map(colData) %>% smartbind(list=.) %>% set_rownames(experiments %>% map(colnames) %>% unlist) %>% cbind(batch),
    metadata=list(batch=batch)
  ))
}

correct_batch_effect<-function(experiment, model, method=c('ComBat','RUV','MNN','BMC'), k){
  log<-experiment@assays$data %>% names %>% switch(log_counts=TRUE, counts=FALSE)
  model.data<-model.frame(model, experiment@colData[all.vars(model)])
  return(SummarizedExperiment(
    assays = list(switch(
      method,
      ComBat = ComBat(experiment@assays$data[[1]], experiment$batch, mod=model.matrix(model, data=model.data)),
      RUV = RUVs(experiment@assays$data[[1]], cIdx=seq_len(nrow(experiment@assays$data[[1]])), k=k, 
                 scIdx=model.data %>% expand.grid %>% apply(1,paste) %>% makeGroups, isLog=log)$normalizedCounts,
      MNN = mnnCorrect(experiment@assays$data[[1]], batch=experiment$batch, k=k)@assays$data$corrected,
      BMC = pamr.batchadjust(list(x=experiment@assays$data[[1]],y=experiment$batch,batchlabels=experiment$batch))$x
    )) %>% set_names(if(log) 'corrected_log_counts' else 'corrected_counts'),
    colData = experiment@colData,
    metadata = experiment@metadata
  ))
}