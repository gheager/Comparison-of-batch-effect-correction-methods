library(magrittr)
library(ExpressionAtlas)

merge.replicates<-function(experimentSummary,group.column,merging){
  if(!is.null(experimentSummary$rnaseq[[group.column]])){
    experimentSummary$rnaseq %>% assays %>% use_series(counts) -> data
    experimentSummary$rnaseq[[group.column]] %>% factor -> group
    merged<-NULL
    pheno<-NULL
    for(g in group %>% levels){
      cat(g)
      if(sum(group==g)>1){
        merged %<>% cbind(
          data[,group==g] %>% apply(1,merging)
        )
      }else{
        merged %<>% cbind(
          data[,group==g]
        )
      }
      if(all(experimentSummary$rnaseq %>% colData %>% extract(group==g,) %>% apply(2,as_mapper(
        ~length(unique(.))==1
      )))){
        pheno%<>%rbind(experimentSummary$rnaseq %>% colData %>% extract(group==g,) %>% apply(2,unique))
      }else{
        pheno%<>%rbind(experimentSummary$rnaseq %>% colData %>% extract(group==g,) %>% extract(1,) %>% data.frame)
        warning('Phenological data are not identical among the group :' %>% paste(g))
      }
    }
    return(newSeqExpressionSet(
      merged,
      phenoData=pheno %>% data.frame
    ))
  }else{
    message('No' %>% paste(group.column))
    return(experimentSummary$rnaseq)
  }
}

angle<-function(u,v){
  u%<>%matrix
  v%<>%matrix
  return((t(u)%*%v/(norm(u,'F')*norm(v,'F')))[1])
}

angles<-function(data,tissue,batch){
  batch%<>%factor
  levels(batch)<-c(0,1)
  common.tissues<-tissue[batch==0] %>% intersect(tissue[batch==1])
  arrows<-NULL
  for(o in common.tissues){
    arrows%<>%cbind(data[,tissue==o&batch==1]-data[,tissue==o&batch==0])
  };colnames(arrows)<-common.tissues
  arrows %>% apply(2,
    function(u) arrows %>% apply(2,
      function(v) angle(u,v)
    )
  )
}

tissue.distances<-function(data,tissue,batch){
  batch%<>%factor
  levels(batch)<-c(0,1)
  common.tissues<-tissue[batch==0] %>% intersect(tissue[batch==1])
  arrows<-NULL
  for(o in common.tissues){
    arrows%<>%cbind(data[,tissue==o&batch==1]-data[,tissue==o&batch==0])
  };colnames(arrows)<-common.tissues
  arrows %>% apply(2,as_mapper(
    ~norm(matrix(.),'F')
  ))
}