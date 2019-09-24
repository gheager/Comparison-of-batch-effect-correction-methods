library(magrittr)
library(purrr)

angledet<-function(U,V){
  U%<>%as.matrix; V%<>%as.matrix
  cbind(U,V) %>% qr %>% qr.R -> new.coord
  cbind(
    new.coord[,1:ncol(U)] %>% qr %>% qr.Q,
    new.coord[,ncol(U)+(1:ncol(V))] %>% qr %>% qr.Q
  ) %>% det %>% abs %>% asin/pi
  ) %>% det %>% asin/pi
}

angles<-function(U,V){
  a<-NULL
  for(i in seq_len(min(ncol(U),ncol(V)))){
    a%<>%c(angledet(U[,1:i],V[,1:i]))
  }
  return(a)
}

eigenangles<-function(data,batch,tissue,scale=FALSE){
  data%<>%t%<>%na.omit%<>%t
eigenangles<-function(data,batch,tissue){
  batch%<>%factor; batch %>% levels -> batches
  ###########################
  angles_batch_vs_all<-list()
  for(b in batches){
    data.all<-NULL; data.batch<-NULL
    for(t in tissue[batch==b] %>% unique){
      data.all%<>%cbind(data[,tissue==t] %>% as.matrix %>% rowMeans)
      data.batch%<>%cbind(data[,tissue==t & batch==b] %>% as.matrix %>% rowMeans)
      t %>% paste('\t') %>% cat
    }
    filter<-rowVars(data.all)!=0 & rowVars(data.batch)!=0
    angles_batch_vs_all[[b]]<-angles(
      data.all[filter,] %>% t %>% prcomp(scale=scale) %$% rotation,
      data.batch[filter,] %>% t %>% prcomp(scale=scale) %$% rotation
      if(sum(tissue==t)>1){
        data.all%<>%cbind(data[,tissue==t] %>% as.matrix %>% rowMeans)
      }else{
        data.all%<>%cbind(data[,tissue==t])
      }
      if(sum(tissue[batch==b]==t)>1){
        data.batch%<>%cbind(data[,tissue==t & batch==b] %>% as.matrix %>% rowMeans)
      }else{
        data.batch%<>%cbind(data[,tissue==t & batch==b])
      }
      t %>% paste('\t') %>% cat
    }
    angles_batch_vs_all[[b]]<-angles(
      data.all %>% t %>% prcomp %$% rotation,
      data.batch %>% t %>% prcomp %$% rotation
    )
    b %>% paste('\n') %>% cat
  }
  ##########################
  angles_inter_batch<-list()
  for(i in 1:(nlevels(batch)-1)){
    for(j in (i+1):nlevels(batch)){
      data.i<-NULL; data.j<-NULL
      for(t in intersect(tissue[batch==batches[i]],tissue[batch==batches[j]])){
        data.i%<>%cbind(data[,tissue==t & batch==batches[i]] %>% as.matrix %>% rowMeans)
        data.j%<>%cbind(data[,tissue==t & batch==batches[j]] %>% as.matrix %>% rowMeans)
        t %>% paste('\t') %>% cat
      }
      if(length(intersect(tissue[batch==batches[i]],tissue[batch==batches[j]]))>ifelse(scale,1,0)){
        filter<-if(scale) rowVars(data.i)!=0 & rowVars(data.j)!=0 else TRUE
        angles_inter_batch[[batches[i]]][[batches[j]]]<-angles(
          data.i[filter,] %>% t %>% prcomp(scale=scale) %$% rotation,
          data.j[filter,] %>% t %>% prcomp(scale=scale) %$% rotation
        )
        angles_inter_batch[[batches[j]]][[batches[i]]]<-angles_inter_batch[[batches[i]]][[batches[j]]]
      if(length(intersect(tissue[batch==batches[i]],tissue[batch==batches[j]]))>0){
        angles_inter_batch[[batches[i]]][[batches[j]]]<-angles(
          data.i %>% t %>% prcomp %$% rotation,
          data.j %>% t %>% prcomp %$% rotation
        )
        #angles_inter_batch[[batches[j]]][[batches[i]]]<-angles_inter_batch[[batches[i]]][[batches[j]]]
      }
      batches[i] %>% paste('vs',batches[j],'\n') %>% cat
    }
  }
  for(b in batches) angles_inter_batch[[b]][[b]]<-0
  return(list(
    batch_vs_all=angles_batch_vs_all,
    inter_batch=angles_inter_batch
  ))
}

#parallel version
# eigenangles<-function(data,batch,tissue){
#   batch%<>%factor; batch %>% levels -> batches
#   angles_batch_vs_all<-list()
#   cl <- detectCores() %>% subtract(1) %>% makeSOCKcluster
#   cl %>% clusterExport(c('angles','angledet','orth','normalise','inner'))
#   cl %>% registerDoSNOW
#   angles_batch_vs_all<-
#     foreach(b=batches,.packages=c('magrittr','purrr','rlang')) %dopar% {
#       data.all<-NULL; data.batch<-NULL
#       for(t in tissue[batch==b] %>% unique){
#         if(sum(tissue==t)>1){
#           data.all%<>%cbind(data[,tissue==t] %>% as.matrix %>% rowMeans)
#         }else{
#           data.all%<>%cbind(data[,tissue==t])
#         }
#         if(sum(tissue[batch==b]==t)>1){
#           data.batch%<>%cbind(data[,tissue==t & batch==b] %>% as.matrix %>% rowMeans)
#         }else{
#           data.batch%<>%cbind(data[,tissue==t & batch==b])
#         }
#       }
#       angles(
#         data.all %>% t %>% prcomp(scale=scale) %$% rotation,
#         data.batch %>% t %>% prcomp(scale=scale) %$% rotation
#         data.all %>% t %>% prcomp %$% rotation,
#         data.batch %>% t %>% prcomp %$% rotation
#       )
#     }
#   cl %>% stopCluster
#   return(angles_batch_vs_all)
# }
