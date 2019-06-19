inner<-function(u,v) (t(as.matrix(u))%*%as.matrix(v))[1]

normalise<-function(u) u %>% as.matrix %>% divide_by(norm(.,'F'))

orth<-function(V){
  V%<>%as.matrix
  ncol(V)->n
  V[,1] %>% matrix %>% norm('F') %>% c(0 %>% rep(n-1)) %>% matrix -> coord
  #1 %>% c(0 %>% rep(n-1)) %>% matrix -> coord
  V[,1]%<>%normalise
  if(n>1){
    for(i in 2:n){
      coord %<>% cbind((t(V[,1:(i-1)])%*%V[,i]) %>% c(1,0 %>% rep(n-i)))
      V[,i]%<>%subtract(V[,1:(i-1)]%*%t(t(V[,i])%*%V[,1:(i-1)]))
      coord[i,i] %<>% multiply_by(V[,i] %>% matrix %>% norm('F'))
      V[,i]%<>%normalise
    }
  }
  return(V %>% structure(coord=coord))
}

anglepp<-function(U,V){
  # U%<>%orth; V%<>%orth
  # cbind(U,V) %>% orth %@% coord %>% det %>% asin/pi
  U%<>%as.matrix; V%<>%as.matrix
  cbind(U,V) %>% orth -> base
  cbind(
    base%@%'coord' %>% extract(,1:ncol(U)) %>% orth,
    base%@%'coord' %>% extract(,ncol(U)+(1:ncol(V))) %>% orth
  ) %>% det %>% asin/pi
}

angles<-function(U,V){
  a<-NULL
  for(i in seq_len(min(ncol(U),ncol(V)))){
    a%<>%c(anglepp(U[,1:i],V[,1:i]))
  }
  return(a)
}

eigenangles<-function(data,batch,tissue){
  batch%<>%factor; batch %>% levels -> batches
  angles_batch_vs_all<-list()
  for(b in batches){
    data.all<-NULL; data.batch<-NULL
    for(t in tissue[batch==b] %>% unique){
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
      t %>% paste('\n') %>% cat
    }
    angles_batch_vs_all[[b]]<-angles(
      data.all %>% t %>% prcomp %$% rotation,
      data.batch %>% t %>% prcomp %$% rotation
    )
    b %>% paste('\n') %>% cat
  }
  return(angles_batch_vs_all)
}

#parallel version
# eigenangles<-function(data,batch,tissue){
#   batch%<>%factor; batch %>% levels -> batches
#   angles_batch_vs_all<-list()
#   cl <- detectCores() %>% subtract(1) %>% makeSOCKcluster
#   cl %>% clusterExport(c('angles','anglepp','orth','normalise','inner'))
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
#         data.all %>% t %>% prcomp %$% rotation,
#         data.batch %>% t %>% prcomp %$% rotation
#       )
#     }
#   cl %>% stopCluster
#   return(angles_batch_vs_all)
# }
