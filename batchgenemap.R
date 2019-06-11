library(magrittr)
library(ggplot2)
library(latex2exp)
library(rlang)

angle<-function(u,v,round=Inf){
  u%<>%matrix
  v%<>%matrix
  return((t(u)%*%v/(norm(u,'F')*norm(v,'F')))[1])# %>% acos/pi) #%>% divide_by(pi) %>% round(round) %>% paste('$\\pi$'))
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

anglepplook<-function(U,V,k=100){
  a<-NULL
  for(i in 1:k){
    U %>% ncol %>% runif(-1,1) -> alpha
    V %>% ncol %>% runif(-1,1) -> beta
    a%<>%c(angle(orth(U)%*%alpha,orth(V)%*%beta))
  }
  return(a %>% abs %>% max %>% acos/pi)
}

anglepplook<-function(U,V,k=100){
  a<-NULL
  for(i in 1:k){
    alpha<-runif(1,-1,1)
    beta<-sqrt(1-alpha^2)
    s<-c(-1,1) %>% sample(1)
    a%<>%c(anglelp(alpha*orth(U)[,1]+s*beta*orth(U)[,2],V))
  }
  return(a %>% abs %>% max %>% acos/pi)
}

anglelp<-function(u,V){
  angle(u,proj(u,V))
}

angleppdet<-function(U,V){
  U%<>%orth; V%<>%orth
  cbind(U,V) %>% det %>% asin/pi
}

U<-rnorm(100*3) %>% matrix(ncol=3)
V<-rnorm(100*3) %>% matrix(ncol=3)

# a<-NULL; b<-NULL
# for(i in 1:100){
#   A<-rnorm(4*2) %>% matrix(ncol=2)
#   B<-rnorm(4*2) %>% matrix(ncol=2)
#   a%<>%c(anglepplook(A,B,100))
#   b%<>%c(angleppdet(A,B))
#   cat(i)
# }
# plot(a,b %>% abs)

angles.matrix<-function(V){
  n<-ncol(V)
  angles<-NULL
  for(i in 1:n){
    for(j in 1:n){
      angles%<>%c(angle(V[,i],V[,j]))
    }
  }
  return(angles %>% matrix(n) %>% set_colnames(colnames(V)) %>% set_rownames(colnames(V)))
}

# batchgeneangles<-function(data,batch){
#   batch%<>%factor; batch %>% levels -> batches
#   data[,batch==batches[1]] %>% t %>% prcomp(rank.=1) %$% rotation -> eig1
#   data[,batch==batches[2]] %>% t %>% prcomp(rank.=1) %$% rotation -> eig2
#   data %>% t %>% prcomp(rank.=1) %$% rotation -> eig3
#   return(list(
#     inter_batch_angle=angle(eig1,eig2),
#     first_angle=angle(eig1,eig3),
#     second_angle=angle(eig2,eig3)
#   ))
# }
batchgeneangles<-function(data,batch){
  batch%<>%factor; batch %>% levels -> batches
  eig<-NULL
  for(b in batches)
    eig %<>% cbind(data[,batch==b] %>% t %>% prcomp(rank.=1) %$% rotation)
  colnames(eig)<-batches
  eig0<-data %>% t %>% prcomp(rank.=1) %$% rotation
  return(list(
    batch_vs_all = eig %>% apply(2,as_mapper(~angle(.,eig0))),
    inter_batch = angles.matrix(eig)
  ))
}
# batchgeneangles2<-function(data,batch,tissue,merge.replicates=mean,global.only=FALSE){
#   batch%<>%factor; batch %>% levels -> batches
#   angles<-NULL
#   if(!global.only){
#     for(b0 in batches){
#       for(b1 in batches){
#         common.tissues<-intersect(tissue[batch==b0],tissue[batch==b1])
#         tissue0<-tissue[batch==b0 & tissue%in%common.tissues]
#         tissue1<-tissue[batch==b1 & tissue%in%common.tissues]
#         data0<-NULL
#         for(t in tissue0){
#           if(sum(tissue0==t)>1){
#             data0%<>%cbind(data[,batch==b0&tissue==t] %>% apply(1,merge.replicates))
#           }else{
#             data0%<>%cbind(data[,batch==b0&tissue==t])
#           }
#         }
#         data1<-NULL
#         for(t in tissue1){
#           if(sum(tissue1==t)>1){
#             data1%<>%cbind(data[,batch==b1&tissue==t] %>% apply(1,merge.replicates))
#           }else{
#             data1%<>%cbind(data[,batch==b1&tissue==t])
#           }
#         }
#         angles%<>%c(angle(
#           data0 %>% t %>% prcomp(rank.=1) %$% rotation,
#           data1 %>% t %>% prcomp(rank.=1) %$% rotation
#         ))
#       }
#     }
#   }
#   common.tissues<-tissue %>% split(batch) %>% purrr::reduce(intersect) #Reduce(intersect,tissue %>% split(batch))
#   eig<-NULL
#   for(b in batches)
#     eig %<>% cbind(data[,batch==b&tissue%in%common.tissues] %>% t %>% prcomp(rank.=1) %$% rotation)
#   eig0<-data[,tissue%in%common.tissues] %>% t %>% prcomp(rank.=1) %$% rotation
#   return(list(
#     batch_vs_all=eig %>% apply(2,angle %>% partial(eig0)),
#     inter_batch=angles %>% matrix(batch %>% nlevels) %>% set_colnames(batches) %>% set_rownames(batches)
#   ))
# }

eigenangles<-function(data,batch,tissue,merge.replicates=mean){
  batch%<>%factor; batch %>% levels -> batches
  angles_batch_vs_all<-list()
  for(b in batches){
    data.all<-NULL; data.batch<-NULL
    for(t in tissue[batch==b] %>% unique){
      if(sum(tissue==t)>1){
        data.all%<>%cbind(data[,tissue==t] %>% apply(1,merge.replicates))
      }else{
        data.all%<>%cbind(data[,tissue==t])
      }
      if(sum(tissue[batch==b]==t)>1){
        data.batch%<>%cbind(data[,tissue==t & batch==b] %>% apply(1,merge.replicates))
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

batchgeneangles2<-function(data,batch,tissue,merge.replicates=mean,global.only=FALSE){
  batch%<>%factor; batch %>% levels -> batches
  angles<-NULL
  if(!global.only){
    for(b0 in batches){
      for(b1 in batches){
        common.tissues<-intersect(tissue[batch==b0],tissue[batch==b1])
        if(length(common.tissues)!=0){
          tissue0<-tissue[batch==b0 & tissue%in%common.tissues]
          tissue1<-tissue[batch==b1 & tissue%in%common.tissues]
          data0<-NULL
          for(t in tissue0){
            if(sum(tissue0==t)>1){
              data0%<>%cbind(data[,batch==b0&tissue==t] %>% apply(1,merge.replicates))
            }else{
              data0%<>%cbind(data[,batch==b0&tissue==t])
            }
          }
          data1<-NULL
          for(t in tissue1){
            if(sum(tissue1==t)>1){
              data1%<>%cbind(data[,batch==b1&tissue==t] %>% apply(1,merge.replicates))
            }else{
              data1%<>%cbind(data[,batch==b1&tissue==t])
            }
          }
          angles%<>%c(angle(
            data0 %>% t %>% prcomp(rank.=1) %$% rotation,
            data1 %>% t %>% prcomp(rank.=1) %$% rotation
          ))
        }else{
          angles%<>%c(NA)
        }
      }
    }; angles%<>%matrix(batch %>% nlevels) %>% set_colnames(batches) %>% set_rownames(batches)
  }
  angles_batch_vs_all<-NULL
  for(b in batches){
    data.all<-NULL; data.batch<-NULL
    for(t in tissue[batch==b]){
      if(sum(tissue==t)>1){
        data.all%<>%cbind(data[,tissue==t] %>% apply(1,merge.replicates))
      }else{
        data.all%<>%cbind(data[,tissue==t])
      }
      if(sum(tissue[batch==b]==t)>1){
        data.batch%<>%cbind(data[,tissue==t & batch==b]) %>% apply(1,merge.replicates)
      }else{
        data.batch%<>%cbind(data[,tissue==t & batch==b])
      }
    }
    angles_batch_vs_all%<>%c(angle(
      data.all %>% t %>% prcomp %$% rotation %>% extract(,1),
      data.batch %>% t %>% prcomp %$% rotation %>% extract(,1)
    ))
  }
  return(list(
    batch_vs_all=angles_batch_vs_all %>% set_names(batches),
    inter_batch='angles' %>% get0
  ))
}


batchgeneladder<-function(data,batch,mapping=identity,pc=1){
  batch%<>%factor; batch %>% levels -> batches
  data[,batch==batches[1]] %>% t %>% prcomp -> pca1
  data[,batch==batches[2]] %>% t %>% prcomp -> pca2
  data %>% t %>% prcomp -> pca3
  ngenes<-nrow(data)
  ggplot()+aes(
    x=c(-1,-1/3,1/3,1) %>% rep(each=ngenes),
    y=c(
      pca3$rotation[,pc],
      pca1$rotation[,pc],
      pca2$rotation[,pc],
      pca3$rotation[,pc]
    ) %>% sapply(mapping),
    colour=c('',batches[1],batches[2],'') %>% rep(each=ngenes),
    group=ngenes %>% seq_len %>% rep(4)
  )+geom_point()+scale_colour_manual(values=c('green','red','blue'),name='batch')+
    geom_line(colour='black',alpha=.2)+
    xlab('')+ylab('weight in eigengene' %>% paste(pc))+
    annotate(geom='text',x=-2/3,y=0,colour='white',label=angle(pca1$rotation[,pc],pca3$rotation[,pc]) %>% TeX(output='character'),parse=TRUE)+
    annotate(geom='text',x=0,y=0,colour='white',label=angle(pca1$rotation[,pc],pca2$rotation[,pc]) %>% TeX(output='character'),parse=TRUE)+
    annotate(geom='text',x=2/3,y=0,colour='white',label=angle(pca2$rotation[,pc],pca3$rotation[,pc]) %>% TeX(output='character'),parse=TRUE)
}

batchgenemap<-function(data,batch,pcs=1:2,xshift=NULL,yshift=NULL){
  ngenes<-nrow(data)
  batch%<>%factor
  batch %>% levels -> batches
  data[,batch==batches[1]] %>% t %>% prcomp -> pca1
  data[,batch==batches[2]] %>% t %>% prcomp -> pca2
  data %>% t %>% prcomp -> pca3
  if(xshift %>% is.null) 
    c(max(abs(pca1$rotation[,pcs[1]]))+max(abs(pca3$rotation[,pcs[1]])),0)*2->xshift
  if(yshift %>% is.null)
    c(0,max(abs(pca2$rotation[,pcs[2]]))+max(abs(pca3$rotation[,pcs[2]])))*2->yshift
  ggplot()+aes(
    x=c(
      pca1$rotation[,pcs[1]],
      pca2$rotation[,pcs[1]],
      pca3$rotation[,pcs[1]]
    )+xshift %>% c(0) %>% rep(each=ngenes),
    y=c(
      pca1$rotation[,pcs[2]],
      pca2$rotation[,pcs[2]],
      pca3$rotation[,pcs[2]]
    )+yshift %>% c(0) %>% rep(each=ngenes),
    colour=c(batches[1],batches[2],'') %>% rep(each=ngenes)
  )+geom_point()+
    scale_colour_manual(values=c('green','red','blue'),name='batch')+
    geom_line(aes(group=c(1:ngenes %>% rep(2),-(1:ngenes)) %>% factor),colour='magenta',alpha=.1)+
    geom_line(aes(group=c(1:ngenes,-(1:ngenes),1:ngenes) %>% factor),colour='orange',alpha=.1)+
    geom_line(aes(group=c(-(1:ngenes),1:ngenes %>% rep(2)) %>% factor),colour='cyan',alpha=.1)+
    xlab('coordinate on eigengene' %>% paste(pcs[1]))+ylab('coordinate on eigengene' %>% paste(pcs[2]))
}