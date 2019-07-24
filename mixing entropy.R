#mixing entropy
library(purrr)
library(magrittr)
library(ggplot2)

mixing_entropy<-function(experiment, k=NaN){
  if(k %>% is.finite){
    experiment@assays$data$corrected %>% t %>% dist %>% as.matrix %>% 
      apply(2,list) %>% map(purrr::simplify) %>% 
      map(order) %>% map(~.[1:k]) %>% map(~experiment@metadata$batch[.]) %>% map(table) %>% map(~./k) %>% 
      map(~-sum(ifelse(.!=0,.*log(.),0))) %>% unlist %>% mean
  }else{
    experiment@assays$data$corrected %>% t %>% dist %>% as.matrix %>% 
      apply(2,list) %>% map(purrr::simplify) %>% 
      map(order) -> o
    l<-NULL
    for(k in 1:dim(experiment)[2]){
      l%<>%c(
        o %>% map(~.[1:k]) %>% map(~experiment@metadata$batch[.]) %>% map(table) %>% map(~./k) %>% 
          map(~-sum(ifelse(.!=0,.*log(.),0))) %>% unlist %>% mean
      )
    }
    l
  }
}


mixing_entropy2<-function(experiment, k=NaN){
  if(k %>% is.finite){
    experiment@assays$data$corrected %>% t %>% dist %>% as.matrix %>% 
      apply(2,list) %>% map(purrr::simplify) %>% 
      map(order) %>% map(~.[1:k]) %>% map(~experiment@metadata$batch[.]) %>% map(table) %>% map(~./k) %>% 
      reduce(rbind) %>% colMeans %>% {-sum(ifelse(.!=0,.*log(.),0))}
  }else{
    experiment@assays$data$corrected %>% t %>% dist %>% as.matrix %>% 
      apply(2,list) %>% map(purrr::simplify) %>% 
      map(order) -> o
    l<-NULL
    for(k in 2:dim(experiment)[2]){
      l%<>%c(
        o %>% map(~.[2:k]) %>% map(~experiment@metadata$batch[.]) %>% map(table) %>% map(~./k) %>% 
          reduce(rbind) %>% colMeans %>% {-sum(ifelse(.!=0,.*log(.),0))}
      )
    }
    l
  }
}

get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas no correction.Rdata')) %>% mixing_entropy -> none_entropy
get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas ComBat.Rdata')) %>% mixing_entropy -> combat_entropy
get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas RUV k=1:10.Rdata')) %>% map(mixing_entropy) -> ruv_entropy
get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas MNN k=1:10.Rdata')) %>% map(mixing_entropy) -> mnn_entropy

get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas no correction.Rdata')) %>% mixing_entropy2 -> none_entropy
get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas ComBat.Rdata')) %>% mixing_entropy2 -> combat_entropy
get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas RUV k=1:10.Rdata')) %>% map(mixing_entropy2) -> ruv_entropy
get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas MNN k=1:10.Rdata')) %>% map(mixing_entropy2) -> mnn_entropy

# entropy.df<-data.frame(matrix(nrow=0,ncol=3),stringsAsFactors = FALSE) %>% set_colnames(c('method','k','entropy'))
# entropy.df%<>%rbind(cbind(method='none',k=0,entropy=none_entropy))
# entropy.df%<>%rbind(cbind(method='ComBat',k=0,entropy=combat_entropy))
# entropy.df%<>%rbind(cbind(method='RUV', ruv_entropy %>% imap(~data.frame(k=.y,entropy=.x)) %>% purrr::reduce(rbind)))
# entropy.df<-rbind(
#   cbind(method='RUV', ruv_entropy %>% imap(~data.frame(k=.y,entropy=.x)) %>% purrr::reduce(rbind)),
#   cbind(method='MNN', mnn_entropy %>% imap(~data.frame(k=.y,entropy=.x)) %>% purrr::reduce(rbind))
# )
ggplot()+
  geom_text(data=ruv_entropy %>% imap(~data.frame(k=.y,entropy=.x)) %>% purrr::reduce(rbind),
            mapping=aes(x=seq_along(none_entropy) %>% rep(10), y=entropy, label=k),colour='skyblue')+
  geom_text(data=mnn_entropy %>% imap(~data.frame(k=.y,entropy=.x)) %>% purrr::reduce(rbind),
            mapping=aes(x=seq_along(none_entropy) %>% rep(10), y=entropy, label=k),colour='green')+
  geom_line(aes(x=seq_along(none_entropy),y=none_entropy),colour='black')+
  geom_line(aes(x=seq_along(combat_entropy),y=combat_entropy),colour='red')#+
  xlim(c(0,20))
  
