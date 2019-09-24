#mixing entropy
library(purrr)
library(magrittr)
library(ggplot2)

mixing_entropy<-function(experiment){
  experiment@assays$data[[1]] %>% t %>% dist %>% as.matrix %>% 
    apply(2,list) %>% map(purrr::simplify) %>% 
    map(order) -> o
  ncol(experiment) %>% seq_len %>% map(
    function(k){
      o %>% map(~.x[1:k]) %>% map(~experiment@metadata$batch[.x]) %>% map(table) %>% map(~.x/k) %>% 
        map(~-sum(ifelse(.x!=0,.x*log(.x),0))) %>% reduce(sum)
    }
  ) %>% unlist
}

get(load('Mouse Expression Atlas/integrated data/uncorrected.Rdata')) %>% mixing_entropy %>% 
  save(file='Mouse Expression Atlas/entropy/uncorrected.Rdata')
get(load('Mouse Expression Atlas/integrated data/ComBat_corrected.Rdata')) %>% mixing_entropy %>% 
  save(file='Mouse Expression Atlas/entropy/ComBat.Rdata')
get(load('Mouse Expression Atlas/integrated data/BMC_corrected.Rdata')) %>% mixing_entropy %>% 
  save(file='Mouse Expression Atlas/entropy/BMC.Rdata')
get(load('Mouse Expression Atlas/integrated data/RUVs_corrected_k=1:20.Rdata')) %>% map(mixing_entropy) %>% 
  save(file='Mouse Expression Atlas/entropy/RUVs_k=1:20.Rdata')
get(load('Mouse Expression Atlas/integrated data/MNN_corrected_k=1:20.Rdata')) %>% map(mixing_entropy) %>% 
  save(file='Mouse Expression Atlas/entropy/MNN_k=1:20.Rdata')

KNN=119;K=20
ggplot(data.frame(
  kNN=c(
    1:KNN %>% rep(K),
    1:KNN %>% rep(K),
    1:KNN,
    1:KNN,
    1:KNN
  ),
  entropy=c(
    get(load('Mouse Expression Atlas/entropy/RUVs_k=1:20.Rdata')) %>% unlist,
    get(load('Mouse Expression Atlas/entropy/MNN_k=1:20.Rdata')) %>% unlist,
    get(load('Mouse Expression Atlas/entropy/ComBat.Rdata')),
    get(load('Mouse Expression Atlas/entropy/BMC.Rdata')),
    get(load('Mouse Expression Atlas/entropy/uncorrected.Rdata'))
  ),
  algorithm=c(
    'RUVs' %>% rep(K*KNN),
    'MNN' %>% rep(K*KNN),
    'ComBat' %>% rep(KNN),
    'BMC' %>% rep(KNN),
    'uncorrected' %>% rep(KNN)
  ),
  k=c(
    1:K %>% rep(each=KNN),
    1:K %>% rep(each=KNN),
    NA %>% rep(KNN),
    NA %>% rep(KNN),
    NA %>% rep(KNN)
  )
))+aes(x=kNN,y=entropy,colour=algorithm,alpha=k,group=paste(algorithm,k))+geom_line()+
  scale_colour_manual(values=c('orange','red','green','blue','black'))

get(load('Human Expression Atlas/integrated data/uncorrected.Rdata')) %>% mixing_entropy %>% 
  save(file='Human Expression Atlas/entropy/uncorrected.Rdata')
get(load('Human Expression Atlas/integrated data/ComBat_corrected.Rdata')) %>% mixing_entropy %>% 
  save(file='Human Expression Atlas/entropy/ComBat.Rdata')
get(load('Human Expression Atlas/integrated data/BMC_corrected.Rdata')) %>% mixing_entropy %>% 
  save(file='Human Expression Atlas/entropy/BMC.Rdata')
get(load('Human Expression Atlas/integrated data/RUVs_corrected_k=1:5.Rdata')) %>% map(mixing_entropy) %>% 
  save(file='Human Expression Atlas/entropy/RUVs_k=1:5.Rdata')
get(load('Human Expression Atlas/integrated data/MNN_corrected_k=1:5.Rdata')) %>% map(mixing_entropy) %>% 
  save(file='Human Expression Atlas/entropy/MNN_k=1:5.Rdata')

KNN=317;K=5
ggplot(data.frame(
  kNN=c(
    1:KNN %>% rep(K),
    1:KNN %>% rep(K),
    1:KNN,
    1:KNN,
    1:KNN
  ),
  entropy=c(
    get(load('Human Expression Atlas/entropy/RUVs_k=1:5.Rdata')) %>% unlist,
    get(load('Human Expression Atlas/entropy/MNN_k=1:5.Rdata')) %>% unlist,
    get(load('Human Expression Atlas/entropy/ComBat.Rdata')),
    get(load('Human Expression Atlas/entropy/BMC.Rdata')),
    get(load('Human Expression Atlas/entropy/uncorrected.Rdata'))
  ),
  algorithm=c(
    'RUVs' %>% rep(K*KNN),
    'MNN' %>% rep(K*KNN),
    'ComBat' %>% rep(KNN),
    'BMC' %>% rep(KNN),
    'uncorrected' %>% rep(KNN)
  ),
  k=c(
    1:K %>% rep(each=KNN),
    1:K %>% rep(each=KNN),
    NA %>% rep(KNN),
    NA %>% rep(KNN),
    NA %>% rep(KNN)
  )
))+aes(x=kNN,y=entropy,colour=algorithm,alpha=k,group=paste(algorithm,k))+geom_line()+
  scale_colour_manual(values=c('orange','red','green','blue','black'))

# mixing_entropy<-function(experiment, k=NaN){
#   if(k %>% is.finite){
#     experiment@assays$data$corrected %>% t %>% dist %>% as.matrix %>% 
#       apply(2,list) %>% map(purrr::simplify) %>% 
#       map(order) %>% map(~.[1:k]) %>% map(~experiment@metadata$batch[.]) %>% map(table) %>% map(~./k) %>% 
#       map(~-sum(ifelse(.!=0,.*log(.),0))) %>% unlist %>% mean
#   }else{
#     experiment@assays$data$corrected %>% t %>% dist %>% as.matrix %>% 
#       apply(2,list) %>% map(purrr::simplify) %>% 
#       map(order) -> o
#     l<-NULL
#     for(k in 1:dim(experiment)[2]){
#       l%<>%c(
#         o %>% map(~.[1:k]) %>% map(~experiment@metadata$batch[.]) %>% map(table) %>% map(~./k) %>% 
#           map(~-sum(ifelse(.!=0,.*log(.),0))) %>% unlist %>% mean
#       )
#     }
#     l
#   }
# }
# 
# 
# mixing_entropy2<-function(experiment, k=NaN){
#   if(k %>% is.finite){
#     experiment@assays$data$corrected %>% t %>% dist %>% as.matrix %>% 
#       apply(2,list) %>% map(purrr::simplify) %>% 
#       map(order) %>% map(~.[1:k]) %>% map(~experiment@metadata$batch[.]) %>% map(table) %>% map(~./k) %>% 
#       reduce(rbind) %>% colMeans %>% {-sum(ifelse(.!=0,.*log(.),0))}
#   }else{
#     experiment@assays$data$corrected %>% t %>% dist %>% as.matrix %>% 
#       apply(2,list) %>% map(purrr::simplify) %>% 
#       map(order) -> o
#     l<-NULL
#     for(k in 2:dim(experiment)[2]){
#       l%<>%c(
#         o %>% map(~.[2:k]) %>% map(~experiment@metadata$batch[.]) %>% map(table) %>% map(~./k) %>% 
#           reduce(rbind) %>% colMeans %>% {-sum(ifelse(.!=0,.*log(.),0))}
#       )
#     }
#     l
#   }
# }


# get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas no correction.Rdata')) %>% mixing_entropy -> none_entropy
# get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas ComBat.Rdata')) %>% mixing_entropy -> combat_entropy
# get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas RUV k=1:10.Rdata')) %>% map(mixing_entropy) -> ruv_entropy
# get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas MNN k=1:10.Rdata')) %>% map(mixing_entropy) -> mnn_entropy
# 
# get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas no correction.Rdata')) %>% mixing_entropy2 -> none_entropy
# get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas ComBat.Rdata')) %>% mixing_entropy2 -> combat_entropy
# get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas RUV k=1:10.Rdata')) %>% map(mixing_entropy2) -> ruv_entropy
# get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas MNN k=1:10.Rdata')) %>% map(mixing_entropy2) -> mnn_entropy
# 
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
  
