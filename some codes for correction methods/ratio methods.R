#Ratio-based methods
#/!\ functions ratioa and ratiog from package 'bapred' have strange behaviours : 
#I don't recommend their use
library(magrittr)
library(purrr)

ratioA<-function(data,batch,reference){
  for(k in batch %>% unique){
    data[,batch==k] %<>% subtract(data[,batch==k & reference] %>% as.matrix %>% rowMeans)
  }
  return(data)
}

ratioG<-function(data,batch,reference){
  for(k in batch %>% unique){
    data[,batch==k] %<>% subtract(data[,batch==k & reference] %>% log %>% as.matrix %>% rowMeans %>% exp)
  }
  return(data)
}

nearest.neighbours.references<-function(data,batch,k){
  batch %<>% factor
  batches<-batch %>% levels
  means<-batches %>% sapply(as_mapper(
    ~rowMeans(data[,batch==.])
  ))
  distances<-data %>% ncol %>% seq_len %>% sapply(as_mapper(
    ~(data[,.]-means[batch[.]]) %>% as.matrix %>% norm('F')
  ))
  return(FALSE %>% rep(ncol(data)) %>% replace(
    batches %>% sapply(as_mapper(
      ~replace(distances,batch!=.,Inf) %>% order %>% extract(1:k)
    )),TRUE
  ) %>% structure(distances=distances))
}
