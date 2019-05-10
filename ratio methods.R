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
    data[,batch==k] %<>% subtract(data[,batch==k & reference] %>% log %>% rowMeans %>% exp)
  }
  return(data)
}

batch<-runif(100,0,2) %>% floor
B<-NULL
for(i in 1:100) B%<>%cbind(rnorm(6000,mean=.1*batch[i]))
reference<-rbinom(100,1,.1) %>% as.logical

gb<-B %>% gPCA(batch)
gb %>% viz_gpca(guided=FALSE)

ggplot(mapping=aes(x=gb$pca$u[,1],colour=batch %>% factor))+geom_density()+geom_point(y=1)

t<-seq(-.3,.3,.005)
x<-gb$pca$u[,1]
x0<-x[batch==0]
x1<-x[batch==1]
ggplot()+
  geom_line(aes(x=t,y=dnorm(t,mean(x0),sd(x0))),colour='red')+
  geom_line(aes(x=t,y=dnorm(t,mean(x1),sd(x1))),colour='blue')+
  geom_point(aes(x=x,colour=batch %>% factor),y=0)

ra<-B %>% ratioA(batch,reference)
gra<-ra %>% gPCA(batch)
gra %>% viz_gpca
gra %>% viz_gpca(guided=FALSE)
gra %>% viz_gpca_contrib


rar<-B %>% ratioA(batch, r)
grar<-rar %>% gPCA(batch)
grar %>% viz_gpca#(guided=FALSE)
ref<-reference.choice(B,batch,2)
r<-ref[,1]|ref[,2]
ggplot(mapping=aes(x=t(rar)%*%gb$pca$v[,1], y=t(rar)%*%gb$pca$v[,2],colour=batch %>% factor))+geom_point()+stat_ellipse()

reference.choice<-function(data,batch,n){
  batch %<>% factor
  batches<-batch %>% levels
  g<-data %>% gPCA(batch)
  LS<-batches %>% map(as_mapper(~g$pca$u[batch==.,1])) %>% map(as_mapper(
    ~list(mean=mean(.),sd=sd(.))
  ))
  batch.likelihoods<-LS %>% sapply(
    function(ls) g$pca$u[,1] %>% sapply(
      function(sample) sample %>% dnorm(ls$mean,ls$sd)
    )
  );colnames(batch.likelihoods)<-batches
  odd.ratios<-batch.likelihoods %>% apply(
    MARGIN=1,
    FUN=function(sample) batches %>% sapply(as_mapper(
      ~sample[.]/(sum(sample)-sample[.])
    ))
  ) %>% t;colnames(odd.ratios)<-batches
  #batch.predictions<-batch.likelihoods %>% apply(1,as_mapper(~colnames(batch.likelihoods)[which.max(.)]))
  return(odd.ratios %>% apply(
    MARGIN=2,
    FUN=function(or) TRUE %>% rep(n) %>% extract(or %>% order(decreasing=TRUE)) %>% replace(is.na(.),FALSE)
  ))
  #return(seq_along(odd.ratios) %in% c(which.max(odd.ratios),which.min(odd.ratios)))
  # densities<-batches %>% map(as_mapper(
  #   ~density(.)
  # ))
}

# reference.choice<-function(data,batch){
#   batches<-batch %>% unique
#   meansvars<-batches %>% lapply(as_mapper(~data[,batch==.])) %>% lapply(as_mapper(
#     ~list(mean=rowMeans(.), var=var(t(.)))
#   ))
#   # means<-batches %>% lapply(as_mapper(~data[,batch==.] %>% rowMeans))
#   # vars<-batches %>% lapply(as_mapper(~data[,batch==.] %>% t %>% var))
#   # batch.probs<-meansvars %>% lapply(as_mapper(
#   #   ~data %>% apply(
#   #     MARGIN=2,
#   #     FUN=function(sample) pmvnorm(sample,mean=.$mean,sigma=.$var)
#   #   )
#   # ))
#   batch.probs<-meansvars %>% lapply(
#     function(mv) data %>% apply(
#       MARGIN=2,
#       FUN=function(sample) sample %>% pmvnorm(mean=mv$mean, sigma=mv$var)
#     )
#   )
#   batch.estimate<-data %>% apply(
#     MARGIN=2,
#     FUN=(~pmvnorm)
#   )
# }