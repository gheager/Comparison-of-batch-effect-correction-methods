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

batch<-runif(100,0,2) %>% floor
B<-NULL
for(i in 1:100) B%<>%cbind(rnorm(6000,mean=.1*batch[i]))
reference<-rbinom(100,1,.1) %>% as.logical

gb<-B %>% gPCA(batch)
gb %>% viz_gpca(guided=FALSE)
gb %>% viz_gpca

ggplot(mapping=aes(x=gb$pca$u[,1],colour=batch %>% factor))+geom_density()+geom_point(y=1)

t<-seq(-.3,.3,.005)
x<-gb$pca$u[,1]
x0<-x[batch==0]
x1<-x[batch==1]
ggplot()+
  geom_line(aes(x=t,y=dnorm(t,mean(x0),sd(x0))),colour='red')+
  geom_line(aes(x=t,y=dnorm(t,mean(x1),sd(x1))),colour='blue')+
  geom_point(aes(x=x,colour=batch %>% factor),y=0)

#reference<-TRUE %>% rep(2) %>% extract(1:100 %>% sample) %>% replace(is.na(.),FALSE)
batch<-runif(100,0,2) %>% floor
B<-NULL
for(i in 1:100) B%<>%cbind(rnorm(6000,mean=.1*batch[i]))
reference<-nearest.neighbours.references(B,batch,40)
batch[which(reference)]
ra<-B %>% ratioA(batch,reference)
gra<-ra %>% gPCA(batch)
#effect of Ratio A on the PCs of raw data
gb %>% viz_gpca(guided=FALSE) %>% add(geom_point(aes(size=reference)))
ggplot(mapping=aes(x=t(ra)%*%gb$pca$v[,1], y=t(ra)%*%gb$pca$v[,2], colour=batch %>% factor))+geom_point(aes(size=reference))+stat_ellipse()
#effect of Ratio A on the PCs of corrected data
ggplot(mapping=aes(x=t(B)%*%gra$pca$v[,1], y=t(B)%*%gra$pca$v[,2], colour=batch %>% factor))+geom_point(aes(size=reference))+stat_ellipse()
gra %>% viz_gpca(guided=FALSE) %>% add(geom_point(aes(size=reference)))
#gPCA
gra %>% viz_gpca %>% add(geom_point(aes(size=reference)))
gra %>% viz_gpca_contrib
#the problem is that in a high dimension problem, samples are very far from each other (moreover since here all features are independent)
#every sample is an outlier in one direction and taking it as reference leads to spread the points in this direction even if it brings them closer in other directions

raw2<-cbind(
  get(load('data/E-MTAB-513-atlasExperimentSummary.Rdata')) %>% use_series(rnaseq) %>% assays %>% use_series(counts),# %>% apply(2,as_mapper(~rbind('mtab3871',.))),
  get(load('data/E-MTAB-3716-atlasExperimentSummary.Rdata')) %>% use_series(rnaseq) %>% assays %>% use_series(counts)# %>% apply(2,as_mapper(~rbind('mtab3716',.)))
)

batch<-c(
  'mtab513' %>% rep(get(load('data/E-MTAB-513-atlasExperimentSummary.Rdata')) %>% use_series(rnaseq) %>% dim %>% extract(2)),
  'mtab3716' %>% rep(get(load('data/E-MTAB-3716-atlasExperimentSummary.Rdata')) %>% use_series(rnaseq) %>% dim %>% extract(2))
)
graw<-raw2 %>% gPCA(batch,scale=T)
reference<-nearest.neighbours.references(raw2,batch,12)#FALSE %>% rep(ncol(raw2)) %>% replace(seq_along(.) %>% sample(2),TRUE) 
batch[which(reference)]#check
graw %>% viz_gpca %>% add(geom_point(aes(size=reference)))
ra<-raw2 %>% ratioG(batch,reference)
gra<-ra %>% gPCA(batch,scale=TRUE)
#effect of Ratio A on the PCs of raw data
graw %>% viz_gpca(guided=FALSE) %>% add(geom_point(aes(size=reference)))
ggplot(mapping=aes(x=t(ra[-graw$na.omit,])%*%graw$pca$v[,1], y=t(ra[-graw$na.omit,])%*%graw$pca$v[,2], colour=batch %>% factor))+geom_point(aes(size=reference))+stat_ellipse()
#effect of ratio A on the PCs of corrected data
ggplot(mapping=aes(x=t(raw2[-gra$na.omit,])%*%gra$pca$v[,1], y=t(raw2[-gra$na.omit,])%*%gra$pca$v[,2], colour=batch %>% factor))+geom_point(aes(size=reference))+stat_ellipse()
gra %>% viz_gpca(guided=FALSE) %>% add(geom_point(aes(size=reference)))
#gPCA of corrected data
gra %>% viz_gpca %>% add(geom_point(aes(size=reference)))

gra %>% viz_gpca_contrib


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

library(mvtnorm)
reference.choice<-function(data,batch,nref,npc=length(batch),PCs=1:npc,maximize=c('likelihood','odds.ratio')){
  batch %<>% factor
  batches<-batch %>% levels
  g<-data %>% gPCA(batch)
  LS<-batches %>% map(as_mapper(~g$pca$u[batch==.,PCs])) %>% map(as_mapper(
    ~list(mean=colMeans(.),sigma=var(.))
  ))
  likelihood<-LS %>% sapply(
    function(ls) g$pca$u[,PCs] %>% apply(MARGIN=1,
      function(sample) sample %>% dmvnorm(mean=ls$mean,sigma=ls$sigma)
    )
  ) %>% `colnames<-`(batches)
  odds.ratio<-likelihood %>% apply(MARGIN=1,
    function(sample) batches %>% sapply(as_mapper(
      ~sample[.]/(sum(sample)-sample[.])
    ))
  ) %>% t %>% `colnames<-`(batches)
  batch.refs<-maximize %>% get %>% apply(MARGIN=2,
    function(x) TRUE %>% rep(nref) %>% extract(x %>% order(decreasing=TRUE)) %>% replace(is.na(.),FALSE)
  )
  return(batch.refs %>% apply(1,any) %>% structure(
    by_batch=batch.refs,
    odds.ratio=odds.ratio,
    likelihood=likelihood
  ))
}

ref<-B %>% reference.choice(batch,10,npc=2,maximize='likelihood')
rar<-B %>% ratioA(batch, ref)
grar<-rar %>% gPCA(batch)
grar %>% viz_gpca(guided=FALSE) %>% add(geom_point(aes(size=ref)))
grar %>% viz_gpca %>% add(geom_point(aes(size=ref)))
ggplot(mapping=aes(x=t(rar)%*%gb$pca$v[,1], y=t(rar)%*%gb$pca$v[,2], colour=batch %>% factor))+geom_point(aes(size=ref))+stat_ellipse()
gb %>% viz_gpca(guided=FALSE) %>% add(geom_point(aes(size=ref)))

#before
ggplot(mapping=aes(x=t(B)%*%grar$pca$v[,1], y=t(B)%*%grar$pca$v[,2], colour=batch %>% factor))+geom_point(aes(size=ref))+stat_ellipse()
#after
grar %>% viz_gpca(guided=FALSE) %>% add(geom_point(aes(size=ref)))
#explanation of the issue : 
#



rarr<-rar %>% ratioA(batch,rar %>% reference.choice(batch,10,npc=10,maximize='odds.ratio'))
grarr<-rarr %>% gPCA(batch)
grarr %>% viz_gpca %>% add(geom_point(aes(size=ref)))

'effect of the number of reference points in ratio A results 15-30.pdf' %>% pdf(width=50,height=20)
gb<-B %>% gPCA(batch)
gb %>% viz_gpca(guided=FALSE) %>% print
for(n in 15:30){
  ref<-B %>% reference.choice(batch,n,npc=2)
  rar<-B %>% ratioA(batch,ref)
  grar<-rar %>% gPCA(batch)
  grid.arrange(ncol=3,
    ggplot(mapping=aes(x=t(rar)%*%gb$pca$v[,1], y=t(rar)%*%gb$pca$v[,2],colour=batch %>% factor))+geom_point(aes(size=ref))+stat_ellipse(),
    grar %>% viz_gpca %>% add(geom_point(aes(size=ref))),
    grar %>% viz_gpca(guided=FALSE) %>% add(geom_point(aes(size=ref)))
  ) %>% print
}
dev.off()


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