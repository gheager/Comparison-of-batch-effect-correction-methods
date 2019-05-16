#tests
library(magrittr)
library(purrr)
library(ExpressionAtlas)
library(plyr)
library(ggplot2)
library(patchwork)

compare.pca<-function(corrected,raw,batch,guided=FALSE){
  if(guided){
    raw.pca<-raw$gpca
    corrected.pca<-corrected$gpca
  }else{
    raw.pca<-raw$pca
    corrected.pca<-corrected$pca
  }
  #grid.arrange(nrow=2,
  wrap_plots(
    raw %>% viz_gpca(guided=guided) + theme(legend.position = 'none'),
    ggplot(mapping=aes(
      x=t(corrected$data)%*%raw.pca$v[,1],
      y=t(corrected$data)%*%raw.pca$v[,2],
      colour=batch
    ))+geom_point()+stat_ellipse() + theme(legend.position = 'none'),
    ggplot(mapping=aes(
      x=t(raw$data)%*%corrected.pca$v[,1],
      y=t(raw$data)%*%corrected.pca$v[,2],
      colour=batch
    ))+geom_point()+stat_ellipse() + theme(legend.position = 'none'),
    corrected %>% viz_gpca(guided=guided) + theme(legend.position = 'none')
  )
  #)
}

atlas<-dir('data') %>% lapply(as_mapper(
  ~get(load(paste0('data/',.)))
))

raw<-atlas %>% map(as_mapper(
  ~.$rnaseq %>% assays %>% use_series(counts)
)) %>% Reduce(cbind,.)

batch <- 1:5 %>% sapply(as_mapper(
  ~rep(.,dim(atlas[[.]]$rnaseq)[2])
)) %>% Reduce(c,.) %>% factor %>% `levels<-`(dir('data'))

tissue<-atlas %>% map(as_mapper(
  ~.$rnaseq$organism_part
)) %>% Reduce(c,.)

ggplot()+geom_bar(aes(x=tissue,colour=batch,fill=batch),width=1,position='dodge')+
  theme(axis.text.x=element_text(angle=45),legend.position = 'none')

#only heart tissue from now
# raw%<>%extract(,tissue=='heart')
# batch%<>%extract(tissue=='heart')

raw2 <- raw %>% extract(,batch==levels(batch)[4]|batch==levels(batch)[5])
batch2 <- batch %>% extract(batch==levels(batch)[4]|batch==levels(batch)[5])

raw2 %<>% log1p
raw2 %>% gPCA(batch2) -> graw
#raw2 %>% gPCA(batch2,nperm=100) -> grawp

library(pamr)
pam<-pamr.batchadjust(list(
  x=raw2,
  batchlabels=batch2
))$x
pam %>% gPCA(batch2) -> gpam
#gpam %>% viz_gpca(guided=FALSE)
gpam %>% compare.pca(graw,batch2)


library(sva)
t<-list(
  mtab4344=get(load('data/E-MTAB-4344-atlasExperimentSummary.Rdata')),
  mtab513=get(load('data/E-MTAB-513-atlasExperimentSummary.Rdata'))
)
all<-NULL
for (i in names(t)) {
  expAcc <- i
  k <- t[[i]]
  exp <- k$rnaseq
  eCounts <- assays(exp)$counts
  samples<-colnames(eCounts)
  average.counts<-technical_replicate_average_gtex(exp,expAcc)
  all <- cbind(all,average.counts)
}
filter <- rowSums(all>10)>=15
filtered <- all[filter,]
filterCols <- colSums(filtered == 0) / nrow(filtered) < 0.90
x <- x[filterCols]
filtered <- filtered[,filterCols]
x<-as.factor(sapply(strsplit(colnames(filtered),split="_"),"[",1))
combat<-ComBat(filtered %>% log1p,x)
combat %>% gPCA(x) -> gcombat
filtered %>% log1p %>% gPCA(x) -> gfiltered
gcombat %>% compare.pca(gfiltered,x)


ref<-raw2 %>% nearest.neighbours.references(batch2,2)
ra<-raw2 %>% ratioA(batch2,ref)
ra %>% gPCA(batch2) -> gra
gra %>% viz_gpca_contrib
gra %>% compare.pca(graw,batch2)


rg<-raw2 %>% ratioG(batch2,ref)
rg %>% gPCA(batch2) -> grg
grg %>% compare.pca(graw,batch2)
