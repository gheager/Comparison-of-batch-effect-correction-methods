#tests on log-transformed filtered real data
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
  grid.arrange(ncol=2,nrow=2,
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

t<-list(
  mtab4344=get(load('data/E-MTAB-4344-atlasExperimentSummary.Rdata')),
  mtab513=get(load('data/E-MTAB-513-atlasExperimentSummary.Rdata')),
  mtab3716=get(load('data/E-MTAB-3716-atlasExperimentSummary.Rdata'))
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
filtered <- filtered[,filterCols]
x<-as.factor(sapply(strsplit(colnames(filtered),split="_"),"[",1))
tissue<-as.factor(sapply(strsplit(colnames(filtered),split="_"),"[",3))


#### LOG-TRANSFORMATION
filtered%<>%log1p

filtered %>% gPCA(x) -> gfiltered
gfiltered %>% viz_gpca_contrib

library(pamr)
pamr.batchadjust(list(x=filtered,batchlabels=x))$x->pam
pam %>% gPCA(x) -> gpam
gpam %>% compare.pca(gfiltered,x)
gpam %>% compare.pca(gfiltered,x,guided=TRUE)
gpam %>% viz_gpca_contrib

library(sva)
ComBat(filtered,x)->combat
combat %>% gPCA(x) -> gcombat
gcombat %>% compare.pca(gfiltered,x)
gcombat %>% compare.pca(gfiltered,x,guided=TRUE)
gcombat %>% viz_gpca_contrib

ComBat(filtered,x,mod=model.matrix(~tissue %>% factor)) -> ccombat
ccombat %>% gPCA(x) -> gccombat
gccombat %>% compare.pca(gfiltered,x)
gccombat %>% compare.pca(gfiltered,x,guided=TRUE)
gccombat %>% viz_gpca_contrib


ref<-nearest.neighbours.references(filtered,x,2)
ratioA(filtered,x,ref)->ra
ra %>% gPCA(x) -> gra
gra %>% compare.pca(gfiltered,x)
gra %>% compare.pca(gfiltered,x,guided=TRUE)
gra %>% viz_gpca_contrib


ratioG(filtered,x,ref)->rg
rg %>% gPCA(x) -> grg
grg %>% compare.pca(gfiltered,x)
grg %>% compare.pca(gfiltered,x,guided=TRUE)
grg %>% viz_gpca_contrib


geneshift.pca(filtered,x,3,.01)->gs
gs %>% gPCA(x) -> ggs
ggs %>% compare.pca(gfiltered,x)
ggs %>% compare.pca(gfiltered,x,guided=TRUE)
ggs %>% viz_gpca_contrib


library(harmony)
HarmonyMatrix(filtered,data.frame(x,tissue),vars_use = 'x',do_pca=F,npcs=14,nclust=2)->hm
#hm %>% {.%*%diag(gfiltered$pca$d[1:14])} %>% t %>% batch_inertia(tissue)
hm%>% batch_inertia(tissue)
ggplot(mapping=aes(x=-hm[,1],y=-hm[,2],colour=x))+geom_point()+stat_ellipse()+geom_line(aes(group=tissue))
gfiltered %>% viz_gpca(guided=FALSE,dims=3:4)+geom_line(aes(group=tissue))#+geom_text(aes(label=tissue),size=3)

hmr<-t(hm%*%diag(gfiltered$pca$d[1:14])%*%t(gfiltered$pca$v[,1:14]))
hm %>% gPCA(x) -> ghm
ghm %>% viz_gpca_contrib
ghm %>% compare.pca(gfiltered,x)
ghm %>% compare.pca(gfiltered,x,guided=TRUE)
batch.inertia(hmr,tissue)$total

#intra-class inertia indices
filtered %>% batch_inertia(tissue)
pam %>% batch_inertia(tissue)
combat %>% batch_inertia(tissue)
ccombat %>% batch_inertia(tissue)
ra %>% batch_inertia(tissue)
rg %>% batch_inertia(tissue)
gs %>% batch_inertia(tissue)


#volumes indices
batch.inertia(filtered,tissue)$cumvars
batch.inertia(pam,tissue)$volume.index
batch.inertia(combat,tissue)$volume.index
batch.inertia(ccombat,tissue)$volume.index
batch.inertia(ra,tissue)$volume.index
batch.inertia(rg,tissue)$volume.index
batch.inertia(gs,tissue)$volume.index
batch.inertia(hm,tissue,transpose=FALSE)$volume.index

ggplot(mapping=aes(x=1:29))+
  geom_line(aes(y=batch.inertia(filtered,tissue)$volume.index %>% log),colour='blue')+
  geom_line(aes(y=batch.inertia(pam,tissue)$volume.index %>% log),colour='red')+
  geom_line(aes(y=batch.inertia(combat,tissue)$volume.index %>% log),colour='green')+
  geom_line(aes(y=batch.inertia(ccombat,tissue)$volume.index %>% log),colour='yellow')+
  geom_line(aes(y=batch.inertia(ra,tissue)$volume.index %>% log),colour='brown')+
  geom_line(aes(y=batch.inertia(rg,tissue)$volume.index %>% log),colour='black')+
  geom_line(aes(y=batch.inertia(gs,tissue)$volume.index %>% log),colour='orange')+
  geom_line(aes(y=batch.inertia(hm,tissue,transpose=FALSE)$volume.index %>% log),colour='purple')+
  xlim(c(0,10))+ylim(c(-25,0))

ggplot(mapping=aes(x=1:14))+
  geom_line(aes(y=batch.inertia(filtered,tissue)$volume.index[1:14] %>% log),colour='blue')+
  geom_line(aes(y=batch.inertia(hm,tissue,transpose=FALSE)$volume.index %>% log),colour='red')

batch_inertia(filtered,tissue)
batch_inertia(pam,tissue)
batch_inertia(combat,tissue)
batch_inertia(ccombat,tissue)
batch_inertia(ra,tissue)
batch_inertia(rg,tissue)
batch_inertia(gs,tissue)
batch_inertia(hm,tissue)

