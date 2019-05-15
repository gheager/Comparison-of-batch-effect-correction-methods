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
tissue<-as.factor(sapply(strsplit(colnames(filtered),split="_"),"[",3))


#### LOG-TRANSFORMATION
filtered%<>%log1p

filtered %>% gPCA(x) -> gfiltered

library(pamr)
pamr.batchadjust(list(x=filtered,batchlabels=x))$x->pam
pam %>% gPCA(x) -> gpam
gpam %>% compare.pca(gfiltered,x)
gpam %>% compare.pca(gfiltered,x,guided=TRUE)

library(sva)
ComBat(filtered,x)->combat
combat %>% gPCA(x) -> gcombat
gcombat %>% compare.pca(gfiltered,x)
gcombat %>% compare.pca(gfiltered,x,guided=TRUE)

ComBat(filtered,x,mod=model.matrix(~tissue %>% factor)) -> ccombat
ccombat %>% gPCA(x) -> gccombat
gccombat %>% compare.pca(gfiltered,x)
gccombat %>% compare.pca(gfiltered,x,guided=TRUE)


ref<-nearest.neighbours.references(filtered,x,2)
ratioA(filtered,x,ref)->ra
ra %>% gPCA(x) -> gra
gra %>% compare.pca(gfiltered,x)
gra %>% compare.pca(gfiltered,x,guided=TRUE)


ratioG(filtered,x,ref)->rg
rg %>% gPCA(x) -> grg
grg %>% compare.pca(gfiltered,x)
grg %>% compare.pca(gfiltered,x,guided=TRUE)



geneshift.pca(filtered,x,3,.001)->gs
gs %>% gPCA(x) -> ggs
ggs %>% compare.pca(gfiltered,x)
ggs %>% compare.pca(gfiltered,x,guided=TRUE)

filtered %>% batch_inertia(tissue)
pam %>% batch_inertia(tissue)
combat %>% batch_inertia(tissue)
ccombat %>% batch_inertia(tissue)
ra %>% batch_inertia(tissue)
rg %>% batch_inertia(tissue)
gs %>% batch_inertia(tissue)
