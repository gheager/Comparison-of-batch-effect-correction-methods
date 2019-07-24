library(magrittr)
library(SummarizedExperiment)
library(ggplot2)

geod74747<-get(load(
  'Mouse Expression Atlas/batch data/E-GEOD-74747-atlasExperimentSummary.Rdata'
))$rnaseq
mtab3725<-get(load(
  'Mouse Expression Atlas/batch data/E-MTAB-3725-atlasExperimentSummary.Rdata'
))$rnaseq

geod74747 %>% assays %$% counts -> geod
mtab3725 %>% assays %$% counts -> mtab
cbind(geod,mtab)->whole; batch<-c('geod' %>% rep(ncol(geod)),'mtab' %>% rep(ncol(mtab))); tissue<-c(geod74747$organism_part,mtab3725$organism_part)

whole %>% log1p %>% t %>% prcomp(rank.=2) -> pca
ggplot()+aes(pca$x[,1],y=pca$x[,2],colour=batch,group=tissue,label=tissue)+geom_text()+geom_line(colour='grey')
whole%<>%log1p
common.tissues<-intersect(tissue[batch=='geod'],tissue[batch=='mtab'])
whole%<>%extract(,tissue%in%common.tissues)
batch%<>%extract(tissue%in%common.tissues)
tissue%<>%extract(tissue%in%common.tissues)

filter<-rowSums(whole[,batch=='geod']!=0)>0 & rowSums(whole[,batch=='mtab']!=0)>0
filtered<-whole[filter,] %>% log1p
#filtered %>% log1p %>% t %>% prcomp(rank.=2) -> pca
#ggplot()+aes(x=pca$x[,1],y=pca$x[,2],colour=batch)+geom_point()+stat_ellipse()+geom_line(aes(group=tissue),colour='grey')

#correlation between coefficients of eigengenes and genes variances with and without scaling before PCA
v<-rowVars(filtered)
p<-filtered %>% t %>% prcomp(rank.=1) %$% rotation
#plot(v,p %>% abs)
ggplot()+aes(x=v,y=p %>% abs)+geom_point()+geom_smooth()

v<-rowVars(filtered)
p<-filtered %>% t %>% prcomp(rank.=1,scale=TRUE) %$% rotation
#plot(v,p %>% abs)
ggplot()+aes(x=v,y=p %>% abs)+geom_point()+geom_smooth()


eigenangles(filtered,batch,tissue)->ang
viz_angles_batch_vs_all(ang)
viz_angles_inter_batch(ang)


library(sva)
filtered %>% ComBat(batch,mod=model.matrix(~tissue)) -> combat
eigenangles(combat,batch,tissue)->angc
viz_angles_batch_vs_all(none=ang,combat=angc,arrow=anga)
viz_angles_inter_batch(none=ang,combat=angc,arrow=anga)

library(RUVSeq)


#ARROW
common.tissues<-intersect(
  tissue[batch=='geod'],
  tissue[batch=='mtab']
)
arrows<-NULL; for(o in common.tissues){
  arrows %<>% cbind(
    filtered[,batch=='geod' & tissue==o]-filtered[,batch=='mtab' & tissue==o]
  )
};colnames(arrows)<-common.tissues
mean.arrow<-rowMeans(arrows)
arw<-filtered
arw[,batch=='mtab'] %<>% add(mean.arrow)
arw %>% eigenangles(batch,tissue)->anga


#viz angles for article
whole[,batch=='geod'] %>% t %>% prcomp(rank.=1) -> pcageod
whole[,batch=='mtab'] %>% t %>% prcomp(rank.=1) -> pcamtab

geod74747 %>% gPCA.integrate('organism_part') -> gpcageod

eigenangles:::angledet(pcageod$rotation,pcamtab$rotation)->alpha
ggplot()+aes(x=pcageod$x,y=0,colour=tissue[batch=='geod'])+geom_point()+
  geom_point(aes(x=pcamtab$x*cospi(alpha),y=pcamtab$x*sinpi(alpha),colour=tissue[batch=='mtab']))+coord_fixed()

((filtered[,batch=='geod']+filtered[,batch=='mtab'])/2) %>% t %>% prcomp(rank.=1) -> pcaall 
filtered %>% t %>% prcomp(rank.=1) -> pcaall
eigenangles:::angledet(pcageod$rotation,pcaall$rotation)->alpha
ggplot()+aes(x=pcageod$x,y=0,colour=tissue[batch=='geod'])+geom_point()+
  geom_point(aes(x=pcaall$x*cospi(alpha),y=pcaall$x*sinpi(alpha),colour=tissue))+coord_fixed()
