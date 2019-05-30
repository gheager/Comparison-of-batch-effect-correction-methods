#TESTS ON MOUSE DATASETS
library(magrittr)
library(ExpressionAtlas)

geod74747<-get(load(
  'Mus musculus/~organism part/E-GEOD-74747-atlasExperimentSummary.Rdata'
))$rnaseq
mtab4644<-get(load(
  'Mus musculus/~organism part/E-MTAB-4644-atlasExperimentSummary.Rdata'
))$rnaseq

whole<-cbind(
  geod74747 %>% assays %>% use_series(counts),
  mtab4644 %>% assays %>% use_series(counts)
)
batch<-c(
  'geod' %>% rep(dim(geod74747)[2]),
  'mtab' %>% rep(dim(mtab4644)[2])
)
tissue<-c(
  geod74747$organism_part,
  mtab4644$organism_part
)
whole %>% log1p %>% t %>% prcomp -> pca
ggplot(mapping=aes(x=pca$x[,1],y=pca$x[,2],colour=batch))+
  geom_point()+stat_ellipse()+
  geom_line(aes(group=tissue))

filter<-rowSums(whole[,batch=='geod']!=0)>0 & rowSums(whole[,batch=='mtab']!=0)>0
filtered<-whole[filter,]
filtered %>% log1p %>% t %>% prcomp -> pca
ggplot(mapping=aes(x=pca$x[,1],y=pca$x[,2],colour=batch))+
  geom_point()+stat_ellipse()+
  geom_line(aes(group=tissue),colour='grey')

#LOG-TRANSFORMATION
filtered %<>% log1p
filtered %>% gPCA(batch,nperm=1000) -> gfiltered
gfiltered %>% viz_gpca
gfiltered %>% viz_gpca(guided=F) + geom_line(aes(group=tissue),colour='grey')+geom_text(aes(label=tissue))
gfiltered %>% viz_gpca_contrib
gfiltered %>% viz_gpca_pvalue

filtered %>% gPCA(batch)->gfa
filtered %>% gPCA(tissue)->gfb
gfa$delta/gfb$delta
gfb %>% viz_gpca()

gfa$gpca$v[,1] %>% angle(gfa$pca$v[,1])
gfb$gpca$v[,1] %>% angle(gfb$pca$v[,1])


library(pamr)
pamr.batchadjust(list(x=filtered,batchlabels=batch))$x->pam
pam %>% gPCA(batch)->gpa
pam %>% gPCA(tissue)->gpb
gpa$delta/gpb$delta

# gpa$gpca$v[,1] %>% angle(gpb$gpca$v[,1])
# gfa$gpca$v[,1] %>% angle(gfb$gpca$v[,1])

gpa$gpca$v[,1] %>% angle(gpa$pca$v[,1])
gpb$gpca$v[,1] %>% angle(gpb$pca$v[,1])

library(sva)
ComBat(filtered,batch)->combat
combat %>% gPCA(batch)->gca
combat %>% gPCA(tissue)->gcb
gca$gpca$v[,1] %>% angle(gca$pca$v[,1]) #%>% abs %>% acos/pi
gcb$gpca$v[,1] %>% angle(gcb$pca$v[,1]) #%>% abs %>% acos/pi

ComBat(filtered,batch,mod=model.matrix(~tissue))->ccombat
ccombat %>% gPCA(batch)->gcca
ccombat %>% gPCA(tissue)->gccb
gcca$gpca$v[,1] %>% angle(gcca$pca$v[,1])/gccb$gpca$v[,1] %>% angle(gccb$pca$v[,1])


#Arrow
angles(filtered,tissue,batch)
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

arw %>% gPCA(batch) -> garwa
arw %>% gPCA(tissue)->garwb
garwa$gpca$v[,1] %>% angle(garwa$pca$v[,1])/garwb$gpca$v[,1] %>% angle(garwb$pca$v[,1])

library(RUVSeq)
RUVs(filtered,
     cIdx=1:nrow(filtered),
     k=1,
     scIdx=makeGroups(tissue),
     # scIdx=c(1,10,
     #         2,11,
     #         3,12,
     #         4,13,
     #         5,14,
     #         7,16,
     #         8,17) %>% matrix(ncol=2,byrow = TRUE),
     isLog=TRUE)->ruvs
ruvs$normalizedCounts %>% gPCA(batch) -> grsa
ruvs$normalizedCounts %>% gPCA(tissue) -> grsb
grsa$gpca$v[,1] %>% angle(grsa$pca$v[,1])/grsb$gpca$v[,1] %>% angle(grsb$pca$v[,1])
grsa %>% viz_gpca(guided=F)+geom_line(aes(group=tissue),colour='grey')
grsb %>% viz_gpca
grsb %>% viz_gpca(guided=F)
grsa %>% compare.pca(gfiltered,batch,tissue)

library(harmony)
HarmonyMatrix(filtered,
              meta_data=data.frame(batch,tissue),
              vars_use = 'batch',
              do_pca=FALSE,nclust=2)->hm
#HarmonyMatrix(filtered,batch,npcs=5,nclust=5)->hm
hm %>% gPCA(batch)->ghma
hm %>% gPCA(tissue)->ghmb
ghma$gpca$v[,1] %>% angle(ghma$pca$v[,1])/ghmb$gpca$v[,1] %>% angle(ghmb$pca$v[,1])

ghm %>% compare.pca(gfiltered,batch,tissue,guided=T)


###
common.tissues<-intersect(
  tissue[batch=='geod'],
  tissue[batch=='mtab']
)
whole %>% extract(,tissue %in% common.tissues) -> whole1
batch %>% extract(tissue %in% common.tissues) -> batch1
tissue %>% extract(tissue %in% common.tissues) -> tissue1

filter<-
  rowSums(whole1[,batch1=='geod']>0)>0 & 
  rowSums(whole1[,batch1=='mtab']>0)>0
whole1%<>%extract(filter,)

#LOG TRANSFORMATION
whole1%<>%log1p

whole1[,batch1=='geod'] %>% t %>% prcomp -> pcb1
whole1[,batch1=='mtab'] %>% t %>% prcomp -> pcb2
whole1 %>% t %>% prcomp -> pcb3
angle(pcb1$rotation[,1],pcb2$rotation[,1])
angle(pcb1$rotation[,1],pcb3$rotation[,1])
angle(pcb2$rotation[,1],pcb3$rotation[,1])

ggplot()+
  geom_point(aes(x=pcb1$rotation[,1],y=pcb1$rotation[,2]),colour='red',pch=4)+
  geom_point(aes(x=pcb2$rotation[,1],y=pcb2$rotation[,2]),colour='blue',pch=4)+
  geom_point(aes(x=pcb3$rotation[,1],y=pcb3$rotation[,2]),colour='green',pch=4)

ggplot()+
  geom_density(aes(x=pcb1$rotation[,2]),colour='red')+
  geom_density(aes(x=pcb2$rotation[,2]),colour='blue')+
  geom_density(aes(x=pcb3$rotation[,2]),colour='green')

#number of genes
ngenes<-nrow(whole1)

ggplot()+
  geom_point(aes(x=pcb1$rotation[,1] %>% abs,y=pcb1$rotation[,2] %>% abs),colour='red')+
  geom_point(aes(x=pcb2$rotation[,1] %>% abs,y=pcb2$rotation[,2] %>% abs),colour='blue')+
  geom_point(aes(x=pcb3$rotation[,1] %>% abs,y=pcb3$rotation[,2] %>% abs),colour='green')+
  geom_line(aes(
    x=c(pcb1$rotation[,1] %>% abs,pcb2$rotation[,1] %>% abs),
    y=c(pcb1$rotation[,2] %>% abs,pcb2$rotation[,2] %>% abs),
    group=1:ngenes %>% rep(2)
  ),colour='grey',alpha=.7)

ggplot()+
  geom_point(aes(x=pcb1$rotation[,1],y=pcb1$rotation[,2]),colour='red')+
  geom_point(aes(x=pcb2$rotation[,1],y=pcb2$rotation[,2]),colour='blue')+
  geom_line(aes(
    x=c(pcb1$rotation[,1],pcb2$rotation[,1]),
    y=c(pcb1$rotation[,2],pcb2$rotation[,2]),
    group=1:ngenes %>% rep(2)
  ),colour='grey',alpha=.7)

ggplot()+
  geom_point(aes(x=pcb1$rotation[,1],y=pcb1$rotation[,2]),colour='red')+
  geom_point(aes(x=pcb2$rotation[,1]+.2,y=pcb2$rotation[,2]+.2),colour='blue')+
  geom_point(aes(x=pcb3$rotation[,1],y=pcb3$rotation[,2]+.2),colour='green')+
  geom_line(aes(
    x=c(pcb1$rotation[,1],pcb2$rotation[,1]+.2),
    y=c(pcb1$rotation[,2],pcb2$rotation[,2]+.2),
    group=1:ngenes %>% rep(2)
  ),colour='magenta',alpha=.2)+
  geom_line(aes(
    x=c(pcb1$rotation[,1],pcb3$rotation[,1]),
    y=c(pcb1$rotation[,2],pcb3$rotation[,2]+.2),
    group=1:ngenes %>% rep(2)
  ),colour='yellow',alpha=.2)+
  geom_line(aes(
    x=c(pcb2$rotation[,1]+.2,pcb3$rotation[,1]),
    y=c(pcb2$rotation[,2]+.2,pcb3$rotation[,2]+.2),
    group=1:ngenes %>% rep(2)
  ),colour='cyan',alpha=.2)

threshold<-.08
ggplot()+
  geom_point(aes(x=pcb1$rotation[,1],y=pcb1$rotation[,2]),colour='red')+
  geom_point(aes(x=pcb2$rotation[,1]+.2,y=pcb2$rotation[,2]+.2),colour='blue')+
  geom_point(aes(x=pcb3$rotation[,1],y=pcb3$rotation[,2]+.2),colour='green')+
  geom_line(aes(
    x=c(
      ifelse(pcb1$rotation[,1] %>% abs>threshold|pcb1$rotation[,2] %>% abs>threshold|
               pcb2$rotation[,1] %>% abs>threshold|pcb2$rotation[,2] %>% abs>threshold,
             pcb1$rotation[,1],0),
      ifelse(pcb1$rotation[,1] %>% abs>threshold|pcb1$rotation[,2] %>% abs>threshold|
               pcb2$rotation[,1] %>% abs>threshold|pcb2$rotation[,2] %>% abs>threshold,
             pcb2$rotation[,1]+.2,0)
    ),
    y=c(
      ifelse(pcb1$rotation[,1] %>% abs>threshold|pcb1$rotation[,2] %>% abs>threshold|
               pcb2$rotation[,1] %>% abs>threshold|pcb2$rotation[,2] %>% abs>threshold,
             pcb1$rotation[,2],0),
      ifelse(pcb1$rotation[,1] %>% abs>threshold|pcb1$rotation[,2] %>% abs>threshold|
               pcb2$rotation[,1] %>% abs>threshold|pcb2$rotation[,2] %>% abs>threshold,
             pcb2$rotation[,2]+.2,0)
    ),
    group=1:ngenes %>% rep(2) %>% factor
  ),colour='green',alpha=.2)#+
  geom_line(aes(
    x=c(pcb1$rotation[,1],pcb3$rotation[,1]),
    y=c(pcb1$rotation[,2],pcb3$rotation[,2]+.2),
    group=1:ngenes %>% rep(2)
  ),colour='yellow',alpha=.2)+
  geom_line(aes(
    x=c(pcb2$rotation[,1]+.2,pcb3$rotation[,1]),
    y=c(pcb2$rotation[,2]+.2,pcb3$rotation[,2]+.2),
    group=1:ngenes %>% rep(2)
  ),colour='cyan',alpha=.2)

ggplot()+aes(
  x=pcb1$rotation[,1],
  y=pcb2$rotation[,1]
)+geom_point()+geom_smooth(method='lm')

pcb1$rotation[,1] %>% cor(pcb2$rotation[,1])

ggplot()+aes(
  x=c(pcb1$rotation[,1],pcb2$rotation[,1]),
  y=c(pcb1$rotation[,2],pcb2$rotation[,2]),
  group=1:ngenes %>% rep(each=2),
  colour=c('red','blue') %>% rep(each=ngenes)
)+geom_point()+geom_line(colour='grey',alpha=.5)

#correspondance map on first PC for genes
# ggplot()+aes(
#   x=c(0,1) %>% rep(each=ngenes),
#   y=c(pcb1$rotation[,1],pcb2$rotation[,1]),
#   colour=c(0,1) %>% rep(each=ngenes)
# )+geom_point()+
#   geom_line(aes(group=1:ngenes %>% rep(2)),colour='grey')+
#   theme(legend.position = 'none')

ggplot()+aes(
  x=c(0,1) %>% rep(each=ngenes),
  y=c(pcb1$rotation[,1],pcb2$rotation[,1]) %>% abs,
  colour=c('geod','mtab') %>% rep(each=ngenes)
)+geom_point()+
  geom_line(aes(group=1:ngenes %>% rep(2)),colour='grey',alpha=.7)+
  theme(legend.position = 'none')

#PCA of first batch
ggplot()+aes(
  x=pcb1$x[,1],
  y=pcb1$x[,2]
)+geom_point(colour='red')+geom_text(aes(label=tissue1[batch1=='geod']))
#PCA of second batch
ggplot()+aes(
  x=pcb2$x[,1],
  y=pcb2$x[,2]
)+geom_point(colour='blue')+geom_text(aes(label=tissue1[batch1=='mtab']))

#PCA of first batch with projection of second batch
pcb1 %>% predict(whole1[,batch1=='mtab'] %>% t) -> proj1
ggplot()+
  aes(
    x=c(pcb1$x[,1],proj1[,'PC1']),
    y=c(pcb1$x[,2],proj1[,'PC2']),
    colour=batch1,
    group=tissue1,
    label=tissue1
  )+geom_point()+geom_line(colour='grey')+geom_text()

#PCA of second batch with projection of first batch
pcb2 %>% predict(whole1[,batch1=='geod'] %>% t) -> proj2
ggplot()+
  aes(
    x=c(pcb2$x[,1],proj2[,'PC1']),
    y=c(pcb2$x[,2],proj2[,'PC2']),
    colour=batch1,
    group=tissue1,
    label=tissue1
  )+geom_point()+geom_line(colour='grey')+geom_text()

ggplot()+aes(
  x=pcb1$rotation[,2] %>% abs,
  y=pcb2$rotation[,1] %>% abs
)+geom_point()+geom_smooth(method='lm')

1:ngenes %>% sapply(as_mapper(
  ~norm(as.matrix(pcb1$rotation[.,1:2]-pcb2$rotation[.,1:2]),'F')
)) %>% summary
1:ngenes %>% sapply(as_mapper(
  ~norm(as.matrix(pcb1$rotation[.,1:2]-pcb3$rotation[.,1:2]),'F')
)) %>% summary
1:ngenes %>% sapply(as_mapper(
  ~norm(as.matrix(pcb3$rotation[.,1:2]-pcb2$rotation[.,1:2]),'F')
)) %>% summary

library(pamr)
pamr.batchadjust(list(x=whole1,batchlabels=batch1))$x->pam
pam[,batch1=='geod'] %>% t %>% prcomp -> pcb1p
pam[,batch1=='mtab'] %>% t %>% prcomp -> pcb2p
pam %>% t %>% prcomp -> pcb3p
angle(pcb1p$rotation[,1],pcb2p$rotation[,1])
angle(pcb1p$rotation[,1],pcb3p$rotation[,1])
angle(pcb3p$rotation[,1],pcb2p$rotation[,1])

ggplot()+
  geom_point(aes(x=pcb1p$rotation[,1],y=pcb1p$rotation[,2]),colour='red')+
  geom_point(aes(x=pcb2p$rotation[,1],y=pcb2p$rotation[,2]),colour='blue')+
  geom_point(aes(x=pcb3p$rotation[,1],y=pcb3p$rotation[,2]),colour='green')

ggplot()+
  geom_density(aes(x=pcb1p$rotation[,2]),colour='red')+
  geom_density(aes(x=pcb2p$rotation[,2]),colour='blue')+
  geom_density(aes(x=pcb3p$rotation[,2]),colour='green')

ggplot()+
  geom_point(aes(x=pcb1p$rotation[,1],y=pcb1p$rotation[,2]),colour='red')+
  geom_point(aes(x=pcb2p$rotation[,1]+.2,y=pcb2p$rotation[,2]+.2),colour='blue')+
  geom_point(aes(x=pcb3p$rotation[,1],y=-pcb3p$rotation[,2]+.2),colour='green')+
  geom_line(aes(
    x=c(pcb1p$rotation[,1],pcb2p$rotation[,1]+.2),
    y=c(pcb1p$rotation[,2],pcb2p$rotation[,2]+.2),
    group=1:ngenes %>% rep(2)
  ),colour='magenta',alpha=.2)+
  geom_line(aes(
    x=c(pcb1p$rotation[,1],pcb3p$rotation[,1]),
    y=c(pcb1p$rotation[,2],-pcb3p$rotation[,2]+.2),
    group=1:ngenes %>% rep(2)
  ),colour='yellow',alpha=.2)+
  geom_line(aes(
    x=c(pcb2p$rotation[,1]+.2,pcb3p$rotation[,1]),
    y=c(pcb2p$rotation[,2]+.2,-pcb3p$rotation[,2]+.2),
    group=1:ngenes %>% rep(2)
  ),colour='cyan',alpha=.2)

ggplot()+
  geom_point(aes(x=pcb1p$rotation[,1] %>% abs,y=pcb1p$rotation[,2] %>% abs),colour='red')+
  geom_point(aes(x=pcb2p$rotation[,1] %>% abs,y=pcb2p$rotation[,2] %>% abs),colour='blue')+
  geom_line(aes(
    x=c(pcb1p$rotation[,1] %>% abs,pcb2p$rotation[,1] %>% abs),
    y=c(pcb1p$rotation[,2] %>% abs,pcb2p$rotation[,2] %>% abs),
    group=1:ngenes %>% rep(2)
  ),colour='grey',alpha=.7)

ggplot()+
  geom_point(aes(x=pcb1$rotation[,1],y=pcb1$rotation[,2]),colour='red')+
  geom_point(aes(x=pcb2$rotation[,1],y=pcb2$rotation[,2]),colour='blue')+
  geom_line(aes(
    x=c(pcb1$rotation[,1],pcb2$rotation[,1]),
    y=c(pcb1$rotation[,2],pcb2$rotation[,2]),
    group=1:ngenes %>% rep(2)
  ),colour='grey',alpha=.7)

1:ngenes %>% sapply(as_mapper(
  ~norm(as.matrix(pcb1p$rotation[.,1:2]-pcb2p$rotation[.,1:2]),'F')
)) %>% summary
1:ngenes %>% sapply(as_mapper(
  ~norm(as.matrix(pcb1p$rotation[.,1:2]-pcb3p$rotation[.,1:2]),'F')
)) %>% summary
1:ngenes %>% sapply(as_mapper(
  ~norm(as.matrix(pcb3p$rotation[.,1:2]-pcb2p$rotation[.,1:2]),'F')
)) %>% summary

library(sva)
ComBat(whole1,batch1)->combat
combat[,batch1=='geod'] %>% t %>% prcomp -> pcb1combat
combat[,batch1=='mtab'] %>% t %>% prcomp -> pcb2combat
combat %>% t %>% prcomp -> pcb3combat
angle(pcb1combat$rotation[,1],pcb2combat$rotation[,1])
angle(pcb1combat$rotation[,1],pcb3combat$rotation[,1])
angle(pcb3combat$rotation[,1],pcb2combat$rotation[,1])

ggplot()+
  geom_point(aes(x=pcb1combat$rotation[,1],y=pcb1combat$rotation[,2]),colour='red')+
  geom_point(aes(x=pcb2combat$rotation[,1],y=pcb2combat$rotation[,2]),colour='blue')+
  geom_point(aes(x=pcb3combat$rotation[,1],y=pcb3combat$rotation[,2]),colour='green')

ggplot()+
  geom_density(aes(x=pcb1combat$rotation[,2]),colour='red')+
  geom_density(aes(x=pcb2combat$rotation[,2]),colour='blue')+
  geom_density(aes(x=pcb3combat$rotation[,2]),colour='green')

ggplot()+
  geom_point(aes(x=pcb1combat$rotation[,1],y=pcb1combat$rotation[,2]),colour='red')+
  geom_point(aes(x=pcb2combat$rotation[,1]+.2,y=pcb2combat$rotation[,2]+.2),colour='blue')+
  geom_point(aes(x=pcb3combat$rotation[,1],y=-pcb3combat$rotation[,2]+.2),colour='green')+
  geom_line(aes(
    x=c(pcb1combat$rotation[,1],pcb2combat$rotation[,1]+.2),
    y=c(pcb1combat$rotation[,2],pcb2combat$rotation[,2]+.2),
    group=1:ngenes %>% rep(2)
  ),colour='magenta',alpha=.2)+
  geom_line(aes(
    x=c(pcb1combat$rotation[,1],pcb3combat$rotation[,1]),
    y=c(pcb1combat$rotation[,2],-pcb3combat$rotation[,2]+.2),
    group=1:ngenes %>% rep(2)
  ),colour='yellow',alpha=.2)+
  geom_line(aes(
    x=c(pcb2combat$rotation[,1]+.2,pcb3combat$rotation[,1]),
    y=c(pcb2combat$rotation[,2]+.2,-pcb3combat$rotation[,2]+.2),
    group=1:ngenes %>% rep(2)
  ),colour='cyan',alpha=.2)

1:ngenes %>% sapply(as_mapper(
  ~norm(as.matrix(pcb1combat$rotation[.,1:2]-pcb2combat$rotation[.,1:2]),'F')
)) %>% summary
1:ngenes %>% sapply(as_mapper(
  ~norm(as.matrix(pcb1combat$rotation[.,1:2]-pcb3combat$rotation[.,1:2]),'F')
)) %>% summary
1:ngenes %>% sapply(as_mapper(
  ~norm(as.matrix(pcb3combat$rotation[.,1:2]-pcb2combat$rotation[.,1:2]),'F')
)) %>% summary


###
ComBat(whole1,batch1,model.matrix(~tissue1))->ccombat
ccombat[,batch1=='geod'] %>% t %>% prcomp -> pcb1ccombat
ccombat[,batch1=='mtab'] %>% t %>% prcomp -> pcb2ccombat
ccombat %>% t %>% prcomp -> pcb3ccombat
angle(pcb1ccombat$rotation[,1],pcb2ccombat$rotation[,1])
angle(pcb1ccombat$rotation[,1],pcb3ccombat$rotation[,1])
angle(pcb3ccombat$rotation[,1],pcb2ccombat$rotation[,1])
