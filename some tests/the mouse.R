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

filtered %>% gPCA(batch,nperm=1000) -> gfiltered
gfiltered %>% viz_gpca + geom_line(aes(group=tissue),colour='grey')
gfiltered %>% viz_gpca_contrib
gfiltered %>% viz_gpca_pvalue

#LOG-TRANSFORMATION
filtered %<>% log1p
filtered %>% gPCA(batch,nperm=1000) -> gfiltered
gfiltered %>% viz_gpca + geom_line(aes(group=tissue),colour='grey')
gfiltered %>% viz_gpca_contrib
gfiltered %>% viz_gpca_pvalue

library(pamr)
pamr.batchadjust(list(x=filtered,batchlabels=batch))$x -> pam
pam %>% gPCA(batch,scaleY=TRUE,nperm=1000) -> gpam
gpam %>% compare.pca(gfiltered,batch,tissue)

library(sva)
ComBat(filtered,batch) -> combat
combat %>% gPCA(batch,nperm=1000) -> gcombat
gcombat %>% compare.pca(gfiltered,batch,tissue)

ComBat(filtered,batch,mod=model.matrix(~tissue)) -> ccombat
ccombat %>% gPCA(batch,nperm=1000) -> gccombat
gccombat %>% compare.pca(gfiltered,batch,tissue)

library(harmony)
HarmonyMatrix(filtered,
              meta_data=data.frame(batch,tissue),
              vars_use = 'batch',
              do_pca=FALSE,nclust=2)->hm
#HarmonyMatrix(filtered,batch,npcs=5,nclust=5)->hm
hm %>% gPCA(batch)->ghm
ghm %>% compare.pca(gfiltered,batch,tissue)

library(RUVSeq)
rowVars(filtered[,batch=='geod']) %>% order %>% head
rowVars(filtered[,batch=='mtab']) %>% order %>% head
ggplot()+
  geom_bar(aes(x=1:nrow(filtered),
               y=1/rowVars(filtered[,batch=='geod'])),
           width=.5,
           stat='identity',
           fill='red')+
  geom_bar(aes(x=1:nrow(filtered),
               y=1/rowVars(filtered[,batch=='mtab'])),
           width=.5,
           stat='identity',
           fill='blue',position=.5)
o<-c(
  rowVars(filtered[,batch=='geod']),
  rowVars(filtered[,batch=='mtab'])
) %>% order %>% mod(nrow(filtered)) %>% add(1)
seen<-NULL;ranking<-data.frame(matrix(ncol=3,nrow=0)) %>% set_colnames(c('gene_number','vargeod','varmtab'))
for(i in o){
  if(i %in% seen){
    ranking%<>%rbind(
      c(i,var(filtered[i,batch=='geod']),var(filtered[i,batch=='mtab']))
    )
  }
  seen%<>%c(i)
}; ranking %<>% set_colnames(c('gene_number','vargeod','varmtab'))

RUVg(filtered,
     ranking$gene_number[1:100],
     k=1,isLog=TRUE)->ruvg
ruvg$normalizedCounts %>% gPCA(batch) -> gruvg
gruvg %>% compare.pca(gfiltered,batch,tissue)

ggplot()+
  geom_density(aes(x=rowVars(filtered[,batch=='geod'])),colour='blue')+
  geom_density(aes(x=rowVars(filtered[,batch=='mtab'])),colour='red')
rowVars(filtered[,batch=='geod']) %>% summary
rowVars(filtered[,batch=='mtab']) %>% summary

RUVg(filtered,
     which(
       rowVars(filtered[,batch=='geod'])<3 & rowVars(filtered[,batch=='mtab'])<3
     ),
     k=1,isLog=TRUE)->ruvg
ruvg$normalizedCounts %>% gPCA(batch) -> gruvg
gruvg %>% compare.pca(gfiltered,batch,tissue)


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

arw %>% gPCA(batch,nperm=1000) -> garw
garw %>% compare.pca(gfiltered,batch,tissue)

#scale on arrow-corrected data
filter<-rowSums(whole[,batch=='geod' & tissue%in%common.tissues]!=0)>0 & rowSums(whole[,batch=='mtab' & tissue%in%common.tissues]!=0)>0
filtered<-whole[filter,]
#LOG-TRANSFORMATION
filtered %<>% log1p
filtered %>% gPCA(batch,nperm=1000) -> gfiltered
gfiltered %>% viz_gpca + geom_line(aes(group=tissue),colour='grey')

#normal arrow method before scaling
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
#scaling
cbind(
  rowSds(arw[,batch=='geod' & tissue%in%common.tissues]),
  rowSds(arw[,batch=='mtab' & tissue%in%common.tissues])
)
arwsd<-arw
center<-rowMeans(arwsd[,batch=='mtab' & tissue%in%common.tissues])
arwsd[,batch=='mtab' & tissue%in%common.tissues] %<>%
  subtract(center) %<>% 
  divide_by(rowSds(.)) %<>% 
  multiply_by(rowSds(arwsd[,batch=='geod' & tissue%in%common.tissues])) %<>%
  add(center)
cbind(
  rowSds(arwsd[,batch=='geod' & tissue%in%common.tissues]),
  rowSds(arwsd[,batch=='mtab' & tissue%in%common.tissues])
)
arwsd %>% gPCA(batch,nperm=1000) -> garwsd
garwsd %>% compare.pca(gfiltered,batch,tissue)

#length of the arrows
rbind(
  tissue.distances(filtered,tissue,batch) %>% summary,
  tissue.distances(pam,tissue,batch) %>% summary,
  tissue.distances(combat,tissue,batch) %>% summary,
  tissue.distances(ccombat,tissue,batch) %>% summary,
  tissue.distances(arw,tissue,batch) %>% summary,
  tissue.distances(arwsd,tissue,batch) %>% summary
)
