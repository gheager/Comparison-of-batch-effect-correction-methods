library(magrittr)
library(stringr)
library(purrr)
library(ExpressionAtlas)

experiments<-list()
experiments[['GEOD30720']]<-get(load('Arabidopsis thaliana/~ecotype/E-GEOD-30720-atlasExperimentSummary.Rdata'))$rnaseq
experiments[['GEOD53197']]<-get(load('Arabidopsis thaliana/~ecotype + organism part/E-GEOD-53197-atlasExperimentSummary.Rdata'))$rnaseq


data<-NULL;batch<-NULL;tissues<-NULL; for(i in experiments %>% seq_along){
  data%<>%cbind(experiments[[i]] %>% assays %$% counts)
  batch%<>%c(names(experiments)[[i]] %>% rep(dim(experiments[[i]])[2]))
  tissues[[i]]<-(experiments[[i]]$ecotype)
};tissue<-Reduce(c,tissues)

batch[batch=='GEOD53197']%<>%paste(experiments$GEOD53197$organism_part)
batch%<>%factor
# intersections<-NULL; for(i in 1:5){
#   for(j in 1:5){
#     for(k in 1:5){
#       intersections%<>%c(tissues[[i]] %>% intersect(tissues[[j]]) %>% intersect(tissues[[k]]) %>% length)
#     }
#   }
# };intersections%<>%array(5 %>% rep(3))

# #common tissues only
# common.tissues<-Reduce(intersect,tissues)
# data%<>%extract(,tissue%in%common.tissues)
# batch%<>%extract(tissue%in%common.tissues)
# tissue%<>%extract(tissue%in%common.tissues)
# 
# #merging replicates
# data[,5] %<>% add(data[,6]) %<>% divide_by(2); data%<>%extract(,-6); batch%<>%extract(-6); tissue%<>%extract(-6)

filter<-Reduce(`&`,batch %>% levels %>% lapply(as_mapper(~rowSums(data[,batch==.])!=0)))
filtered<-data[filter,]

filtered %>% log1p %>% gPCA(batch,nperm=100) -> gfiltered
gfiltered %>% viz_gpca(guided=FALSE) + geom_line(aes(group=tissue),colour='grey')

filtered%<>%log1p
library(pamr)
filtered %>% list(x=.,batchlabels=batch) %>% pamr.batchadjust %$% x -> pam
pam %>% gPCA(batch) -> gpam
gpam %>% viz_gpca(guided=FALSE) + geom_line(aes(group=tissue),colour='grey')
gpam %>% compare.pca(gfiltered,batch,tissue)

library(sva)
filtered %>% ComBat(batch,mod=model.matrix(~tissue)) -> combat
combat %>% gPCA(batch,nperm=100) -> gcombat
gcombat %>% viz_gpca(guided = FALSE) + geom_line(aes(group=tissue),colour='grey')
gcombat %>% compare.pca(gfiltered,batch,tissue)

library(RUVSeq)
RUVs(filtered,
     rownames(filtered),
     k=20,
     makeGroups(tissue),
     isLog=TRUE)->ruvs
ruvs$normalizedCounts %>% gPCA(batch,nperm=100) -> gruvs
gruvs %>% compare.pca(gfiltered,batch,tissue)

batchgeneangles2(filtered,batch,tissue)
batchgeneangles2(pam,batch,tissue)
batchgeneangles2(combat,batch,tissue)
ruvs$normalizedCounts %>% batchgeneangles2(batch,tissue)

cc<-combat
m<-rowMeans(cc[,batch=='MTAB4644'])
cc[,batch=='MTAB4644']%<>%
  subtract(m)%<>%
  divide_by(rowSds(.))%<>%
  multiply_by(rowSds(cc[,batch=='GEOD74747']))%<>%
  add(m)
cc %>% gPCA(batch)->gcc

m<-rowMeans(cc[,batch=='MTAB3725'])
cc[,batch=='MTAB3725']%<>%
  subtract(m)%<>%
  divide_by(rowSds(.))%<>%
  multiply_by(rowSds(cc[,batch=='GEOD74747']))%<>%
  add(m)
cc %>% gPCA(batch,nperm=100)->gcc

gcc %>% viz_gpca(guided = FALSE) + geom_line(aes(group=tissue),colour='grey')

batchgeneangles(filtered,batch)
batchgeneangles(pam,batch)
batchgeneangles(combat,batch)
batchgeneangles(cc,batch)


#common tissues only
common.tissues<-Reduce(intersect,tissues)
combat %>% extract(,tissue%in%common.tissues) -> combat.common
batch %>% extract(tissue%in%common.tissues) -> batch.common
tissue %>% extract(tissue%in%common.tissues) -> tissue.common

#merging replicates
combat.common[,5] %<>% add(combat.common[,6]) %<>% divide_by(2); combat.common%<>%extract(,-6); batch.common%<>%extract(-6); tissue.common%<>%extract(-6)

batchgeneangles(filtered,batch)
batchgeneangles(pam,batch)
batchgeneangles(combat.common,batch.common)
batchgeneangles(cc,batch)
