library(magrittr)
library(stringr)
library(purrr)
library(ExpressionAtlas)

experiments<-list();for(filename in dir("Mus musculus/~organism part + developmental stage")[1:2]){
  experiments[[filename %>% str_split("-") %>% extract2(1) %>% extract(2:3) %>% paste(collapse="")]] <- get(load(paste0("Mus musculus/~organism part + developmental stage/",filename)))$rnaseq
}

data<-NULL;batch<-NULL;tissues<-NULL; for(i in experiments %>% seq_along){
  data%<>%cbind(experiments[[i]] %>% assays %$% counts)
  batch%<>%c(names(experiments)[[i]] %>% rep(dim(experiments[[i]])[2]))
  tissues[[i]]<-(experiments[[i]]$organism_part)
}; batch%<>%factor;tissue<-Reduce(c,tissues)

# intersections<-NULL; for(i in 1:5){
#   for(j in 1:5){
#     for(k in 1:5){
#       intersections%<>%c(tissues[[i]] %>% intersect(tissues[[j]]) %>% intersect(tissues[[k]]) %>% length)
#     }
#   }
# };intersections%<>%array(5 %>% rep(3))

common.tissues<-Reduce(intersect,tissues)
data%<>%extract(,tissue%in%common.tissues)
batch%<>%extract(tissue%in%common.tissues)
tissue%<>%extract(tissue%in%common.tissues)

#merging replicates
data[,5] %<>% add(data[,6]) %<>% divide_by(2); data%<>%extract(,-6); batch%<>%extract(-6); tissue%<>%extract(-6)

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
RUVg(filtered,
)

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
