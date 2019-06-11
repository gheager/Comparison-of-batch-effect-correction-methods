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
  mtab513=get(load('data/E-MTAB-513-atlasExperimentSummary.Rdata'))#,
  #mtab3716=get(load('data/E-MTAB-3716-atlasExperimentSummary.Rdata'))
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

filtered %>% save(file='filtered4344+513.Rdata')

#### LOG-TRANSFORMATION
filtered%<>%log1p

filtered %>% gPCA(x,nperm=1000) -> gfiltered
gfiltered %>% viz_gpca_contrib
gfiltered %>% viz_gpca_pvalue

library(pamr)
pamr.batchadjust(list(x=filtered,batchlabels=x))$x->pam
pam %>% gPCA(x,nperm=1000) -> gpam
gpam %>% compare.pca(gfiltered,x)
gpam %>% compare.pca(gfiltered,x,guided=TRUE)
gpam %>% viz_gpca_contrib

library(sva)
ComBat(filtered,x)->combat
combat %>% gPCA(x,nperm=1000) -> gcombat
gcombat %>% compare.pca(gfiltered,x)
gcombat %>% compare.pca(gfiltered,x,guided=TRUE)
gcombat %>% viz_gpca_contrib
gcombat %>% viz_gpca_pvalue

ComBat(filtered,x,mod=model.matrix(~tissue %>% factor)) -> ccombat
ccombat %>% gPCA(x,nperm=1000) -> gccombat
gccombat %>% compare.pca(gfiltered,x)
gccombat %>% compare.pca(gfiltered,x,guided=TRUE)
gccombat %>% viz_gpca_contrib
gccombat %>% viz_gpca_pvalue


ref<-nearest.neighbours.references(filtered,x,2)
ratioA(filtered,x,ref)->ra
ra %>% gPCA(x,nperm=1000) -> gra
gra %>% compare.pca(gfiltered,x)
gra %>% compare.pca(gfiltered,x,guided=TRUE)
gra %>% viz_gpca_contrib
gra %>% viz_gpca_pvalue

ratioG(filtered,x,ref)->rg
rg %>% gPCA(x,nperm=1000) -> grg
grg %>% compare.pca(gfiltered,x)
grg %>% compare.pca(gfiltered,x,guided=TRUE)
grg %>% viz_gpca_contrib
grg %>% viz_gpca_pvalue

geneshift.pca(filtered,x,3,.01)->gs
gs %>% gPCA(x,nperm=1000) -> ggs
ggs %>% compare.pca(gfiltered,x)
ggs %>% compare.pca(gfiltered,x,guided=TRUE)
ggs %>% viz_gpca_contrib
ggs %>% viz_gpca_pvalue

library(harmony)
HarmonyMatrix(filtered,data.frame(x,tissue),vars_use = 'x',do_pca=F,npcs=14,nclust=2)->hm
hm %>% gPCA(x,nperm=1000) -> ghm
ghm %>% compare.pca(gfiltered,x)
ghm %>% compare.pca(gfiltered,x,guided=TRUE)
ghm %>% viz_gpca_contrib
ghm %>% viz_gpca_pvalue



#hm %>% {.%*%diag(gfiltered$pca$d[1:14])} %>% t %>% batch_inertia(tissue)
hm%>% batch_inertia(tissue)
ggplot(mapping=aes(x=-hm[,1],y=-hm[,2],colour=x))+geom_point()+stat_ellipse()+geom_line(aes(group=tissue))

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
batch.inertia(filtered,tissue)$volume.index
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
  geom_line(aes(y=batch.inertia(arw,tissue)$volume.index %>% log),colour='purple')+
  xlim(c(0,28))+ylim(c(-75,0))
ggplot(mapping=aes(x=1:5))+
  geom_point(aes(y=batch.inertia(filtered,tissue)$volume.index[1:5] %>% {-log(.)}),colour='blue')+
  geom_point(aes(y=batch.inertia(pam,tissue)$volume.index[1:5] %>% {-log(.)}),colour='red')+
  geom_point(aes(y=batch.inertia(combat,tissue)$volume.index[1:5] %>% {-log(.)}),colour='green')+
  geom_point(aes(y=batch.inertia(ccombat,tissue)$volume.index[1:5] %>% {-log(.)}),colour='yellow')+
  geom_point(aes(y=batch.inertia(ra,tissue)$volume.index[1:5] %>% {-log(.)}),colour='brown')+
  geom_point(aes(y=batch.inertia(rg,tissue)$volume.index[1:5] %>% {-log(.)}),colour='black')+
  geom_point(aes(y=batch.inertia(gs,tissue)$volume.index[1:5] %>% {-log(.)}),colour='orange')+
  geom_point(aes(y=batch.inertia(hm,tissue)$volume.index[1:5] %>% {-log(.)}),colour='pink')+
  geom_point(aes(y=batch.inertia(arw,tissue)$volume.index[1:5] %>% {-log(.)}),colour='purple')


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

###additional stuff
gfiltered %>% viz_gpca(guided=FALSE,dims=c(1,4))+geom_line(aes(group=tissue))+geom_text(aes(label=tissue),size=3)

angle<-function(u,v){
  u%<>%matrix
  v%<>%matrix
  return(
    (t(u)%*%v/(norm(u,'F')*norm(v,'F')))[1]
  )
}

common.tissues<-intersect(tissue[x=='mtab4344'],tissue[x=='mtab513'])
arrows<-NULL
for(h in common.tissues){
  arrows%<>%cbind(filtered[,x=='mtab4344'&tissue==h]-filtered[,x=='mtab513'&tissue==h])
};colnames(arrows)<-common.tissues
arrows %>% rowMeans -> arrow
arw<-filtered
arw[,x=='mtab513']%<>%add(arrow)
arw %>% gPCA(x,nperm=1000) -> garw
garw %>% viz_gpca+geom_line(aes(group=tissue))+geom_text(aes(label=tissue),size=3)
garw %>% compare.pca(gfiltered,x)
garw %>% compare.pca(gfiltered,x,guided=TRUE)
garw %>% viz_gpca_contrib
garw %>% viz_gpca_pvalue


methods<-c(
  'filtered',
  'pam',
  'combat',
  'ccombat',
  'ra',
  'rg',
  'gs',
  'hm',
  'arw'
)
titles<-c(
  'Before\n correction',
  'BMC',
  'ComBat (without\n biological covariate)',
  'ComBat (with\n biological covariate)',
  'Ratio A',
  'Ratio G',
  'Geneshift',
  'Harmony',
  'Arrow'
)
batch_inertias<-methods %>% sapply(as_mapper(
  ~batch_inertia(get(.),tissue)
))
ggplot()+
  geom_bar(aes(x=methods %>% seq_along %>% factor,y=batch_inertias),stat='identity')+
  scale_x_discrete(labels=titles)+
  xlab('Method')+ggtitle('Intra-class inertia index')+
  theme(axis.text.x = element_text(angle = 90))

ggplot()+
  geom_bar(aes(x=methods %>% seq_along %>% factor,y=log(log(batch_inertias)+1)+1),stat='identity')+
  scale_x_discrete(labels=titles)+
  xlab('Method')+ggtitle('Intra-class inertia index (log scale)')+
  theme(axis.text.x = element_text(angle = 90))

deltas<-methods %>% sapply(as_mapper(
  ~get(paste0('g',.))$delta
))
pvals<-methods %>% sapply(as_mapper(
  ~get(paste0('g',.))$p.value
))
ggplot()+
  geom_bar(aes(x=methods %>% seq_along %>% factor,y=deltas),stat='identity')+
  scale_x_discrete(labels=titles)+
  xlab('Method')+ggtitle('Delta statistic')+
  theme(axis.text.x = element_text(angle = 90))



gpcatissues<-methods %>% lapply(as_mapper(
  ~gPCA(get(.),tissue)
))
tissuesdeltas<-gpcatissues %>% sapply(as_mapper(
  ~.$delta
))
ggplot()+
  geom_bar(aes(x=methods %>% seq_along %>% factor,y=deltas/tissuesdeltas),stat='identity')+
  scale_x_discrete(labels=titles)+
  xlab('Method')+ggtitle('Deltas ratio')+
  theme(axis.text.x = element_text(angle = 90))


batch.inertias<-methods %>% lapply(as_mapper(
  ~batch.inertia(get(.),tissue)
))
lengths<-batch.inertias %>% sapply(as_mapper(
  ~.$volume.index[1]
)) %>% set_names(methods)
areas<-batch.inertias %>% sapply(as_mapper(
  ~.$volume.index[2]
)) %>% set_names(methods)
volumes<-batch.inertias %>% sapply(as_mapper(
  ~.$volume.index[3]
)) %>% set_names(methods)
volumes4<-batch.inertias %>% sapply(as_mapper(
  ~.$volume.index[4]
)) %>% set_names(methods)
volumes<-function(k) batch.inertias %>% sapply(as_mapper(
  ~.$volume.index[k]
)) %>% set_names(methods)

for(k in 1:29) plot(log(volumes(k)),main=k)
