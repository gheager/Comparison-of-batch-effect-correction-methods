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

filtered %>% save(file='filtered4344+513+3716.Rdata')
filtered=get(load('filtered4344+513+3716.Rdata'))
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
geneshift.pca(gs,x[.==levels(.)[2]|.==levels(.)[3]],3,.01)->gs
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

gfiltered %>% viz_gpca(guided=FALSE)+geom_line(aes(group=tissue))#+geom_text(aes(label=tissue),size=3)

common.tissues<-intersect(tissue[x=='mtab4344'],tissue[x=='mtab513']) %>% intersect(tissue[x=='mtab3716'])

arw<-filtered
arrows<-NULL
for(h in common.tissues){
  arrows%<>%cbind(filtered[,x=='mtab4344'&tissue==h]-filtered[,x=='mtab513'&tissue==h])
};colnames(arrows)<-common.tissues
arrows %>% rowMeans -> arrow
arw[,x=='mtab513']%<>%add(arrow)

biological.replicates<-NULL
for(h in common.tissues){
  biological.replicates%<>%cbind(rowMeans(filtered[,x=='mtab3716'&tissue==h]))
};colnames(biological.replicates)<-common.tissues

arrows<-NULL
for(h in common.tissues){
  arrows%<>%cbind(filtered[,x=='mtab4344'&tissue==h]-rowMeans(filtered[,x=='mtab3716'&tissue==h]))
};colnames(arrows)<-common.tissues
arrows %>% rowMeans -> arrow
arw[,x=='mtab3716']%<>%add(arrow)


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



batch.inertias<-methods %>% lapply(as_mapper(
  ~batch.inertia(get(.),tissue)
))
volumes<-function(k) batch.inertias %>% sapply(as_mapper(
  ~.$volume.index[k]
)) %>% set_names(methods)
volumes(1)
'volumes.pdf' %>% pdf
for(k in 1:29) plot(log(volumes(k)),main=k) %>% print
dev.off()
