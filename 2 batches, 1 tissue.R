t<-list(
  mtab2836=get(load('data/E-MTAB-2836-atlasExperimentSummary.Rdata')),
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

filtered%<>%extract(,tissue=='heart')
x%<>%extract(tissue=='heart')
tissue%<>%extract(tissue=='heart')
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

