library(ExpressionAtlas)
setwd('C:/Users/guill/Desktop/RUVg tests')
source("generic_functions.R")
library(magrittr)
library(factoextra)
library(FactoMineR)
library(parallel)
library(doSNOW)
## get gtex expreriments from exression atlas
#atlasData <- getAtlasData(c("E-MTAB-5214", "E-MTAB-513", "E-MTAB-2836", "E-MTAB-3716", "E-MTAB-4344"))
#previous method doesn't work for me
atlasData=list(
  #'mtab5214'=get(load("E-MTAB-5214-atlasExperimentSummary.Rdata")),
  'mtab513'=get(load("E-MTAB-513-atlasExperimentSummary.Rdata")),
  'mtab2836'=get(load("E-MTAB-2836-atlasExperimentSummary.Rdata")),
  'mtab3716'=get(load("E-MTAB-3716-atlasExperimentSummary.Rdata")),
  'mtab4344'=get(load("E-MTAB-4344-atlasExperimentSummary.Rdata")),
  'mtab3871'=get(load("E-MTAB-3871-atlasExperimentSummary.Rdata"))
)
t<-atlasData
all <- vector()
#################################################
for (i in names(t)) {
  expAcc <- i
  k <- t[[i]]
  exp <- k$rnaseq
  eCounts <- assays(exp)$counts
  colnames(eCounts) %<>% paste0('_') %<>% paste0(i) %<>% paste0('_') %<>% paste0(exp %>% colData %>% use_series('organism_part'))
  samples<-colnames(eCounts)
  all <- cbind(all,eCounts)
}
## filtering low expression signals
filter <- rowSums(all>10)>=15
filtered <- all[filter,]
filterCols <- colSums(filtered == 0) / nrow(filtered) < 0.90
filtered <- filtered[,filterCols]

methods<-c(
  'filtered',
  'pam',
  'combat',
  'ccombat',
  'rg',
  'ra',
  'ruvg'
)
setwd('c:/users/guill/desktop/ruvg tests/batch_effect_data')
for(method in methods) eval(parse(text=paste0(
  method,"<-get(load('",method,".Rdata'))"
)))

x<-as.factor(sapply(strsplit(colnames(filtered),split="_"),"[",2))#experiments
tissue<-as.factor(sapply(strsplit(colnames(filtered),split="_"),"[",3))

filtered %>% save(file='filtered.Rdata')
#raw data
pca_raw<-PCA(t(filtered) %>% add(1) %>% log,graph=FALSE)
pdf('gtexw_raw.pdf', width=18, height=18)
fviz_pca_ind(pca_raw,col.ind=x,geom='point')
fviz_pca_ind(pca_raw,col.ind=x,geom='point',addEllipses = TRUE)
dev.off()

# Y<-x %>% unique %>% sapply(function(.)x==.)
# pca<-PCA(t(filtered),graph=FALSE)
# gpca<-PCA(t(Y)%*%t(filtered),graph=FALSE)
# var(t(filtered)%*%gpca$svd$V[,1])/var(t(filtered)%*%pca$svd$V[,1])
# gp<-gPCA.batchdetect(t(filtered),x %>% as.numeric,nperm=100,scaleY=TRUE)
# gp$delta

# X<-filtered %>% t %>% scale
# Y<-x %>% unique %>% sapply(function(.)x==.)
# Ym<-t(t(Y)/colSums(Y))
# sv<-svd(X)
# gsv<-svd(t(Ym)%*%X)
# #PCA on the whole data
# ggplot()+aes(x=sv$u[,1],y=sv$u[,2],colour=x)+geom_point()+stat_ellipse()+
#   labs(title='PCA on the whole dataset',x='PC1 (eigengene 1)',y='PC2 (eigengene 2)')
# #PCA on batch-aggregated data (gPCA)
# ggplot()+aes(x=gsv$u[,1],y=gsv$u[,2],colour=x %>% unique)+geom_point(size=5)+
#   labs(title='PCA on the batch-aggregated data (gPCA)',x="PC'1",y="PC'2")
# #projection of the whole data on the PC of batch-aggregated data
# ggplot()+
#   aes(x=X%*%gsv$v[,1]/gsv$d[1],
#       y=X%*%gsv$v[,2]/gsv$d[2],
#       colour=x)+geom_point()+stat_ellipse()+
#   geom_point(aes(x=gsv$u[,1],y=gsv$u[,2],colour=x %>% unique),size=5)+
#   labs(title='Projection of the dataset on the batch-aggregated eigengenes',x="PC'1",y="PC'2")
# 
# var(X%*%gsv$v[,1])/var(X%*%sv$v[,1])


#pamr
library(pamr)
pamobject<-list(
  'x'=filtered,
  'batchlabels'=x
)
pam <- pamr.batchadjust(pamobject)$x
#pamset <- newSeqExpressionSet(pam %>% floor, phenoData = data.frame(tissue, row.names=colnames(filtered)))
pca_pam<-PCA(t(pam) %>% add(1-pam %>% min) %>% log,graph=FALSE)
pdf('gtexw_after_pamr.pdf', width=18, height=18)
fviz_pca_ind(pca_pam,col.ind=x,geom='point')
fviz_pca_ind(pca_pam,col.ind=x,geom='point',addEllipses = TRUE)
dev.off()
pca_pam<-PCA(pam %>% t,graph=FALSE)
gpca_pam<-PCA(t(Y)%*%t(pam),graph=FALSE)
var(t(pam)%*%gpca_pam$svd$V[,1])/var(t(pam)%*%pca_pam$svd$V[,1])
gp_pam<-gPCA.batchdetect(t(pam),x %>% as.numeric,nperm=100)
gp_pam$delta


#ComBat
library(sva)
combat<-ComBat(filtered,x)
pca_combat<-PCA(t(combat) %>% add(1-combat %>% min) %>% log,graph=FALSE)
pdf('gtexw_after_combat.pdf',width=18,height=18)
fviz_pca_ind(pca_combat,col.ind=x,geom='point')
fviz_pca_ind(pca_combat,col.ind=x,geom='point',addEllipses = TRUE)
dev.off()
pca_combat<-PCA(combat %>% t,graph=FALSE)
gpca_combat<-PCA(t(Y)%*%t(combat),graph=FALSE)
var(t(combat)%*%gpca_combat$svd$V[,1])/var(t(combat)%*%pca_combat$svd$V[,1])
gp_combat<-gPCA.batchdetect(t(combat),x %>% as.numeric,nperm=100)
gp_combat$delta

#ComBat with regression on tissues
library(sva)
ccombat<-ComBat(filtered, batch=x, mod=model.matrix(~tissue %>% factor), prior.plots = TRUE)
pca_ccombat<-PCA(t(ccombat),graph=FALSE)
log_pca_ccombat<-PCA(t(ccombat) %>% add(1-ccombat %>% min) %>% log,graph=FALSE)
pdf('gtex_after_ccombat_with_outcome.pdf',width=18,height=18)
fviz_pca_ind(pca_ccombat,col.ind=x,geom='point')
fviz_pca_ind(log_pca_ccombat,col.ind=x,geom='point')
fviz_pca_ind(pca_ccombat,col.ind=x,geom='point',addEllipses = TRUE)
fviz_pca_ind(log_pca_ccombat,col.ind=x,geom='point',addEllipses = TRUE)
fviz_pca_ind(log_pca_ccombat,col.ind=x,geom='text',addEllipses = TRUE)
dev.off()
pca_ccombat<-PCA(ccombat %>% t,graph=FALSE)
gpca_ccombat<-PCA(t(Y)%*%t(ccombat),graph=FALSE)
var(t(ccombat)%*%gpca_ccombat$svd$V[,1])/var(t(ccombat)%*%pca_ccombat$svd$V[,1])
gp_ccombat<-gPCA.batchdetect(t(ccombat),x %>% as.numeric,nperm=100)
gp_ccombat$delta

#ratio G
library(bapred)
rg<-ratiog(filtered %>% t, x %>% as.numeric %>% factor)$xadj %>% t
pca_rg<-PCA(t(rg) %>% add(1-rg %>% min) %>% log,graph=FALSE)
pdf('gtexw_after_ratiog.pdf',width=18,height=18)
fviz_pca_ind(pca_rg,col.ind=x,geom='point')
fviz_pca_ind(pca_rg,col.ind=x,geom='point',addEllipses = TRUE)
fviz_pca_ind(pca_rg,col.ind=x,geom='text',addEllipses = TRUE)
dev.off()
pca_rg<-PCA(rg %>% t,graph=FALSE)
gpca_rg<-PCA(t(Y)%*%t(rg),graph=FALSE)
var(t(rg)%*%gpca_rg$svd$V[,1])/var(t(rg)%*%pca_rg$svd$V[,1])
gp_rg<-gPCA.batchdetect(t(rg),x %>% as.numeric,nperm=100)
gp_rg$delta

#ratio A
library(bapred)
ra<-ratioa(filtered %>% t, x %>% as.numeric %>% factor)$xadj %>% t
pca_ra<-PCA(t(ra) %>% add(1-ra %>% min) %>% log,graph=FALSE)
pdf('gtexw_after_ratioa.pdf',width=18,height=18)
fviz_pca_ind(pca_ra,col.ind=x,geom='point')
fviz_pca_ind(pca_ra,col.ind=x,geom='point',addEllipses = TRUE)
dev.off()
pca_ra<-PCA(ra %>% t,graph=FALSE)
gpca_ra<-PCA(t(Y)%*%t(ra),graph=FALSE)
var(t(ra)%*%gpca_ra$svd$V[,1])/var(t(ra)%*%pca_ra$svd$V[,1])
gp_ra<-gPCA.batchdetect(t(ra),x %>% as.numeric,nperm=100)
gp_ra$delta

#RUVg
## coefficient of Variation
co.var <- function(x) ( 100*apply(x,1,sd)/rowMeans(x) ) 
## coefficient of variation across all the sammples
cov.allGenes<-na.omit(co.var(as.matrix(filtered)))
# Using threshold (1% qauntile) for identification of non vriable genes that has least cov.
## identyfing number of genes changed acros several cov thresholds range
cov.range<-seq(range(cov.allGenes)[1], range(cov.allGenes)[2], by = 10)
ngenes<-matrix()
for (i in seq_along(cov.range)){
  ngenes[i]<- sum((cov.allGenes<=cov.range[i])*1)
}
# Using threshold (1% qauntile) for identification of non vriable genes that has least cov.
leastVar.genes<-rownames(as.matrix(sort(cov.allGenes[cov.allGenes < quantile(cov.range, c(.01))[[1]]])))
library(RUVSeq)
ruvg_<-RUVg(filtered,leastVar.genes,1)
ruvg<-ruvg_ %>% use_series(normalizedCounts)
ruvg %>% save(file='RUVg_data.Rdata')


### batch effect detection in the datasets with gPCA (guided PCA)

cumulative_gpca <- function(data,batch){
  X<-data %>% t %>% scale %>% t %>% na.omit %>% t
  batch%<>%factor
  'Computing PCA\n' %>% cat
  sv<-svd(X)
  'Computing gPCA\n' %>% cat
  Y<-batch %>% unique %>% sapply(function(.)batch==.) %>% t %>% divide_by(colSums(t(.))) %>% t
  gsv<-svd(t(Y)%*%X)
  PCu<-sv$d[1]^2/sum(sv$d^2)
  PCg<-var(X%*%gsv$v[,1])[1]/sum(diag(var(X%*%gsv$v)))
  (PCg-PCu)/PCg
}

compare_gpca<-function(data,batch,center=TRUE,scale=FALSE,log=FALSE,nperm=0,filename=NULL,filename_prefix=NULL,path=getwd()){
  X<-data %>% t %>% scale(center=center,scale=scale) %>% t %>% na.omit %>% t
  batch %<>% factor
  if(log) X %<>% subtract(X %>% min) %<>% add(1) %<>% log
  Y<-batch %>% unique %>% sapply(function(.)batch==.)
  Ym<-t(t(Y)/colSums(Y))
  'Computing PCA\n' %>% cat
  sv<-svd(X)
  'Computing gPCA\n' %>% cat
  gsv<-svd(t(Ym)%*%X)
  delta<-var(X%*%gsv$v[,1])[1]/var(X%*%sv$v[,1])[1] #var(sv$u[,1]*sv$d[1])#sum(colVars(X%*%sv$v))
  'delta statistic :' %>% paste(delta,'\n') %>% cat
  if(!(filename %>% is.null)){
    'Printing PCA plots\n' %>% cat
    #graphs printing in a pdf file
    path %>% paste0('/',filename_prefix,'_',filename,'.pdf') %>% pdf
    #pdf(getwd() %>% paste0('/',filename_prefix,filename,'.pdf'),width=18,height=18)
    #PCA on the whole data
    print(ggplot()+aes(x=sv$u[,1],y=sv$u[,2],colour=batch)+geom_point()+stat_ellipse()+
            labs(title='PCA on the whole dataset',x='PC1 (eigengene 1)',y='PC2 (eigengene 2)'))
    #PCA on batch-aggregated data (gPCA)
    print(ggplot()+aes(x=gsv$u[,1],y=gsv$u[,2],colour=batch %>% unique)+geom_point(size=5)+
            labs(title='PCA on the batch-aggregated data (gPCA)',x="PC'1",y="PC'2"))
    #projection of the whole data on the PCs of batch-aggregated data
    print(ggplot()+
            aes(x=X%*%gsv$v[,1]/gsv$d[1],
                y=X%*%gsv$v[,2]/gsv$d[2],
                colour=batch)+geom_point()+stat_ellipse()+
            geom_point(aes(x=gsv$u[,1],y=gsv$u[,2],colour=batch %>% unique),size=5)+
            labs(title='Projection of the dataset on the batch-aggregated eigengenes',x="PC'1",y="PC'2"))
    dev.off()
    'PCA plots printed at ' %>% paste0(path,'/',filename_prefix,'_',filename,'.pdf\n') %>% message
  }
  #estimation of the p-value
  'Estimating p-value\n' %>% cat
  cl <- detectCores() %>% subtract(1) %>% makeSOCKcluster
  cl %>% clusterExport(c('%>%','colVars'))
  cl %>% registerDoSNOW
  pb <- txtProgressBar(min=1, max=nperm, style=3)
  delta_perm <- foreach(i=seq_len(nperm),
                        .options.snow=list(progress=function(n)setTxtProgressBar(pb,n)),
                        .combine='c') %dopar% {
                          Y<-batch %>% unique %>% sapply(function(.)batch %>% sample==.)
                          Ym<-t(t(Y)/colSums(Y))
                          gsv<-svd(t(Ym)%*%X)
                          var(X%*%gsv$v[,1])[1]
                        }
  pb %>% close
  cl %>% stopCluster
  #normalization
  #delta_perm %<>% divide_by(sum(colVars(sv$u%*%diag(sv$d))))
  delta_perm %<>% divide_by(var(X%*%sv$v[,1])[1]) #divide_by(var(sv$u[,1]*sv$d[1]))
  
  # delta_perm<-NULL
  # delta_perm<-(1:nperm) %>% parSapply(cl,.,
  #                                     function(i){
  #                                       Y<-batch %>% unique %>% sapply(function(.)batch %>% sample==.)
  #                                       Ym<-t(t(Y)/colSums(Y))
  #                                       gsv<-svd(t(Ym)%*%X)
  #                                       cat(i)
  #                                       return(colVars(X%*%gsv$v[,1])/sum(colVars(X%*%sv$v)))
  #                                     })
  # for(i in 1:nperm){
  # Y<-batch %>% unique %>% sapply(function(.)batch %>% sample==.)
  # Ym<-t(t(Y)/colSums(Y))
  # gsv<-svd(t(Ym)%*%X)
  # delta_perm %<>% c(colVars(X%*%gsv$v[,1])/sum(colVars(X%*%sv$v)))
  # cat(i)
  # }
  PCu<-var(sv$u[,1]*sv$d[1])/sum(diag(var(sv$u*diag(sv$d))))
  p.value<-mean(delta<delta_perm)
  return(list(
    'delta'=delta,
    'p.value'=p.value,
    'delta_perm'=delta_perm
  ))
  #return((var(X%*%gsv$v[,1])/var(X%*%sv$v[,1]))[1])
}

# library(parallel)
# cl<-makeCluster(detectCores()-1)
# clusterExport(cl,c('add','subtract','%>%','%<>%','colVars'))
# clusterExport(cl,objects())
# #clusterExport(cl,'compare_gpca')

methods<-c(
  'filtered',
  'pam',
  'combat',
  'ccombat',
  'rg',
  'ra',
  'ruvg'
)
for(method in methods) method %>% get %>% save(file=method %>% paste0('_data.Rdata'))
gpca_comparison <- methods %>% sapply(
  function(data)compare_gpca(data %>% get,x,filename=data,filename_prefix='gPCA',scale=FALSE,log=FALSE,nperm=100)
)
gpca_comparison_scaled <- methods %>% sapply(
  function(data)compare_gpca(data %>% get,x,scale=TRUE,log=FALSE,nperm=100)
)
gpca_comparison_logged <- methods %>% sapply(
  function(data)compare_gpca(data %>% get,x,filename=data,filename_prefix='gPCAlogged',scale=TRUE,log=TRUE,nperm=100)
)

'gPCA.pdf' %>% pdf
for(i in 1:7) gpca_comparison_scaled[,i] %>% view_delta_distribution %>% add(ggtitle(methods[i])) %>% print
dev.off()

sv<-svd(filtered %>% t %>% scale)
vars<-1:324 %>% sapply(
  function(i)var(sv$u[,i]*sv$d[i])
)
vars%<>%divide_by(sum(.))
vars%<>%divide_by(vars[1])
# save(gpca_comparison, file='gpca_comparison.Rdata')
# save(gpca_comparison_scaled, file='gpca_comparison_scaled.Rdata')
# save(gpca_comparison_logged, file='gpca_comparison_logged.Rdata')
# 
# write.table(gpca_comparison,file="gpca_comparison.txt",sep="\t")
# write.table(gpca_comparison_scaled,file="gpca_comparison_scaled.txt",sep="\t")
# write.table(gpca_comparison_logged,file="gpca_comparison_logged.txt",sep="\t")
# 
# ggplot()+aes(x=1,y=sv$v[,1])+geom_bar(stat='identity')
# sv$v[,1] %>% order %>% {colnames(X)[.[1:10]]}
# 
# 
# 
# nobatch_data<-filtered[,x=='mtab2836']
# fake_batch<-rbinom(n=ncol(nobatch_data),size=1,prob=.5)
# test_gpca_nobatch<-compare_gpca(nobatch_data,fake_batch,nperm=20,filename='fake',filename_prefix='gpca_')
# compare_gpca(filtered,x,'test','gpca_')
# 
# compare_gpca(filtered,x,'test_filtered','gpca_',nperm=100)
# 

setwd('batch_effect_comparison')
gpca_comparison_logged <- get(load('gpca_comparison_logged.Rdata'))
'delta distributions logged.pdf' %>% pdf
for(i in gpca_comparison_logged %>% seq_along){
  print(ggplot()+geom_density(aes(x=gpca_comparison_logged[[i]]$delta_perm),fill='black')+
          geom_label(aes(label=gpca_comparison_logged[[i]]$p.value,x=gpca_comparison_logged[[i]]$delta,y=0),colour='red')+
          labs(title=methods[i])+
          geom_line(aes(x=seq(gpca_comparison_logged[[i]]$delta_perm %>% min,gpca_comparison_logged[[i]]$delta_perm %>% max,length.out=100),
                        y=dnorm(seq(gpca_comparison_logged[[i]]$delta_perm %>% min,gpca_comparison_logged[[i]]$delta_perm %>% max,length.out=100),
                                mean=gpca_comparison_logged[[i]]$delta_perm %>% mean,
                                sd=gpca_comparison_logged[[i]]$delta_perm %>% sd)),
                    colour='red'))
}
dev.off()

gpca_comparison_scaled <- get(load('gpca_comparison_scaled.Rdata'))
'delta distributions scaled.pdf' %>% pdf
for(i in gpca_comparison_scaled %>% seq_along){
  print(ggplot()+geom_density(aes(x=gpca_comparison_scaled[[i]]$delta_perm),fill='black')+
          geom_label(aes(label=gpca_comparison_scaled[[i]]$p.value,x=gpca_comparison_scaled[[i]]$delta,y=0),colour='red')+
          labs(title=methods[i]))
}
dev.off()

gpca_comparison <- get(load('gpca_comparison.Rdata'))
'delta distributions (non scaled non logged).pdf' %>% pdf
for(i in gpca_comparison %>% seq_along){
  print(ggplot()+geom_density(aes(x=gpca_comparison[[i]]$delta_perm),fill='black')+
          geom_label(aes(label=gpca_comparison[[i]]$p.value,x=gpca_comparison[[i]]$delta,y=0),colour='red')+
          labs(title=methods[i]))
}
dev.off()


ggplot()+stat_density(aes(x=gpca_comparison_logged[[1]]$delta_perm),fill='black')

X<-filtered %>% t %>% scale %>% t %>% na.omit %>% t
X %<>% subtract(X %>% min) %<>% add(1) %<>% log
sv<-X %>% svd
variance.parts<-sv$d %>% seq_along %>% sapply(function(i)var(sv$u[,i]*sv$d[i])) %>% divide_by(sum(.))

normality <- gpca_comparison_logged %>% seq_along %>% sapply(function(i)gpca_comparison_logged[[i]]$delta_perm %>% shapiro.test %>% use_series(p.value)); names(normality)<-methods

#test with gaussian variables
batch<-rbinom(100,1,.5)
gauss<-NULL
for(i in seq_len(1000)) gauss %<>% rbind(rnorm(100,mean=ifelse(batch,0,5),sd=ifelse(batch,1,1)))
a<-compare_gpca(gauss,batch,nperm=1000)
a %>% view_delta_distribution

b<-gPCA.batchdetect(gauss %>% t, batch, scaleY)

view_delta_distribution <- function(gpca){
  ggplot()+geom_density(aes_string(x=gpca$delta_perm),colour='blue',fill='blue')+
    geom_label(aes_string(x=gpca$delta, y=0, label=gpca$p.value), colour='blue')
}

setwd('c:/users/guill/desktop/ruvg tests/comparisons')
comparisons<-get(load('gpca_comparison_logged.Rdata'))
'delta distributions.pdf' %>% pdf
for(i in comparisons %>% seq_along){
  comparisons[[i]] %>% view_delta_distribution %>% add(ggtitle(methods[i])) %>% print
}
dev.off()


#correlation between inertia and delta p-value
pvals<-comparisons %>% sapply(function(gpca)gpca$p.value)
deltas<-comparisons %>% sapply(function(gpca)gpca$delta)
plot(inertias,pvals)
plot(inertias,deltas)
