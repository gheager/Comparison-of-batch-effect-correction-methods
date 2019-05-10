library(ExpressionAtlas)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(parallel)
library(doSNOW)
library(gridExtra)

### data import from previous savings
methods<-c(
  'raw',
  'pam',
  'combat',
  'ccombat',
  'rg',
  'ra',
  'ruvg'
)
setwd('c:/users/medion/desktop/Comparison-of-batch-effect-correction-methods/corrected data')
for(method in methods) eval(parse(text=paste0(
  method,"<-get(load('",method,".Rdata'))"
)))

x<-as.factor(sapply(strsplit(colnames(raw),split="_"),"[",2))#experiments
tissue<-as.factor(sapply(strsplit(colnames(raw),split="_"),"[",3))


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

gPCA<-function(data,batch,center=TRUE,scale=FALSE,log=FALSE,scaleY=FALSE,nperm=0){
  if(log) data %<>% subtract(min(.)) %<>% add(1) %<>% log1p
  X<-data %>% t %>% scale(center,scale) %>% t %>% na.omit %>% t
  batch %<>% factor
  Y<-batch %>% unique %>% sapply(function(.)batch==.)
  if(scaleY) Y %<>% t %<>% divide_by(colSums(t(.))) %<>% t #Ym<-t(t(Y)/colSums(Y))
  'Computing PCA\n' %>% cat
  sv<-svd(X)
  'Computing gPCA\n' %>% cat
  gsv<-svd(t(Y)%*%X)
  gpca<-list(
    u=X%*%gsv$v %*% diag(1/gsv$d),
    v=gsv$v,
    d=gsv$d
  )
  PC.variances<-colVars(sv$u%*%diag(sv$d))
  gPC.variances<-colVars(X%*%gsv$v)
  variance.part<-gPC.variances[1]/sum(PC.variances)
  'part of variance from gPC1 :' %>% paste(variance.part,'\n') %>% cat
  delta<-gPC.variances[1]/PC.variances[1]
  'delta statistic :' %>% paste(delta,'\n') %>% cat
  cumdelta<-cumsum(gPC.variances)/cumsum(PC.variances[gPC.variances %>% seq_along])
  'cumulative delta statistics :\n' %>% cat
  'delta_' %>% paste0(cumdelta %>% seq_along,'=',cumdelta,'\n') %>% cat
  variance.ranks<-gPC.variances %>% sapply(function(v)sum(v<PC.variances))
  'variance ranks :' %>% paste(variance.ranks,'\n') %>% cat
  #variance.ranks<-gPC.variances %>% cumsum %>% sapply(function(v)sum(v<(PC.variances %>% rev %>% cumsum %>% rev))+1)
  if(nperm!=0){
    #estimation of the p-value
    'Estimating p-value\n' %>% cat
    cl <- detectCores() %>% subtract(1) %>% makeSOCKcluster
    cl %>% clusterExport(c('%>%','%<>%','extract','colVars'))
    cl %>% registerDoSNOW
    pb <- txtProgressBar(min=1, max=nperm, style=3)
    delta_perm <- foreach(i=seq_len(nperm),
                          .options.snow=list(progress=function(n)setTxtProgressBar(pb,n)),
                          .combine='c') %dopar% {
                            Y %<>% extract(Y %>% nrow %>% seq_len %>% sample,)
                            batch %>% unique %>% sapply(function(.)batch %>% sample==.)
                            gsv<-svd(t(Y)%*%X)
                            var(X%*%gsv$v[,1])[1]
                          }
    pb %>% close
    cl %>% stopCluster
    #normalization
    #delta_perm %<>% divide_by(sum(colVars(sv$u%*%diag(sv$d))))
    delta_perm %<>% divide_by(var(X%*%sv$v[,1])[1]) #divide_by(var(sv$u[,1]*sv$d[1]))
    
    PCu<-var(sv$u[,1]*sv$d[1])/sum(diag(var(sv$u*diag(sv$d))))
    p.value<-mean(delta<delta_perm)
  }
  return(list(
    'variance.part'=variance.part,
    'PC.variances'=PC.variances,
    'gPC.variances'=gPC.variances,
    'delta'=delta,
    'cumdelta'=cumdelta,
    'variance.ranks'=variance.ranks,
    'p.value'='p.value' %>% get0,
    'delta_perm'='delta_perm' %>% get0,
    'pca'=sv,
    'batch.pca'=gsv,
    'gpca'=gpca,
    'batch'=batch
  ))
}

viz_gpca<-function(gpca,dims=1:2,guided=TRUE){
  pca<-if(guided) gpca$gpca else gpca$pca
  ggplot()+aes(x=pca$u[,dims[1]],y=pca$u[,dims[2]],colour=gpca$batch)+
    geom_point()+stat_ellipse()+
    xlab(if(guided) paste0('gPC',dims[1],'~PC',gpca$variance.ranks[dims[1]]) else paste0('PC',dims[1]))+
    ylab(if(guided) paste0('gPC',dims[2],'~PC',gpca$variance.ranks[dims[2]]) else paste0('PC',dims[2]))+
    labs(colour='batch')
}

# viz_gpca_contrib<-function(gpca,transformation=identity,end='max'){
#   end%<>%as.character%<>%switch(max=gpca$variance.ranks %>% max+1, all=gpca$PC.variances %>% length, end %>% as.numeric)
#   ranks.plot<-ggplot()+
#     geom_bar(aes_string(y=gpca$PC.variances[1:end] %>% transformation,x=1:end),stat='identity',width=1)+
#     geom_bar(aes_string(y=gpca$gPC.variances %>% transformation,x=gpca$variance.ranks,fill=gpca$gPC.variances %>% seq_along %>% factor),stat='identity',width=1,position='dodge')+
#     xlim(c(0,end+1))+
#     theme(legend.position='none')
#   endc<-gpca$gPC.variances %>% length
#   cumulative.plot<-ggplot()+
#     geom_bar(aes(y=gpca$PC.variances[1:endc] %>% cumsum %>% transformation,x=1:endc),stat='identity',width=1)+
#     geom_bar(aes(y=gpca$gPC.variances %>% cumsum %>% transformation,x=gpca$gPC.variances %>% seq_along,fill=gpca$gPC.variances %>% seq_along %>% factor),stat='identity',width=1,position='dodge')+
#     theme(legend.position='none')
#   grid.arrange(ranks.plot,cumulative.plot,ncol=2)
# }

viz_gpca_contrib<-function(gpca,transformation=identity,end='max'){
  end%<>%as.character%<>%switch(max=gpca$variance.ranks %>% max+1, all=gpca$PC.variances %>% length, end %>% as.numeric)
  ranks.plot<-ggplot()+
    geom_bar(aes_string(y=gpca$PC.variances[1:end] %>% transformation,x=1:end),stat='identity',width=1)+
    geom_bar(aes_string(y=gpca$gPC.variances %>% transformation,x=gpca$variance.ranks,fill=gpca$gPC.variances %>% seq_along %>% factor),stat='identity',width=1,position='dodge')+
    xlim(c(0,end+1))+
    theme(legend.position='none')+xlab('PCs')+ylab('Parts of variance')
  endc<-gpca$gPC.variances %>% length
  cumulative.plot<-ggplot(mapping=aes(x=gpca$gPC.variances %>% seq_along %>% factor))+
    geom_bar(aes(y=gpca$PC.variances[1:endc] %>% cumsum %>% transformation),stat='identity',width=1)+
    geom_bar(aes(y=gpca$gPC.variances %>% cumsum %>% transformation,fill=gpca$gPC.variances %>% seq_along %>% factor),stat='identity',width=1,position='dodge')+
    theme(legend.position='none')+xlab('PCs')+ylab('Cumulative parts of variance')+
    geom_text(aes(label=gpca$cumdelta %>% round(2),
                  y=gpca$gPC.variances %>% cumsum %>% transformation %>% divide_by(2)))
  # cumulative.plot<-ggplot()+
  #   geom_bar(aes(y=gpca$PC.variances[1:endc] %>% cumsum %>% transformation,x=1:endc),stat='identity',width=1)+
  #   geom_bar(aes(y=gpca$gPC.variances %>% cumsum %>% transformation,x=gpca$gPC.variances %>% seq_along,fill=gpca$gPC.variances %>% seq_along %>% factor),stat='identity',width=1,position='dodge')+
  #   theme(legend.position='none')
  grid.arrange(ranks.plot,cumulative.plot,ncol=2)
}

gPCAs<-methods %>% lapply(function(method)method %>% get %>% gPCA(x))

setwd('~/Comparison-of-batch-effect-correction-methods')
'gPCA contributions.pdf' %>% pdf
for(i in gPCAs %>% seq_along){
  gPCAs[[i]] %>% viz_gpca_contrib %>% add(ggtitle(methods[i])) %>% print
}
dev.off()

'gPCA log10 contributions.pdf' %>% pdf
for(i in gPCAs %>% seq_along){
  gPCAs[[i]] %>% viz_gpca_contrib(log10) %>% add(ggtitle(methods[i])) %>% print
}
dev.off()


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


gpca_comparison <- methods %>% sapply(
  function(data)compare_gpca(data %>% get,x,filename=data,filename_prefix='gPCA',scale=FALSE,log=FALSE,nperm=100)
)
gpca_comparison_scaled <- methods %>% sapply(
  function(data)compare_gpca(data %>% get,x,scale=TRUE,log=FALSE,nperm=100)
)
gpca_comparison_logged <- methods %>% sapply(
  function(data)compare_gpca(data %>% get,x,filename=data,filename_prefix='gPCAlogged',scale=TRUE,log=TRUE,nperm=100)
)

setwd(~/Comparison-of-batch-effect-correction-methods)
'gPCA.pdf' %>% pdf
for(i in methods %>% seq_along) gpca_comparison_scaled[,i] %>% view_delta_distribution %>% add(ggtitle(methods[i])) %>% print
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
