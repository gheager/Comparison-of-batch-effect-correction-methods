library(ExpressionAtlas)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(parallel)
library(doSNOW)

### data import from previous savings
methods<-c(
  'filtered',
  'pam',
  'combat',
  'ccombat',
  'rg',
  'ra',
  'ruvg'
)
setwd('corrected data')
for(method in methods) eval(parse(text=paste0(
  method,"<-get(load('",method,".Rdata'))"
)))

x<-as.factor(sapply(strsplit(colnames(filtered),split="_"),"[",2))#experiments
tissue<-as.factor(sapply(strsplit(colnames(filtered),split="_"),"[",3))


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
