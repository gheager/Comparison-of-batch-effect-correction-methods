library(magrittr)
library(parallel)
library(doSNOW)
library(gridExtra)

gPCA<-function(data,batch,center=TRUE,scale=FALSE,log=FALSE,scaleY=FALSE,nperm=0){
  if(log) data %<>% subtract(min(.)) %<>% log1p
  X<-data %>% t %>% scale(center,scale) %>% t %>% na.omit %>% t
  batch %<>% factor
  Y<-batch %>% unique %>% sapply(function(.)batch==.)
  if(scaleY) Y %<>% t %<>% divide_by(colSums(t(.))) %<>% t
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
  if(nperm!=0){
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
    delta_perm %<>% divide_by(var(X%*%sv$v[,1])[1])
    PCu<-var(sv$u[,1]*sv$d[1])/sum(diag(var(sv$u*diag(sv$d))))
    p.value<-mean(delta<delta_perm)
    'p-value :' %>% paste(p.value,'\n') %>% cat
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
    'data'=data,
    'batch'=batch,
    'na.omit'=X %>% attr('na.action')
  ))
}

viz_gpca<-function(gpca,dims=1:2,guided=TRUE){
  pca<-if(guided) gpca$gpca else gpca$pca
  ggplot()+aes(x=pca$u[,dims[1]]*pca$d[dims[1]],y=pca$u[,dims[2]]*pca$d[dims[2]],colour=gpca$batch)+
    geom_point()+stat_ellipse()+
    xlab(if(guided) paste0('gPC',dims[1],'~PC',gpca$variance.ranks[dims[1]]) else paste0('PC',dims[1]))+
    ylab(if(guided) paste0('gPC',dims[2],'~PC',gpca$variance.ranks[dims[2]]) else paste0('PC',dims[2]))+
    labs(colour='batch')
}

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
  grid.arrange(ranks.plot,cumulative.plot,ncol=2)
}

viz_gpca_pvalue<-function(gpca){
  if(is.null(gpca$p.value)){
    stop('No p-value computed')
  }else{
    ggplot()+
      geom_density(aes(x=gpca$delta_perm),colour='black',fill='black')+
      geom_point(aes(x=gpca$delta,y=0),colour='red')+
      geom_text(aes(x=gpca$delta,y=.1,label=gpca$p.value),colour='red')
  }
}

compare.pca<-function(corrected,raw,batch,tissue,guided=FALSE){
  if(guided){
    raw.pca<-raw$gpca
    corrected.pca<-corrected$gpca
  }else{
    raw.pca<-raw$pca
    corrected.pca<-corrected$pca
  }
  grid.arrange(ncol=2,nrow=2,
               raw %>% viz_gpca(guided=guided) + geom_line(aes(group=tissue),colour='grey') + theme(legend.position = 'none'),
               ggplot(mapping=aes(
                 x=t(corrected$data)%*%raw.pca$v[,1],
                 y=t(corrected$data)%*%raw.pca$v[,2],
                 colour=batch
               ))+geom_point()+stat_ellipse() + geom_line(aes(group=tissue),colour='grey') + theme(legend.position = 'none'),
               ggplot(mapping=aes(
                 x=t(raw$data)%*%corrected.pca$v[,1],
                 y=t(raw$data)%*%corrected.pca$v[,2],
                 colour=batch
               ))+geom_point()+stat_ellipse() + geom_line(aes(group=tissue),colour='grey') + theme(legend.position = 'none'),
               corrected %>% viz_gpca(guided=guided) + geom_line(aes(group=tissue),colour='grey') + theme(legend.position = 'none')
  )
}
