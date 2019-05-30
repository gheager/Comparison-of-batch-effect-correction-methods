library(magrittr)
library(ggplot2)
library(latex2exp)

angle<-function(u,v){
  u%<>%matrix
  v%<>%matrix
  return((t(u)%*%v/(norm(u,'F')*norm(v,'F')))[1] %>% acos %>% divide_by(pi) %>% round(3) %>% paste('$\\pi$'))
}

batchgeneladder<-function(data,batch,mapping=identity,pc=1){
  batch%<>%factor; batch %>% levels -> batches
  data[,batch==batches[1]] %>% t %>% prcomp -> pca1
  data[,batch==batches[2]] %>% t %>% prcomp -> pca2
  data %>% t %>% prcomp -> pca3
  ngenes<-nrow(data)
  ggplot()+aes(
    x=c(-1,-1/3,1/3,1) %>% rep(each=ngenes),
    y=c(
      pca3$rotation[,pc],
      pca1$rotation[,pc],
      pca2$rotation[,pc],
      pca3$rotation[,pc]
    ) %>% sapply(mapping),
    colour=c('',batches[1],batches[2],'') %>% rep(each=ngenes),
    group=ngenes %>% seq_len %>% rep(4)
  )+geom_point()+scale_colour_manual(values=c('green','red','blue'),name='batch')+
    geom_line(colour='black',alpha=.2)+
    xlab('')+ylab('weight in eigengene' %>% paste(pc))+
    annotate(geom='text',x=-2/3,y=0,colour='white',label=angle(pca1$rotation[,pc],pca3$rotation[,pc]) %>% TeX(output='character'),parse=TRUE)+
    annotate(geom='text',x=0,y=0,colour='white',label=angle(pca1$rotation[,pc],pca2$rotation[,pc]) %>% TeX(output='character'),parse=TRUE)+
    annotate(geom='text',x=2/3,y=0,colour='white',label=angle(pca2$rotation[,pc],pca3$rotation[,pc]) %>% TeX(output='character'),parse=TRUE)
}

batchgenemap<-function(data,batch,pcs=1:2,xshift=NULL,yshift=NULL){
  ngenes<-nrow(data)
  batch%<>%factor
  batch %>% levels -> batches
  data[,batch==batches[1]] %>% t %>% prcomp -> pca1
  data[,batch==batches[2]] %>% t %>% prcomp -> pca2
  data %>% t %>% prcomp -> pca3
  if(xshift %>% is.null) 
    c(max(abs(pca1$rotation[,pcs[1]]))+max(abs(pca3$rotation[,pcs[1]])),0)*2->xshift
  if(yshift %>% is.null)
    c(0,max(abs(pca2$rotation[,pcs[2]]))+max(abs(pca3$rotation[,pcs[2]])))*2->yshift
  ggplot()+aes(
    x=c(
      pca1$rotation[,pcs[1]],
      pca2$rotation[,pcs[1]],
      pca3$rotation[,pcs[1]]
    )+xshift %>% c(0) %>% rep(each=ngenes),
    y=c(
      pca1$rotation[,pcs[2]],
      pca2$rotation[,pcs[2]],
      pca3$rotation[,pcs[2]]
    )+yshift %>% c(0) %>% rep(each=ngenes),
    colour=c(batches[1],batches[2],'') %>% rep(each=ngenes)
  )+geom_point()+
    scale_colour_manual(values=c('green','red','blue'),name='batch')+
    geom_line(aes(group=c(1:ngenes %>% rep(2),-(1:ngenes)) %>% factor),colour='magenta',alpha=.1)+
    geom_line(aes(group=c(1:ngenes,-(1:ngenes),1:ngenes) %>% factor),colour='orange',alpha=.1)+
    geom_line(aes(group=c(-(1:ngenes),1:ngenes %>% rep(2)) %>% factor),colour='cyan',alpha=.1)+
    xlab('coordinate on eigengene' %>% paste(pcs[1]))+ylab('coordinate on eigengene' %>% paste(pcs[2]))
}