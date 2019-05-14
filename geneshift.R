#GENESHIFT
library(magrittr)
library(purrr)
library(ggplot2)
library(ExpressionAtlas)

vector.shift<-function(v,k){
  if(k==0) return(v)
  if(k>0) return(c(0 %>% rep(k),v)[seq_along(v)])
  if(k<0) return(c(v[-(1:-k)],0 %>% rep(-k)))
}

# best.matching<-function(p,q){
#   if(length(p$y)!=length(q$y)) stop('p & q densities should have the same length')
#   d<-length(p$y)
#   step<-p$x[2]-p$x[1]
#   inner.products<-NULL
#   for(k in -d:d) inner.products%<>%c(sum(p$y*(q$y %>% vector.shift(k)))*step)
#   return(list(
#     'offset'=step*(-d:d)[which.max(inner.products)],
#     'matching.score'=max(inner.products),
#     inner.products=inner.products,
#     step=step,
#     x=step*(-d:d)
#   ))
# }

best.matching<-function(p,q,step){
  rightward<-NULL
  qq<-q
  repeat{
    qq$x %<>% add(step); if(min(qq$x)>max(p$x)) break
    rightward %<>% c(
      integrate(
        as_mapper(~approxfun(p)(.)*approxfun(qq)(.)),
        lower=max(min(p$x),min(qq$x)),
        upper=min(max(p$x),max(qq$x))
      )$value
    )
  }
  leftward<-NULL
  qq<-q
  repeat{
    qq$x %<>% subtract(step); if(max(qq$x)<min(p$x)) break
    leftward %<>% c(
      integrate(
        as_mapper(~approxfun(p)(.)*approxfun(qq)(.)),
        lower=max(min(p$x),min(qq$x)),
        upper=min(max(p$x),max(qq$x))
      )$value
    )
  }
  xoffset<-c(
    -step*rev(seq_along(leftward)),
    0,
    step*seq_along(rightward)
  )
  scores<-c(
    leftward %>% rev,
    integrate(as_mapper(~approxfun(p)(.)*approxfun(q)(.)),
              lower=max(min(p$x),min(q$x)),upper=min(max(p$x),max(q$x)))$value,
    rightward
  )
  return(list(
    offset=xoffset[which.max(scores)],
    matching.score=max(scores),
    xoffset=xoffset,
    scores=scores
  ))
}

b<-best.matching(A[1,batch==0] %>% density,A[1,batch==1] %>% density)

#with gradient ascent
#/!\ the inner product is not a concave function of the shift parameter...
#when the distributions are multimodal, this function could even have several local maxima...
inner.product_shift.derivative<-function(x0,x1){
  
}
best.matching.shift<-function(x0,x1){
  derivative<-o
}

#to be parallelized...
#or think to a gradient ascent on the shift parameter (maybe an exact resolution exists for some kernel choice but it would be astonishing...)
geneshift<-function(data,batch,step){
  data %>% apply(
    MARGIN=1,
    FUN=as_mapper(~best.matching(.[batch==0] %>% density, .[batch==1] %>% density,step))
  ) %>% transpose %>% lapply(simplify2array)
}

#idea : geneshift algorithm shifts the marginal densities (i.e. of the original covariates)
#why not shifting densities of other axes (principal components for example) ?
#also why only shifting ? why not scaling or doing other moments transformations ? (shifting is equivalent to transform the expectation i.e. 1st moment of the distribution, we could also act on variance, skewness, kurtosis...)
#finally should we keep using inner product to compare densities ? Kullback-Leibler divergence (relative entropy) might be a good alternative...

pca<-B %>% gPCA(batch) %>% use_series(pca)
pca$v[1:10,1:5]
C<-NULL
for(i in 1:5) C%<>%cbind(rnorm(6000,mean=.1*batch[i]))
projCB<-t(t(C)%*%pca$v)
CC<-t(t(projCB)%*%t(pca$v))

ggplot()+geom_density(aes(x=pca$u[,1],colour=batch %>% factor))
ggplot()+geom_density(aes(x=pca$u[,2],colour=batch %>% factor))
ggplot()+geom_density(aes(x=pca$u[,3],colour=batch %>% factor))

ggplot()+geom_density(aes(x=B[1,],colour=batch %>% factor))
ggplot()+geom_density(aes(x=B[2,],colour=batch %>% factor))
ggplot()+geom_density(aes(x=B[3,],colour=batch %>% factor))


#geneshift performed on principal components with partial reconstruction
geneshift.pca<-function(data,batch,npc,step){
  data %<>% t %<>% scale(scale=FALSE)
  data %>% svd -> pca
  pca$u[,1:npc] %>% t %>% geneshift(batch,step) -> correction
  #ppca<-pca
  #ppca$u[batch==1,1:npc] %<>% t %<>% add(correction$offset) %<>% t
  pca$u[batch==1,1:npc] %<>% t %<>% add(correction$offset) %<>% t
  #avec la ligne precedente decommentee
  return(t(
    pca$u%*%diag(pca$d)%*%t(pca$v)
  ) %>% structure(
    geneshift=correction,
    new.u=pca$u[,1:npc],
    batch=batch %>% factor
  ))
}

viz_geneshift<-function(gs,pc){
  ggplot(mapping=aes_string(
    x=gs %>% attr('new.u') %>% extract(,pc),
    colour=gs %>% attr('batch')
  ))+geom_density()
}

batch.density<-function(x,batch,offset){
  x[batch==1]%<>%add(offset)
  # p<-x[batch==0] %>% density
  # q<-x[batch==1] %>% density
  ggplot(mapping=aes_string(
    x=x, colour=batch %>% factor
  ))+geom_density()#+geom_text(aes(x=0,y=1,label=matching.score(p,q)))
}

matching.score<-function(p,q){
  step<-p$x[2]-p$x[1]
  sum(step*p$y*q$y)
}

gb<-B %>% gPCA(batch)
gb %>% viz_gpca_contrib

batch%<>%factor

'geneshift 3PCs.pdf' %>% pdf
g<-geneshift.pca(B,batch,20,.01)
gg<-g %>% gPCA(batch)
gb %>% viz_gpca_contrib %>% print
gg %>% viz_gpca_contrib %>% print
#gg %>% viz_gpca %>% add(ggtitle('corrected data')) %>% print
#effect on the PCs of corrected data
#before
ggplot(mapping=aes(
  x=t(B)%*%gg$pca$v[,1],
  y=t(B)%*%gg$pca$v[,2],
  colour=batch
))+geom_point()+stat_ellipse() + ggtitle('PCs of corrected data before correction') %>% print
#after
gg %>% viz_gpca(guided=FALSE) + ggtitle('PCs of corrected data after correction') %>% print

#effect on the PCs of raw data
#before
gb %>% viz_gpca(guided=FALSE) + ggtitle('PCs of raw data before correction') %>% print
#after
ggplot(mapping=aes(
  x=t(g)%*%gb$pca$v[,1],
  y=t(g)%*%gb$pca$v[,2],
  colour=batch
))+geom_point()+stat_ellipse() + ggtitle('PCs of raw data after correction') %>% print
dev.off()

#effect on the gPCs of raw data
#before
gb %>% viz_gpca + ggtitle('gPCs of raw data before correction') %>% print
#after
ggplot(mapping=aes(
  x=t(g)%*%gb$gpca$v[,1],
  y=t(g)%*%gb$gpca$v[,2],
  colour=batch
))+geom_point()+stat_ellipse() + ggtitle('PCs of raw data after correction') %>% print

#effect on the gPCs of corrected data
#before
gg %>% viz_gpca + ggtitle('gPCs of raw data before correction') %>% print
#after
ggplot(mapping=aes(
  x=t(B)%*%gg$gpca$v[,1],
  y=t(B)%*%gg$gpca$v[,2],
  colour=batch
))+geom_point()+stat_ellipse() + ggtitle('gPCs of corrected data after correction') %>% print



library(patchwork)
wrap_plots(
  ggplot()+geom_density(aes(x=t(g)%*%gb$pca$v[,1],colour=batch)),
  plot_spacer(),
  ggplot(mapping=aes(
    x=t(g)%*%gb$pca$v[,1],
    y=t(g)%*%gb$pca$v[,2],
    colour=batch
  ))+geom_point()+stat_ellipse(),
  ggplot()+geom_density(aes(x=t(g)%*%gb$pca$v[,2],colour=batch))+coord_flip()
)
