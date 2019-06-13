#GENESHIFT
library(magrittr)
library(purrr)
library(ggplot2)
library(ExpressionAtlas)

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

geneshift<-function(data,batch,step){
  batch%<>%factor
  batches<-levels(batch)
  data %>% apply(
    MARGIN=1,
    FUN=as_mapper(~best.matching(.[batch==batches[1]] %>% density, .[batch==batches[2]] %>% density,step))
  ) %>% transpose %>% lapply(simplify2array)
}

geneshift.pca<-function(data,batch,npc,step){
  batch%<>%factor
  batches<-levels(batch)
  data %<>% t %<>% scale(scale=FALSE)
  data %>% svd -> pca
  pca$u[,1:npc] %>% t %>% geneshift(batch,step) -> correction
  pca$u[batch==batches[2],1:npc] %<>% t %<>% add(correction$offset) %<>% t
  return(t(
    pca$u%*%diag(pca$d)%*%t(pca$v)
  ) %>% structure(
    geneshift=correction,
    new.u=pca$u[,1:npc],
    batch=batch %>% factor
  ))
}
