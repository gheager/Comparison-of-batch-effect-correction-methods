#GENESHIFT
library(magrittr)
library(purrr)

vector.shift<-function(v,k){
  if(k==0) return(v)
  if(k>0) return(c(0 %>% rep(k),v)[seq_along(v)])
  if(k<0) return(c(v[-(1:-k)],0 %>% rep(-k)))
}

best.matching<-function(p,q){
  if(length(p$y)!=length(q$y)) stop('p & q densities should have the same length')
  d<-length(p$y)
  step<-p$x[2]-p$x[1]
  inner.products<-NULL
  for(k in -d:d) inner.products%<>%c(sum(p$y*(q$y %>% vector.shift(k)))*step)
  return(list(
    'offset'=step*(-d:d)[which.max(inner.products)],
    'matching.score'=max(inner.products),
    inner.products=inner.products
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
geneshift<-function(data,batch){
  data %>% apply(
    MARGIN=1,
    FUN=as_mapper(~best.matching(.[batch==0] %>% density, .[batch==1] %>% density))
  ) %>% transpose %>% lapply(simplify)
}

#idea : geneshift algorithm shifts the marginal densities (i.e. of the original covariates)
#why not shifting densities of other axes (principal components for example) ?
#also why only shifting ? why not scaling or doing other moments transformations ? (shifting is equivalent to transform the expectation i.e. 1st moment of the distribution, we could also act on variance, skewness, kurtosis...)
#finally should we keep using inner product to compare densities ? Kullback-Leibler divergence (relative entropy) might be a good alternative...