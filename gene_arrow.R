library(purrr)
library(magrittr)

gene.arrow.2batches<-function(data,batch,model,weights=rep(1/2,2)){
  batch%<>%factor; if(nlevels(batch)!=2) stop('There should be only 2 batches.') else levels(batch)<-0:1
  intersect(model[batch==0],model[batch==1])->intersection
  arrows<-NULL
  for(m in intersection){
    arrows%<>%cbind(rowMeans(data[,model==m & batch==1])-rowMeans(data[,model==m & batch==0]))
  }
  data[,batch==1]%<>%subtract(rowMeans(arrows)*weights[2])
  data[,batch==0]%<>%add(rowMeans(arrows)*weights[1])
  return(data)
}

gene.arrow<-function(data,batch,model){
  batch%<>%factor
  if(nlevels(batch)==2){
    return(gene.arrow.2batches(data,batch,model))
  }else{
    batch %>% levels -> batches
    return(gene.arrow.2batches(
      gene.arrow.2batches()
    ))
  }
}