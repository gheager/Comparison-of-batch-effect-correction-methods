### mean intra-class inertia function written using dplyr grammar
### (very slow as dplyr's group_by function creates quite an abstract structure, which is not scalable with the size of our data)
library(dplyr)
batch_correction_inertia <- function(data, outcome, batch){
  data %>% t %>% data.frame %>% cbind(outcome, batch) %>% 
    group_by(outcome) %>% summarise_at(vars(-batch),var) %>% extract(,-1) %>% na.omit %>% 
    rowMeans %>% sqrt %>% mean
}

### same function written from scratch (much faster)
library(magrittr)
library(matrixStats)
batch_inertia <- function(data, outcome){
  outcome%<>%factor
  vars<-NULL
  for(x in outcome %>% levels){
    if(sum(outcome==x)>1){
      vars%<>%c(rowVars(data[,which(outcome==x)],na.rm=TRUE) %>% mean %>% sqrt)
    }
  }
  return(vars %>% mean)
}

inertias <- methods %>% sapply(
  function(method) method %>% get %>% batch_inertia(tissue)
); names(inertias)<-methods
batch_inertia(filtered, tissue)
batch_inertia(pam, tissue)
batch_inertia(combat, tissue)
batch_inertia(ccombat,tissue)
batch_inertia(rg, tissue)
batch_inertia(ra, tissue)
batch_inertia(ruvg, tissue)
#to be compared with such graphics:
ggplot()+aes(x=pca_combat$ind$coord[,1],y=pca_combat$ind$coord[,2],colour=tissue)+
  geom_point()+stat_ellipse(level=0.95)+theme(legend.position='none')
ggplot()+aes(x=pca_ra$ind$coord[,1],y=pca_ra$ind$coord[,2],colour=tissue)+
  geom_point()+stat_ellipse(level=0.95)+theme(legend.position='none')

