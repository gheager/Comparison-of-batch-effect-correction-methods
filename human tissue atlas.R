library(magrittr)
library(stringr)
library(purrr)
library(ggplot2)
library(ExpressionAtlas)
library(igraph)

batch_data<-"Human tissue atlas/batch data/"

experiments<-list();for(filename in dir(batch_data)){
  experiments[[filename %>% str_split("-") %>% extract2(1) %>% extract(2:3) %>% paste(collapse="")]] <- get(load(paste0(batch_data,filename)))$rnaseq
}

genes<-experiments %>% map(rownames)
common.genes<-genes %>% purrr::reduce(intersect)

data<-NULL;batch<-NULL;tissue<-NULL; for(i in experiments %>% seq_along){
  data%<>%cbind(experiments[[i]] %>% assays %$% counts %>% extract(common.genes,))
  batch%<>%c(names(experiments)[[i]] %>% rep(dim(experiments[[i]])[2]))
  tissue%<>%c(experiments[[i]]$organism_part)
}; batch%<>%factor

tissues <- tissue %>% split(batch)
intersections<-NULL
for(i in tissues %>% seq_along){
  for(j in tissues %>% seq_along){
    intersections%<>%c(length(intersect(tissues[[i]],tissues[[j]])))
  }
};intersections%<>%matrix(length(tissues))%<>%set_colnames(names(tissues))%<>%set_rownames(names(tissues))
graph_from_adjacency_matrix(intersections) %>% plot

#no isolated experiment

# data %>% gPCA(batch) -> gdata
# gdata %>% viz_gpca(guided=FALSE)
# 
# data %>% log1p %>% gPCA(batch) -> gdata
# gdata %>% viz_gpca(guided=FALSE)

data%<>%log1p

#filter<-TRUE %>% rep(nrow(data))
#filter <- data %>% t %>% data.frame %>% split(tissue) %>% map(~colSums(.)!=0) %>% purrr::reduce(`&`)
filter <- data %>% t %>% data.frame %>% split(batch) %>% map(~colSums(.)!=0) %>% purrr::reduce(`&`)
filtered<-data[filter,]
#filtered[filtered==0]<-1

library(pamr)
filtered %>% list(x=.,batchlabels=batch) %>% pamr.batchadjust %$% x -> pam
eigenangles(pam,batch,tissue)->angpam

pam %>% t %>% prcomp(rank.=2) -> pampca
ggplot()+aes(x=pampca$x[,1],y=pampca$x[,2],colour=tissue)+geom_point()+stat_ellipse()+
  scale_colour_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000') %>% rep(3))


library(sva)
filtered %>% ComBat(batch,model.matrix(~tissue)) -> combat
combat %>% t %>% prcomp(rank.=2) -> pcacombat

ggplot()+aes(x=pcacombat$x[,1],y=pcacombat$x[,2],colour=tissue)+geom_point()+stat_ellipse()+
  scale_colour_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000') %>% rep(3))

eigenangles(combat,batch,tissue)->angcomb
