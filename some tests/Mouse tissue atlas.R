library(magrittr)
library(stringr)
library(purrr)
library(ggplot2)
library(ExpressionAtlas)
library(igraph)

batch_data<-"Mouse Expression Atlas/batch data/"

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

intersections<-NULL;for(i in tissues %>% seq_along){
  for(j in tissues %>% seq_along){
    intersections%<>%c(length(intersect(tissues[[i]],tissues[[j]])))
  }
};intersections%<>%matrix(length(tissues))%<>%set_colnames(names(tissues))%<>%set_rownames(names(tissues))
# intersections<-list()
# for(i in 1:(length(tissues)-1)){
#   l<-list()
#   for(j in (i+1):length(tissues)){
#     l[[names(tissues)[[j]]]]<-intersect(tissues[[i]],tissues[[j]])
#   }
#   intersections[[names(tissues)[[i]]]]<-l
# }
graph_from_adjacency_matrix(intersections) %>% plot
#"Mouse Expression Atlas/intersections graph.png" %>% png; graph_from_adjacency_matrix(intersections) %>% plot; dev.off()

experiments[c('GEOD45278','GEOD44366','ERAD169')]<-NULL #removing isolated experiments

genes<-experiments %>% map(rownames)
common.genes<-genes %>% purrr::reduce(intersect)

data<-NULL;batch<-NULL;tissue<-NULL; for(i in experiments %>% seq_along){
  data%<>%cbind(experiments[[i]] %>% assays %$% counts %>% extract(common.genes,))
  batch%<>%c(names(experiments)[[i]] %>% rep(dim(experiments[[i]])[2]))
  tissue%<>%c(experiments[[i]]$organism_part)
}; batch%<>%factor

data %>% gPCA(batch) -> gdata
gdata %>% viz_gpca(guided=FALSE)

data %>% log1p %>% gPCA(batch) -> gdata
gdata %>% viz_gpca(guided=FALSE)

data%<>%log1p

#data %>% log1p %>% batchgeneangles(batch) -> anglesdata #not really meaningful as datasets have very different samples

#filter <- rowSums(data)!=0
filter <- data %>% t %>% data.frame %>% split(batch) %>% map(~colSums(.)!=0) %>% purrr::reduce(`&`)
filtered <- data %>% extract(filter,)
filtered %>% gPCA(batch) -> gfilt
gfilt %>% viz_gpca(guided=FALSE)+geom_line(aes(group=tissue),colour='black')
filtered %>% batchgeneangles(batch) -> anglesfilt

library(pamr)
filtered %>% list(x=.,batchlabels=batch) %>% pamr.batchadjust %$% x -> pam

library(sva)
filtered %>% ComBat(batch,model.matrix(~tissue)) -> combat
combat %>% gPCA(batch) -> gcombat
gcombat %>% viz_gpca(guided=FALSE)#+geom_line(aes(group=tissue),colour='black')

library(RUVSeq)
filtered %>% RUVs(cIdx=rownames(filtered),k=1,scIdx=makeGroups(tissue),isLog=TRUE) -> ruv

filtered %>% t %>% prcomp(rank.=2) -> filtpca
ggplot()+aes(x=filtpca$x[,1],y=filtpca$x[,2],colour=tissue)+geom_point()+stat_ellipse()+
  scale_colour_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')[1:14])#scale_colour_manual(values=c('red','orange','yellow','green','blue','aliceblue','purple','pink','black','grey','brown','magenta','cyan','white'))

pam %>% t %>% prcomp(rank.=2) -> pampca
ggplot()+aes(x=pampca$x[,1],y=pampca$x[,2],colour=tissue)+geom_point()+stat_ellipse()+
  scale_colour_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')[1:14])#scale_colour_manual(values=c('red','orange','yellow','green','blue','aliceblue','purple','pink','black','grey','brown','magenta','cyan','white'))

combat %>% t %>% prcomp -> combatpca
ggplot()+aes(x=combatpca$x[,1],y=combatpca$x[,2],colour=tissue)+geom_point()+stat_ellipse()+
  scale_colour_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')[1:14])#scale_colour_manual(values=c('red','orange','yellow','green','blue','aliceblue','purple','pink','black','grey','brown','magenta','cyan','white'))
"Mouse Expression Atlas/ComBat PCA.png" %>% png; ggplot()+aes(x=combatpca$x[,1],y=combatpca$x[,2],colour=tissue)+geom_point()+stat_ellipse()+
  scale_colour_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')[1:14]); dev.off()

ruv$normalizedCounts %>% t %>% prcomp(rank.=2) -> ruvpca
ggplot()+aes(x=ruvpca$x[,1],y=ruvpca$x[,2],colour=tissue)+geom_point()+stat_ellipse()+
  scale_colour_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')[1:14])#scale_colour_manual(values=c('red','orange','yellow','green','blue','aliceblue','purple','pink','black','grey','brown','magenta','cyan','white'))

combat %>% gPCA(tissue,scaleY=TRUE) -> gtcombat
combat %>% gPCA(tissue) -> gtcombat
#gtcombat %>% viz_gpca(guided=FALSE)
gtcombat %>% viz_gpca(dims=c(1,3)) + scale_colour_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')[1:14])#scale_colour_manual(values=c('red','orange','yellow','green','blue','aliceblue','purple','pink','black','grey','brown','magenta','cyan','white'))

filtered %>% eigenangles(batch,tissue) -> angfilt
combat %>% eigenangles(batch,tissue) -> angcomb
pam %>% eigenangles(batch,tissue) -> angpam
ruv$normalizedCounts %>% eigenangles(batch,tissue) -> angruv

angfilt %>% save(file="Mouse Expression Atlas/angfilt.Rdata")
angcomb %>% save(file="Mouse Expression Atlas/angcomb.Rdata")
angpam %>% save(file="Mouse Expression Atlas/angpam.Rdata")
angruv %>% save(file="Mouse Expression Atlas/angruv.Rdata")

list(filtered=angfilt,combat=angcomb,pam=angpam,ruv=angruv) %>% transpose -> ang
'Mouse Expression Atlas/ang.pdf' %>% pdf;for(b in ang){
  plot(
    ggplot()+aes(x=seq_along(b$filtered))+
      geom_point(aes(y=b$filtered),colour='red')+
      geom_point(aes(y=b$combat),colour='green')+
      geom_point(aes(y=b$pam),colour='blue')+
      geom_point(aes(y=b$ruv),colour='orange')+scale_y_log10()
  )
};dev.off()

combat %>% save(file="Mouse tissue atlas/Mouse Expression Atlas across tissues (ComBat).Rdata")

# combat %>% batchgeneangles(batch) -> anglescombat
# combat %>% batchgeneangles2(batch,tissue,global.only = TRUE) -> angles2combat
# 
# filtered %>% batchgeneangles2(batch,tissue,global.only = TRUE)->angles2filt
# pam %>% batchgeneangles2(batch,tissue,global.only = TRUE)->angles2pam

library(gg3D)
ggplot()+aes(x=combatpca$x[,1],y=combatpca$x[,2],z=combatpca$x[,3],colour=tissue)+theme_void()+axes_3D()+stat_3D()+
  scale_colour_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')[1:14])#scale_colour_manual(values=c('red','orange','yellow','green','blue','aliceblue','purple','pink','black','grey','brown','magenta','cyan','white'))
