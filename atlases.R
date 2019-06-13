library(purrr)
library(magrittr)

#Human Expression Atlas
batch_data<-"Human Expression Atlas/batch data/"
experiments<-list();for(filename in dir(batch_data)){
  experiments[[filename %>% stringr::str_split("-") %>% extract2(1) %>% extract(2:3) %>% paste(collapse="")]] <- get(load(paste0(batch_data,filename)))$rnaseq
}
batch<-NULL;tissue<-NULL; for(i in experiments %>% seq_along){
  batch%<>%c(names(experiments)[[i]] %>% rep(dim(experiments[[i]])[2]))
  tissue%<>%c(experiments[[i]]$organism_part)
}; batch%<>%factor
tissues <- tissue %>% split(batch)
intersections<-NULL; for(i in tissues %>% seq_along){
  for(j in tissues %>% seq_along){
    intersections%<>%c(length(intersect(tissues[[i]],tissues[[j]])))
  }
};intersections%<>%matrix(length(tissues))%<>%set_colnames(names(tissues))%<>%set_rownames(names(tissues))
igraph::graph_from_adjacency_matrix(intersections) %>% plot

experiments %>% remove.batch.effect(list=.,model=~organism_part,method='none')->Human.Expression.Atlas.nocorrection
experiments %>% remove.batch.effect(list=.,model=~organism_part,method='combat')->Human.Expression.Atlas.ComBat
experiments %>% remove.batch.effect(list=.,model=~organism_part,method='ruv')->Human.Expression.Atlas.RUV

Human.Expression.Atlas.nocorrection %>% save(file='Human Expression Atlas/Human Expression Atlas no correction.Rdata')
Human.Expression.Atlas.ComBat %>% save(file='Human Expression Atlas/Human Expression Atlas ComBat.Rdata')
Human.Expression.Atlas.RUV %>% save(file='Human Expression Atlas/Human Expression Atlas RUV.Rdata')

Human.Expression.Atlas.nocorrection %>% eigenangles.summaryexperiment -> angnone
Human.Expression.Atlas.ComBat %>% eigenangles.summaryexperiment -> angcomb
Human.Expression.Atlas.RUV %>% eigenangles.summaryexperiment -> angruv

list(none=angnone,combat=angcomb,ruv=angruv) %>% transpose -> ang; ang %>% save(file='Human Expression Atlas/angles.Rdata')
'Human Expression Atlas/angles.pdf' %>% pdf;for(b in ang){
  plot(
    ggplot()+aes(x=seq_along(b$none))+
      geom_point(aes(y=b$none),colour='red')+
      geom_point(aes(y=b$combat),colour='green')+
      geom_point(aes(y=b$ruv),colour='blue')
  )
};dev.off()

#Mouse Expression Atlas
batch_data<-"Mouse Expression Atlas/batch data/"
experiments<-list();for(filename in dir(batch_data)){
  experiments[[filename %>% stringr::str_split("-") %>% extract2(1) %>% extract(2:3) %>% paste(collapse="")]] <- get(load(paste0(batch_data,filename)))$rnaseq
}
batch<-NULL;tissue<-NULL; for(i in experiments %>% seq_along){
  batch%<>%c(names(experiments)[[i]] %>% rep(dim(experiments[[i]])[2]))
  tissue%<>%c(experiments[[i]]$organism_part)
}; batch%<>%factor
tissues <- tissue %>% split(batch)
intersections<-NULL; for(i in tissues %>% seq_along){
  for(j in tissues %>% seq_along){
    intersections%<>%c(length(intersect(tissues[[i]],tissues[[j]])))
  }
};intersections%<>%matrix(length(tissues))%<>%set_colnames(names(tissues))%<>%set_rownames(names(tissues))
igraph::graph_from_adjacency_matrix(intersections) %>% plot

experiments[c('GEOD44366','GEOD45278','ERAD169')]<-NULL#removing isolated experiments
batch<-NULL;tissue<-NULL; for(i in experiments %>% seq_along){
  batch%<>%c(names(experiments)[[i]] %>% rep(dim(experiments[[i]])[2]))
  tissue%<>%c(experiments[[i]]$organism_part)
}; batch%<>%factor
tissues <- tissue %>% split(batch)
intersections<-NULL; for(i in tissues %>% seq_along){
  for(j in tissues %>% seq_along){
    intersections%<>%c(length(intersect(tissues[[i]],tissues[[j]])))
  }
};intersections%<>%matrix(length(tissues))%<>%set_colnames(names(tissues))%<>%set_rownames(names(tissues))
'Mouse Expression Atlas/body intersections graph.png' %>% png;igraph::graph_from_adjacency_matrix(intersections) %>% plot;dev.off()

experiments %>% remove.batch.effect(list=.,model=~organism_part,method='none')->Mouse.Expression.Atlas.nocorrection
experiments %>% remove.batch.effect(list=.,model=~organism_part,method='combat')->Mouse.Expression.Atlas.ComBat
experiments %>% remove.batch.effect(list=.,model=~organism_part,method='ruv')->Mouse.Expression.Atlas.RUV

Mouse.Expression.Atlas.nocorrection %>% save(file='Mouse Expression Atlas/Mouse Expression Atlas no correction.Rdata')
Mouse.Expression.Atlas.ComBat %>% save(file='Mouse Expression Atlas/Mouse Expression Atlas ComBat.Rdata')
Mouse.Expression.Atlas.RUV %>% save(file='Mouse Expression Atlas/Mouse Expression Atlas RUV.Rdata')

Mouse.Expression.Atlas.nocorrection %>% eigenangles.summaryexperiment -> angnone
Mouse.Expression.Atlas.ComBat %>% eigenangles.summaryexperiment -> angcomb
Mouse.Expression.Atlas.RUV %>% eigenangles.summaryexperiment -> angruv

list(none=angnone,combat=angcomb,ruv=angruv) %>% transpose -> ang; ang %>% save(file='Mouse Expression Atlas/angles.Rdata')
'Mouse Expression Atlas/angles.pdf' %>% pdf;for(b in ang){
  plot(
    ggplot()+aes(x=seq_along(b$none))+
      geom_point(aes(y=b$none),colour='red')+
      geom_point(aes(y=b$combat),colour='green')+
      geom_point(aes(y=b$ruv),colour='blue')
  )
};dev.off()

#Corn Expression Atlas
batch_data<-"Corn Expression Atlas/batch data/"
experiments<-list();for(filename in dir(batch_data)){
  experiments[[filename %>% stringr::str_split("-") %>% extract2(1) %>% extract(2:3) %>% paste(collapse="")]] <- get(load(paste0(batch_data,filename)))$rnaseq
}
batch<-NULL;tissue<-NULL; for(i in experiments %>% seq_along){
  batch%<>%c(names(experiments)[[i]] %>% rep(dim(experiments[[i]])[2]))
  tissue%<>%c(experiments[[i]]$organism_part)
}; batch%<>%factor
tissues <- tissue %>% split(batch)
intersections<-NULL; for(i in tissues %>% seq_along){
  for(j in tissues %>% seq_along){
    intersections%<>%c(length(intersect(tissues[[i]],tissues[[j]])))
  }
};intersections%<>%matrix(length(tissues))%<>%set_colnames(names(tissues))%<>%set_rownames(names(tissues))
'Corn Expression Atlas/intersections graph.png' %>% png
igraph::graph_from_adjacency_matrix(intersections) %>% plot
dev.off()

experiments %>% remove.batch.effect(list=.,model=~organism_part,method='none')->Corn.Expression.Atlas.nocorrection
experiments %>% remove.batch.effect(list=.,model=~organism_part,method='combat')->Corn.Expression.Atlas.ComBat
experiments %>% remove.batch.effect(list=.,model=~organism_part,method='ruv')->Corn.Expression.Atlas.RUV

Corn.Expression.Atlas.nocorrection %>% save(file='Corn Expression Atlas/Corn Expression Atlas no correction.Rdata')
Corn.Expression.Atlas.ComBat %>% save(file='Corn Expression Atlas/Corn Expression Atlas ComBat.Rdata')
Corn.Expression.Atlas.RUV %>% save(file='Corn Expression Atlas/Corn Expression Atlas RUV.Rdata')

Corn.Expression.Atlas.nocorrection %>% eigenangles.summaryexperiment -> angnone
Corn.Expression.Atlas.ComBat %>% eigenangles.summaryexperiment -> angcomb
Corn.Expression.Atlas.RUV %>% eigenangles.summaryexperiment -> angruv

list(none=angnone,combat=angcomb,ruv=angruv) %>% transpose -> ang; ang %>% save(file='Corn Expression Atlas/angles.Rdata')
get(load('Corn Expression Atlas/angles.Rdata'))->ang
'Corn Expression Atlas/angles.pdf' %>% pdf;for(b in ang){
  plot(
    ggplot()+aes(x=seq_along(b$none))+
      geom_point(aes(y=b$none),colour='red')+
      geom_point(aes(y=b$combat),colour='green')+
      geom_point(aes(y=b$ruv),colour='blue')
  )
};dev.off()
