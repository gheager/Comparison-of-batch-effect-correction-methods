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

RUV.Atlases<-list();for(k in 1:10){
  cat(k)
  RUV.Atlases[[k]]<-experiments %>% remove.batch.effect(list=.,model=~organism_part,method='ruv',k=k)
}
RUV.Atlases %>% save(file='Mouse Expression Atlas/RUV.Atlases.Rdata')

RUV.Atlases %>% map(eigenangles.summaryexperiment) -> angruv
angruv %>% save(file='Mouse Expression Atlas/angruv_k=1:10.Rdata')

'Mouse Expression Atlas/angles RUV k=1:10.pdf' %>% pdf;for(b in angruv %>% transpose){
  plot(
    ggplot()+aes(x=seq_along(b$none))+
      geom_point(aes(y=b$none),colour='red')+
      geom_point(aes(y=b$combat),colour='green')+
      geom_point(aes(y=b$ruv),colour='blue')
  )
};dev.off()

get(load('Mouse Expression Atlas/angruv_k=1:10.Rdata'))->angruv
#angruv %>% transpose %>% map(~data.frame(matrix(unlist(.),ncol=length(.[[1]])))) -> angruv
'Mouse Expression Atlas/angles RUV k=1:10.pdf' %>% pdf;for(b in angruv %>% transpose){
  plot(
    ggplot()+aes(
      x=b %>% map(seq_along) %>% unlist,
      y=b %>% unlist,
      colour=1:10 %>% rep(each=b[[1]] %>% length)
    )+geom_point()+scale_colour_gradient(low='blue',high='red')
  )
};dev.off()
