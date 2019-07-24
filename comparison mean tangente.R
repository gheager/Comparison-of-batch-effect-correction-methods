ruv %>% save(file='RUV human k=1:10.Rdata')
mnn %>% save(file='MNN human k=1:10.Rdata')

ruv<-get(load('RUV k=1:10.Rdata'))
mnn<-get(load('MNN k=1:10.Rdata'))

gruv<-ruv %>% map(do_gPCA)
gruv %>% map(~.$delta) %>% map(~.[1]) %>% unlist %>% plot

gmnn<-mnn %>% map(do_gPCA)
gmnn %>% map(~.$delta) %>% map(~.[1]) %>% unlist %>% plot

gcombat<-do_gPCA(combat)

ggplot(data.frame(
  k=1:5 %>% rep(2),
  method=c('MNN','RUV') %>% rep(each=5),
  delta1=c(gruv %>% map(~.$delta) %>% map(~.[1]) %>% unlist, gmnn %>% map(~.$delta) %>% map(~.[1]) %>% unlist)
))+geom_point(aes(x=k,y=delta1,colour=method))+
  geom_hline(yintercept = gcombat$delta[1],colour='red')+annotate('text',label='ComBat',x=2,y=gcombat$delta[1]+.01,colour='red')+
  geom_hline(yintercept = gall$delta[1],colour='black')+annotate('text',label='no correction',x=2,y=gall$delta[1]+.01,colour='black')+
  scale_colour_manual(values=c('green','blue'))+ylim(c(0,1))+
  ggtitle('Delta statistic for Homo sapiens')

ggplot(data.frame(
  k=1:10 %>% rep(2),
  method=c('MNN','RUV') %>% rep(each=10),
  rank1=c(gruv %>% map(~.$ranks) %>% map(~.[1]) %>% unlist, gmnn %>% map(~.$ranks) %>% map(~.[1]) %>% unlist)
))+geom_point(aes(x=k,y=rank1,colour=method))+
  geom_hline(yintercept = gcombat$ranks[1],colour='red')+
  geom_hline(yintercept = gall$ranks[1],colour='black')+
  scale_colour_manual(values=c('green','blue'))+
  ggtitle('Rank for Mus musculus')

library(eigenangles)

all<-integrate.experiments(list=experiments, method='none', model=~organism_part)

combat<-integrate.experiments(list=experiments, method='combat', model=~organism_part)

batch.modification<-function(z){
  all$batch %>% levels %>% map(
    function(batch) eigenangles:::angledet(
      all %>% assays %>% use_series(corrected) %>% extract(,all$batch==batch) %>% t %>% prcomp(rank=1) %$% rotation,
      z %>% assays %>% use_series(corrected) %>% extract(,z$batch==batch) %>% t %>% prcomp(rank=1) %$% rotation
    )
  ) %>% unlist %>% tanpi %>% mean %>% atan/pi
}

batch.modification<-function(z){
  batch %>% levels %>% map(
    function(b) eigenangles:::angledet(
      all %>%  extract(,batch==b) %>% t %>% prcomp(rank=1) %$% rotation,
      z %>% extract(,batch==b) %>% t %>% prcomp(rank=1) %$% rotation
    )
  ) %>% unlist %>% tanpi %>% mean %>% atan/pi
}

batch.modification(combat) -> combat.modif
ruv %>% map(batch.modification) %>% unlist -> ruv.modif
mnn %>% map(batch.modification) %>% unlist -> mnn.modif

ruv.modif %>% plot

ggplot(data.frame(
  k=1:5 %>% rep(2),
  method=c('MNN','RUV') %>% rep(each=5),
  modif=c(mnn.modif,ruv.modif)
))+geom_point(aes(x=k,y=modif,colour=method))+
  geom_hline(yintercept = combat.modif, colour='red')+annotate('text',label='ComBat',x=2,y=combat.modif,colour='red')+
  geom_hline(yintercept=0,colour='black')+annotate('text',label='PC1 of original batch',x=5,y=0,colour='black')+
  xlim(c(0,5))+ylim(c(0,2))+coord_polar(theta='y')+ggtitle('Angles of batch modification for Homo sapiens')+scale_colour_manual(values=c('green','blue'))#+
  #scale_x_log10()

ggplot(data.frame(
  k=1:10 %>% rep(2),
  method=c('MNN','RUV') %>% rep(each=10),
  meantan=c(mnn.meantan,ruv.meantan)
))+geom_point(aes(x=k,y=meantan,colour=method))+
  geom_hline(yintercept = combat.meantan, colour='red')+annotate('text',label='ComBat',x=4,y=combat.meantan,colour='red')+
  geom_hline(yintercept=all.meantan,colour='brown')+annotate('text',label='no correction',x=5,y=all.meantan,colour='brown')+
  geom_hline(yintercept=0,colour='black')+annotate('text',label='complete dataset',x=10,y=0,colour='black')+
  xlim(c(0,10))+ylim(c(0,2))+coord_polar(theta='y')+scale_colour_manual(values=c('green','blue'))+ggtitle('Angles of batch integration for Mus musculus')#+
  #scale_x_log10()

ggplot(data.frame(
  k=1:20 %>% rep(2),
  method=c('MNN','RUV') %>% rep(each=20),
  meantan=c(mnn.meantan+mnn.modif,ruv.meantan+mnn.modif)
))+geom_point(aes(x=k,y=meantan,colour=method))+
  geom_hline(yintercept = combat.meantan+combat.modif, colour='red')+annotate('text',label='ComBat',x=5,y=combat.meantan,colour='red')+
  geom_hline(yintercept=all.meantan,colour='black')+annotate('text',label='PC1 of original batch',x=20,y=0,colour='black')+
  ylim(c(0,2))+coord_polar(theta='y')+scale_x_log10()+ggtitle('Angles of batch integration')+scale_colour_manual(values=c('green','blue'))


# ggplot(data.frame(
#   k=1:10 %>% rep(2),
#   method=c('MNN','RUV') %>% rep(each=10),
#   meantan=c(mnn.meantan,ruv.meantan)
# ))+geom_segment(aes(x=k,y=meantan,colour=method),xend=0,yend=0)+
#   geom_hline(yintercept = combat.meantan, colour='brown')+annotate('text',label='ComBat',x=5,y=combat.meantan,colour='brown')+
#   ylim(c(0,2))+coord_radar(theta='y')+scale_x_log10()

eigenangles:::angledet(
  all %>% assays %>% use_series(corrected) %>% extract(,all$batch=='GEOD43721') %>% t %>% prcomp(rank=1) %$% rotation,
  combat %>% assays %>% use_series(corrected) %>% extract(,combat$batch=='GEOD43721') %>% t %>% prcomp(rank=1) %$% rotation
)

ruv<-list()
for(k in 1:5){ruv[[k]]<-integrate.experiments(list=experiments, method="ruv", model = ~organism_part, k=k);cat(k)}
ruvang<-ruv %>% map(do_eigenangles %>% partial(group='organism_part'))

ruvang %>% map(~.$batch_vs_all) %>% map(function(x) x %>% map(~.[1])) %>% map(unlist) %>% map(tanpi) %>% map(mean) %>% map(atan) %>% map(~./pi) %>% unlist->ruv.meantan
ruvang %>% map(~.$batch_vs_all) %>% map(function(x) x %>% map(~.[1])) %>% map(unlist) %>% map(mean) %>% unlist->ruv.mean

mnn<-list()
for(k in 11:20){mnn[[k]]<-integrate.experiments(list=experiments, method="mnn", model = ~organism_part, k=k);cat(k)}
mnnang<-mnn %>% map(do_eigenangles %>% partial(group='organism_part'))

mnnang %>% map(~.$batch_vs_all) %>% map(function(x) x %>% map(~.[1])) %>% map(unlist) %>% map(tanpi) %>% map(mean) %>% map(atan) %>% map(~./pi) %>% unlist->mnn.meantan
mnnang %>% map(~.$batch_vs_all) %>% map(function(x) x %>% map(~.[1])) %>% map(unlist) %>% map(mean) %>% unlist->mnn.mean

angc<-combat %>% do_eigenangles('organism_part')
angc$batch_vs_all %>% map(~.[1]) %>% unlist %>% tanpi %>% mean %>% atan/pi -> combat.meantan
ang$batch_vs_all %>% map(~.[1]) %>% unlist %>% tanpi %>% mean %>% atan/pi -> all.meantan
angc$batch_vs_all %>% map(~.[1]) %>% unlist %>% mean->combat.mean

all %>% do_eigenangles('organism_part')->ang
ang$batch_vs_all %>% map(~.[1]) %>% unlist %>% tanpi %>% mean %>% atan/pi -> none.meantan
ang$batch_vs_all %>% map(~.[1]) %>% unlist %>% mean->none.mean

ggplot()+aes(x=1:10)+
  geom_point(aes(y=mnn.meantan),colour='green')+
  geom_point(aes(y=ruv.meantan),colour='blue')+
  geom_hline(aes(yintercept=combat.meantan),colour='red')+
  geom_hline(aes(yintercept=none.meantan),colour='black')+
  ylab('meantan')+scale_x_discrete(name='k',breaks=1:10,labels=1:10,limits=1:10)

ggplot()+aes(x=1:10)+
  geom_point(aes(y=mnn.mean),colour='green')+
  geom_point(aes(y=ruv.mean),colour='blue')+
  geom_hline(aes(yintercept=combat.mean),colour='red')+
  geom_hline(aes(yintercept=none.mean),colour='black')+
  ylab('mean')+scale_x_discrete(name='k',breaks=1:10,labels=1:10,limits=1:10)

