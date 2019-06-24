colours20=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
library(SummarizedExperiment)
#Human Expression Atlas
get(load('Corn Expression Atlas/integrated data/Corn Expression Atlas RUV.Rdata')) -> Human.Expression.Atlas.ComBat

Human.Expression.Atlas.ComBat %>% assays %$% counts %>% t %>% prcomp(rank.=2)->humanpca
ggplot()+aes(
  x=humanpca$x[,1],
  y=humanpca$x[,2],
  colour=Human.Expression.Atlas.ComBat$organism_part,
  label=Human.Expression.Atlas.ComBat$organism_part
)+geom_text(size=3) + scale_colour_manual(values=colours20 %>% rep(4))+#stat_ellipse()+
  theme(legend.position = 'none')+ggtitle('PCA Corn Expression Atlas, RUVs')

Human.Expression.Atlas.ComBat %>% gPCA.integrate('organism_part') -> ghuman
ghuman %>% viz_gpca + scale_colour_manual(values=colours20 %>% rep(3))+
  theme(legend.position = 'none')+ggtitle('PCA Human Expression Atlas, no correction')

#Mouse Expression Atlas
get(load('Mouse Expression Atlas/integrated data/Mouse Expression Atlas RUV.Rdata')) -> Human.Expression.Atlas.ComBat
get(load('Mouse Expression Atlas/Mouse Expression Atlas RUV.Rdata')) -> Mouse.Expression.Atlas.RUV
get(load('Mouse Expression Atlas/Mouse Expression Atlas no correction.Rdata')) -> Mouse.Expression.Atlas.none
Mouse.Expression.Atlas.ComBat %>% eigenangles.summaryexperiment -> angcomb

Mouse.Expression.Atlas.ComBat %>% gPCA.integrate('organism_part') -> gmouse
gmouse %>% viz_gpca + scale_colour_manual(values=colours20)

#Corn Expression Atlas
get(load('Corn Expression Atlas/integrated data/Corn Expression Atlas RUV.Rdata')) -> Human.Expression.Atlas.ComBat

Corn.Expression.Atlas.ComBat %>% assays %$% counts %>% t %>% prcomp(rank.=2)->Cornpca
ggplot()+aes(
  x=Cornpca$x[,1],
  y=Cornpca$x[,2],
  colour=Corn.Expression.Atlas.ComBat$organism_part,
  label=Corn.Expression.Atlas.ComBat$organism_part
)+geom_text(size=3)+stat_ellipse() + scale_colour_manual(values=colours20 %>% rep(4))+theme(legend.position = 'none')

Corn.Expression.Atlas.ComBat %>% gPCA.integrate('organism_part') -> gcorn
gcorn %>% viz_gpca + scale_colour_manual(values=colours20 %>% rep(4))
