viz_angles<-function(...,filename){
  list(...) %>% transpose -> angles
  filename %>% pdf
  for(batch in angles){
    plot(ggplot()+xlim(c(0,1.25))+ylim(c(0,1))+
      aes(
        x=batch %>% unlist %>% cospi,
        y=batch %>% unlist %>% sinpi,
        colour=batch %>% imap(~rep(..2,each=length(..1))) %>% unlist %>% factor
      )+
      geom_segment(xend=0,yend=0)+
      geom_text(label=batch %>% unlist %>% round(3) %>% paste('$\\pi$') %>% TeX,parse=TRUE,nudge_x=.1)+
      geom_arc(aes(x0=0,y0=0,r=1,start=0,end=pi/2),colour='black',inherit.aes = FALSE)+coord_fixed())
  };dev.off()
}

ang<-get(load('Human Expression Atlas/angles.Rdata'))
ang%<>%transpose
angnone<-ang$none
angcomb<-ang$comb
angruv<-ang$ruv
viz_angles(none=angnone,combat=angcomb,ruv=angruv,filename='Human Expression Atlas/angles viz.pdf')


# viz_angles<-function(...){
#   list(...) %>% transpose %>% map(
#     function(a){
#       ggplot()+xlim(c(0,1.25))+ylim(c(0,1))+
#         aes(
#           x=a %>% unlist %>% cospi,
#           y=a %>% unlist %>% sinpi,
#           colour=a %>% imap(~rep(..2,each=length(..1))) %>% unlist %>% factor
#         )+
#         geom_segment(xend=0,yend=0)+
#         geom_text(label=a %>% unlist %>% round(3) %>% paste('$\\pi$') %>% TeX,parse=TRUE,nudge_x=.1)+
#         geom_arc(aes(x0=0,y0=0,r=1,start=0,end=pi/2),colour='black',inherit.aes = FALSE)+coord_fixed()
#     }
#   ) %>% arrangeGrob(grobs=.)
# }

# viz_angles<-function(...){
#   list(...) %>% transpose -> angles
#   ggplot()+xlim(c(0,1.25))+ylim(c(0,1))+
#     aes(
#       x=angles %>% unlist %>% cospi,
#       y=angles %>% unlist %>% sinpi,
#       colour=angles %>% map(imap %>% partial(...=,~rep(..2,each=length(..1)))) %>% unlist %>% factor
#     )+
#     geom_segment(xend=0,yend=0)+
#     geom_text(label=angles %>% unlist %>% round(3) %>% paste('$\\pi$') %>% TeX,parse=TRUE,nudge_x=.1)+
#     geom_arc(aes(x0=0,y0=0,r=1,start=0,end=pi/2),colour='black',inherit.aes = FALSE)+coord_fixed()+
#     facet_grid(angles %>% imap(~rep(..2,each=length(unlist(..1)))) %>% unlist %>% vars)
# }
