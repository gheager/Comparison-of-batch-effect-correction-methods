v<-rnorm(1743) %>% divide_by(norm(matrix(.),'F')) * 2

get(load('combat.Rdata')) -> combat;
rightward1<-NULL;rightward2<-NULL;
for(i in 1:50){
  combat[,batch1=='geod']%<>%add(v)
  rightward1%<>%c(batchgeneangles(combat,batch1)$first_angle)
  rightward2%<>%c(batchgeneangles(combat,batch1)$second_angle)
}
get(load('combat.Rdata')) -> combat;
leftward1=NULL;leftward2<-NULL;
for(i in 1:50){
  combat[,batch1=='geod']%<>%subtract(v)
  leftward1%<>%c(batchgeneangles(combat,batch1)$first_angle)
  leftward2%<>%c(batchgeneangles(combat,batch1)$second_angle)
}
get(load('combat.Rdata')) -> combat;
angle01<-batchgeneangles(combat,batch1)$first_angle
angle02<-batchgeneangles(combat,batch1)$second_angle


#(c(leftward1 %>% rev,angle01,rightward1) %>% abs %>% acos/pi) %>% plot
(ggplot()+aes(x=-50:50)+
    geom_point(aes(y=c(leftward1 %>% rev,angle01,rightward1) %>% abs %>% acos/pi),colour='red')+
    geom_point(aes(y=c(leftward2 %>% rev,angle02,rightward2) %>% abs %>% acos/pi),colour='blue')) %>% plot
