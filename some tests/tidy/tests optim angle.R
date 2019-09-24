v<-rnorm(1743) %>% divide_by(norm(matrix(.),'F')) * 2

pamr.batchadjust(list(x=whole1,batchlabels=batch1))$x->pam;
rightward1<-NULL;rightward2<-NULL;
for(i in 1:50){
  pam[,batch1=='geod']%<>%add(v)
  rightward1%<>%c(batchgeneangles(pam,batch1)$first_angle)
  rightward2%<>%c(batchgeneangles(pam,batch1)$second_angle)
}
pamr.batchadjust(list(x=whole1,batchlabels=batch1))$x->pam;
leftward1=NULL;leftward2<-NULL;
for(i in 1:50){
  pam[,batch1=='geod']%<>%subtract(v)
  leftward1%<>%c(batchgeneangles(pam,batch1)$first_angle)
  leftward2%<>%c(batchgeneangles(pam,batch1)$second_angle)
}
pamr.batchadjust(list(x=whole1,batchlabels=batch1))$x->pam;
angle01<-batchgeneangles(pam,batch1)$first_angle
angle02<-batchgeneangles(pam,batch1)$second_angle

(ggplot()+aes(x=-50:50)+
  geom_point(aes(y=c(leftward1 %>% rev,angle01,rightward1) %>% abs %>% acos/pi),colour='red')+
  geom_point(aes(y=c(leftward2 %>% rev,angle02,rightward2) %>% abs %>% acos/pi),colour='blue')) %>% plot

# par(mfrow=c(1,2))
# c(leftward1 %>% rev,angle01,rightward1) %>% plot
# c(leftward2 %>% rev,angle02,rightward2) %>% plot