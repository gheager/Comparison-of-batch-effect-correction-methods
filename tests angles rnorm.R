nsamples<-200
ngenes<-1700
fake_batch<-c(0,3,4,20,100) %>% sample(nsamples,TRUE)
batch<-rep(0,5) %>% sample(nsamples,TRUE)
tissue<-c('brain','liver','lung') %>% sample(nsamples,TRUE)
brain.profile<-0 %>% rep(ngenes) %>% inset(1:ngenes %>% sample(5),10)
liver.profile<-0 %>% rep(ngenes) %>% inset(1:ngenes %>% sample(5),20)
lung.profile<-0 %>% rep(ngenes) %>% inset(1:ngenes %>% sample(5),-10)
tissue.profile<-cbind(brain.profile,liver.profile,lung.profile) %>% set_colnames(c('brain','liver','lung'))
tissues<-tissue.profile[,tissue]
data<-rnorm(ngenes*nsamples,mean=(batch+t(tissues)) %>% t) %>% matrix(nrow=ngenes)

data %>% eigenangles(fake_batch,tissue)->ang
viz_angles(ang=ang,filename = 'viz math fake batch.pdf')
