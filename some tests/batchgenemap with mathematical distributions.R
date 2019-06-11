ngenes=1743
A<-rnorm(ngenes*14,mean=5*pcb1$sdev %>% rep(each=ngenes),sd=pcb1$sdev) %>% matrix(nrow=ngenes)
batchgeneladder(A,batch1,abs)
batchgenemap(A,batch1)
batchgenemap(A,batch1,xshift=c(0,0),yshift=c(0,0))
A %>% t %>% prcomp -> pca3
A[,batch1=='geod'] %>% t %>% prcomp -> pca1
A[,batch1=='mtab'] %>% t %>% prcomp -> pca2
angle(pca1$rotation[,1],pca2$rotation[,1])
angle(pca1$rotation[,1],pca3$rotation[,1])
angle(pca2$rotation[,1],pca3$rotation[,1])

pca1$rotation[,1] %>% order %>% head
pca2$rotation[,1] %>% order %>% head
cor.test(pca1$rotation[,1],pca2$rotation[,1],method='kendall')
cor.test(pca3$rotation[,1],pca2$rotation[,1],method='kendall')
cor.test(pca1$rotation[,1],pca3$rotation[,1],method='kendall')

ggplot()+aes(
  x=pcb1$rotation[,1],
  y=pcb2$rotation[,1]
)+geom_point()+geom_smooth(method='lm')
