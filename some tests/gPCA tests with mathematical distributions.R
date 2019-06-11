library(ExpressionAtlas)
library(magrittr)
library(ggplot2)

### generation of data following the same distribution, gPCA performed with artificial batches
A<-rnorm(6000*100) %>% matrix(ncol=100)#~100 samples, 6000 genes
batch<-runif(100,0,2) %>% floor
g<-A %>% gPCA(batch,scale=TRUE)
g %>% viz_gpca
g %>% viz_gpca(guided=FALSE)
g %>% viz_gpca_contrib

Asharp<-A %>% cbind(rnorm(6000,sd=5))
batchsharp<-batch %>% c(runif(1,0,2) %>% floor)
gs<-Asharp %>% gPCA(batchsharp,scale=TRUE)
gs %>% viz_gpca_contrib

sdgeo<-100000
for(i in 1:50) sdgeo%<>%c(sdgeo[i]*.75)
Ageo<-rnorm(6000*100,sd=sdgeo) %>% matrix(ncol=100)
batch<-runif(100,0,2) %>% floor
ggeo<-Ageo %>% gPCA(batch)
ggeo %>% viz_gpca_contrib(end='all')

Adouble<-Ageo[1:2000,] %>% rbind(.,.)
batch<-runif(100,0,2) %>% floor
gdouble<-Adouble %>% gPCA(batch,scale=TRUE)
gdouble %>% viz_gpca_contrib(end='all')

gr<-raw %>% gPCA(x)
gr %>% viz_gpca_contrib(end=50)

grs<-raw %>% gPCA(x,scale=TRUE)
grs %>% viz_gpca_contrib(end=50)

### generation of data following 5 different distributions (5 batches), gPCA performed according to these 5 batches
batch<-runif(100,0,2) %>% floor
B<-NULL
for(i in 1:100) B%<>%cbind(rnorm(6000,mean=5+batch[i]))
gb<-B %>% gPCA(batch)
gb %>% viz_gpca(guided=FALSE)
gb %>% viz_gpca_contrib(log1p)

library(sva)
Bcombat<-ComBat(B,batch)
gbcombat<-Bcombat %>% gPCA(batch)
gbcombat %>% viz_gpca(guided=FALSE,dims=c(100,99))
gbcombat %>% viz_gpca
gbcombat %>% viz_gpca_contrib(end='all')

library(pamr)
pamobject<-list(x=B,batchlabels=batch)
Bpam<-pamr.batchadjust(pamobject)$x
gbpam<-Bpam %>% gPCA(batch)
gbpam %>% viz_gpca(guided=FALSE)
gbpam %>% viz_gpca
gbpam %>% viz_gpca_contrib(end='all')

library(bapred)
Bra<-ratioa(B %>% t,batch %>% as.numeric %>% add(1) %>% factor)$xadj %>% t
gbra<-Bra %>% gPCA(batch)
gbra %>% viz_gpca(guided=FALSE)
gbra %>% viz_gpca
gbra %>% viz_gpca_contrib(end='all')

library(bapred)
Brg<-ratiog(B %>% t,batch %>% as.numeric %>% add(1) %>% factor)$xadj %>% t
gbrg<-Brg %>% gPCA(batch)
gbrg %>% viz_gpca(guided=FALSE)
gbrg %>% viz_gpca
gbrg %>% viz_gpca_contrib(end='all')


### same data (5 batches), gPCA performed according to 5 fake batches
fake.batch<-runif(25,0,5) %>% floor
gbf<-B %>% gPCA(fake.batch)
gbf %>% viz_gpca_contrib(log10,end=10)

### same data (5 batches), gPCA performed according to 2 fake batches
fake.batch.2<-runif(324,0,2) %>% floor
gbf2<-B %>% gPCA(fake.batch.2)
gbf2 %>% viz_gpca_contrib(end=10)

### generated data following 2 different distributions (2 batches), gPCA performed according to these 2 batches
batch2<-rbinom(324,1,.7)
C<-NULL
for(i in 1:324) C%<>%cbind(rnorm(29912,batch2[i]))
gc<-C %>% gPCA(batch2)
gc %>% viz_gpca_contrib

### same data, gPCA performed according to 2 fake batches
fake.batch2<-rbinom(324,1,.7)
gcf<-C %>% gPCA(fake.batch2)
gcf %>% viz_gpca_contrib

### generated data following 3 different distributions (2 batches), gPCA performed according to these 3 batches
batch3<-rbinom(324,2,.7)
D<-NULL
for(i in 1:324) D%<>%cbind(rnorm(29912,batch3[i]))
gd<-D %>% gPCA(batch3)
gd %>% viz_gpca_contrib

### same data, gPCA performed according to 3 fake batches
fake.batch3<-rbinom(324,2,.7)
gdf<-D %>% gPCA(fake.batch3)
gdf %>% viz_gpca_contrib(log1p)

### 2 batches again but this time, less difference between the two distributions (the mean difference is a sparse vector)
sparse<-rbinom(6000,1,.01)
batch2<-rbinom(100,1,.7)
E<-NULL
for(i in 1:100) E%<>%cbind(rnorm(6000,batch2[i]*sparse))
ge<-E %>% gPCA(batch2)
ge %>% viz_gpca
ge %>% viz_gpca(guided=FALSE)
ge %>% viz_gpca_contrib#(log1p)
## fake batches
fake.batch2<-rbinom(100,1,.7)
gef<-E %>% gPCA(fake.batch2)
gef %>% viz_gpca
gef %>% viz_gpca(guided=FALSE)
gef %>% viz_gpca_contrib

### sparse with 3 batches
sparse<-rbinom(29912,1,.01)
G<-NULL
for(i in 1:324) G%<>%cbind(rnorm(29912,batch3[i]*sparse))
gg<-G %>% gPCA(batch3)
gg %>% viz_gpca_contrib(log1p)
## fake batches
ggf<-G %>% gPCA(fake.batch2)
ggf %>% viz_gpca_contrib(log1p)


### test with genetic data, no batch
setwd('c:/users/medion/desktop/Comparison-of-batch-effect-correction-methods/data')
raw1<-get(load('data/E-MTAB-4344-atlasExperimentSummary.Rdata'))$rnaseq %>% assays %>% use_series(counts)
fake.batch<-rbinom(raw1 %>% ncol,1,.6)
gr1<-raw1 %>% gPCA(fake.batch,scale=TRUE)
gr1 %>% viz_gpca_contrib
gr1 %>% viz_gpca
gr1 %>% viz_gpca(guided=FALSE)

### test with genetic data, no batch, variance filtering on genes
#observation of the variance gene by gene
ggplot()+geom_density(aes(x=raw1 %>% t %>% scale %>% colVars %>% log1p))+ggtitle('distribution of log(1+variance_gene)')
#
variance.filter<-(raw1 %>% rowVars %>% log1p<.07)&(raw1 %>% rowVars %>% log1p>.03)
raw1ft<-raw1[variance.filter,]
ggplot()+geom_density(aes(x=raw1ft %>% rowVars %>% log1p))
fake.batch<-rbinom(raw1ft %>% ncol,1,.6)
grft<-raw1ft %>% gPCA(fake.batch)
grft %>% viz_gpca_contrib#(log1p)


setwd('data')
atlas2<-list(mtab3716=get(load('E-MTAB-3716-atlasExperimentSummary.Rdata')),
             mtab4344=get(load('E-MTAB-4344-atlasExperimentSummary.Rdata')))
raw2<-NULL
for(mtab in atlas2) raw2%<>%cbind(mtab$rnaseq %>% assays %>% use_series(counts))
x2<-NULL
for(i in atlas2 %>% seq_along) x2%<>%c(i %>% rep(dim(atlas2[[i]]$rnaseq)[2]))
gr2<-raw2 %>% gPCA(x2,scale=TRUE,nperm=1000)
gr2 %>% viz_gpca_contrib
gr2 %>% viz_gpca
gr2 %>% viz_gpca(guided=FALSE)

fake.x2<-x2 %>% sample
gr2f<-raw2 %>% gPCA(fake.x2,scale=TRUE,nperm=1000)
gr2f %>% viz_gpca_contrib#(log1p)
gr2f %>% viz_gpca
gr2f %>% viz_gpca(guided=FALSE)


graw<-raw %>% gPCA(x)
graw %>% viz_gpca_contrib
graw %>% viz_gpca
graw %>% viz_gpca(guided=FALSE)


gpam<-pam %>% gPCA(x)
gpam %>% viz_gpca_contrib
gpam %>% viz_gpca
gpam %>% viz_gpca(guided=FALSE)


gcombat<-combat %>% gPCA(x)
gcombat %>% viz_gpca_contrib
gcombat %>% viz_gpca
gcombat %>% viz_gpca(guided=FALSE)

gccombat<-ccombat %>% gPCA(x)
gccombat %>% viz_gpca_contrib
gccombat %>% viz_gpca
gccombat %>% viz_gpca(guided=FALSE)

gra<-ra %>% gPCA(x)
gra %>% viz_gpca_contrib
gra %>% viz_gpca
gra %>% viz_gpca(guided=FALSE)

grg<-rg %>% gPCA(x)
grg %>% viz_gpca_contrib
grg %>% viz_gpca
grg %>% viz_gpca(guided=FALSE)

gruvg<-ruvg %>% gPCA(x)
gruvg %>% viz_gpca_contrib
gruvg %>% viz_gpca
gruvg %>% viz_gpca(guided=FALSE)

library(pamr)
pamobject<-list(
  'x'=raw2,
  'batchlabels'=x2
)
pam <- pamr.batchadjust(pamobject)$x
gpam<-pam %>% gPCA(x2)
gpam %>% viz_gpca_contrib
gpam %>% viz_gpca
gpam %>% viz_gpca(guided=FALSE)

library(sva)
combat<-ComBat(raw2,x2 %>% factor)

library(bapred)
ra<-ratioa(raw2 %>% t, x2 %>% as.numeric %>% factor)$xadj %>% t
gra<-ra %>% gPCA(x2,scale=TRUE)
gra %>% viz_gpca
gra %>% viz_gpca(guided=FALSE)
gra %>% viz_gpca_contrib



gr<-raw %>% gPCA(x)
gr %>% viz_gpca_contrib(log1p)
grf<-raw %>% gPCA(batch)
grf %>% viz_gpca_contrib(log1p)

setwd('data')
atlas5<-list(mtab2836=get(load('E-MTAB-2836-atlasExperimentSummary.Rdata')),
             mtab3716=get(load('E-MTAB-3716-atlasExperimentSummary.Rdata')),
             mtab4344=get(load('E-MTAB-4344-atlasExperimentSummary.Rdata')),
             mtab513=get(load('E-MTAB-513-atlasExperimentSummary.Rdata')),
             mtab3871=get(load('E-MTAB-3871-atlasExperimentSummary.Rdata')))
raw5<-NULL
for (i in names(atlas5)) {
  expAcc <- i
  k <- atlas5[[i]]
  exp <- k$rnaseq
  eCounts <- assays(exp)$counts
  samples<-colnames(eCounts)
  average.counts<-technical_replicate_average_gtex(exp,expAcc)
  raw5 <- cbind(all,average.counts)
}
#for(mtab in atlas5) raw5%<>%cbind(mtab$rnaseq %>% technical_replicate_average_gtex('mtab'))
  #raw5%<>%cbind(mtab$rnaseq %>% assays %>% use_series(counts))
x5<-NULL
for(i in atlas5 %>% seq_along) x5%<>%c(i %>% rep(dim(atlas5[[i]]$rnaseq)[2]))
gr5<-raw5 %>% gPCA(x5,scale=TRUE)
gr5 %>% viz_gpca_contrib
gr5 %>% viz_gpca
gr5 %>% viz_gpca(guided=FALSE)



combat<-ComBat(raw5,x5)
