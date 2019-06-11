# library(bbknn)
# outs <- bbknn(t(filtered),batch,compute_pca="python",nPcs=17)
library(reticulate)
use_python("python3")
anndata<-import('anndata',convert=FALSE)
bbknn<-import('bbknn',convert=FALSE)
sc<-import('scanpy.api',convert=FALSE)
#sp<-import('scipy',convert=FALSE)


filtered %>% t %>% prcomp %>% use_series(x) %>% r_to_py -> pca
adata<-anndata$AnnData(X=pca,obs=batch %>% r_to_py)
sc$tl$pca(adata)
umap0<-sc$tl$umap(adata)
#adata$obsm$X_pca <- pca
bbknn$bbknn(adata,batch_key=0)
#b=bbknn$bbknn_pca_matrix(pca,batch_list = batch)
sc$tl$umap(adata)
umap = py_to_r(adata$obsm$X_umap)
ggplot()+aes(
  x=umap[,1],
  y=umap[,2],
  colour=tissue
)+geom_point()
