data_matrix <- filtered#matrix(rnorm(1000), 20, 50)
#batch <- c(rep(1, 20), rep(2, 20), rep(3, 10))
#outs <- bbknn(data_matrix, batches,pca=FALSE,compute_pca=FALSE)

bbknn<-function(data_matrix, batch, pca = TRUE, compute_pca = "python", 
          nPcs = NULL) 
{
  if (!is.matrix(data_matrix)) {
    warning("matrix expected for data_matrix")
    data_matrix <- as.matrix(data_matrix)
  }
  if (is.null(nPcs)) 
    nPcs <- min(50, nrow(data_matrix), ncol(data_matrix))
  if (nPcs > nrow(data_matrix)) {
    warning("number of genes less than nPcs")
    print(paste("using", nrow(data_matrix), "components"))
    nPcs <- nrow(data_matrix)
  }
  anndata <- reticulate::import("anndata", convert = FALSE)
  bbknn <- reticulate::import("bbknn", convert = FALSE)
  sc <- reticulate::import("scanpy.api", convert = FALSE)
  if (is.character(batch)) 
    batch <- as.factor(batch)
  if (is.factor(batch)) 
    batch <- as.numeric(batch)
  if (is.numeric(batch)) 
    batch <- as.integer(batch)
  if (pca) {
    if (compute_pca == "python") {
      pca <- sc$pp$pca(t(data_matrix))
    }
    else if (compute_pca != "python") {
      print("test")
      pca <- reticulate::r_to_py(prcomp(t(data_matrix))$x[1:nPcs, 
                                                         ])
    }
    adata <- anndata$AnnData(X = pca, obs = batch)
    sc$tl$pca(adata)
    adata$obsm$X_pca <- pca
  }
  else {
    adata <- anndata$AnnData(X = t(data_matrix), obs = batch)
    sc$tl$pca(adata)
  }
  bbknn$bbknn(adata, batch_key = 0)
  corrected_matrix <- t(reticulate::py_to_r(adata$X))
  sc$tl$pca(adata)
  pca <- reticulate::py_to_r(adata$obsm$X_pca)
  sc$tl$tsne(adata)
  tsne <- reticulate::py_to_r(adata$obsm$X_tsne)
  sc$tl$umap(adata)
  umap <- reticulate::py_to_r(adata$obsm$X_umap)
  output <- list(corrected_matrix = corrected_matrix, pca = pca, 
                 tsne = tsne, umap = umap)
  return(output)
}
