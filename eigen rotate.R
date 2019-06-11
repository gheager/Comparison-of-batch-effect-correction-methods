inner<-function(u,v) (t(as.matrix(u))%*%as.matrix(v))[1]

normalise<-function(u) u %>% as.matrix %>% divide_by(norm(.,'F'))

# orth<-function(V){
#   V%<>%as.matrix
#   if(ncol(V)==1){
#     return(normalise(V))
#   }else{
#     ncol(V) -> n
#     V[,-n] %>% orth -> U
#     V[,n] %>% matrix -> u
#     for(i in 1:(n-1)){
#       u%<>%subtract(inner(u,U[,i]) * U[,i])
#     }
#     u%<>%normalise
#     cat(n %>% paste('\t'))
#     return(cbind(U,u))
#   }
# }

orth<-function(V){
  V%<>%as.matrix
  ncol(V)->n
  V[,1] %>% matrix %>% norm('F') %>% c(0 %>% rep(n-1)) %>% matrix -> coord
  #1 %>% c(0 %>% rep(n-1)) %>% matrix -> coord
  V[,1]%<>%normalise
  if(n>1){
    for(i in 2:n){
      coord %<>% cbind((t(V[,1:(i-1)])%*%V[,i]) %>% c(1,0 %>% rep(n-i)))
      V[,i]%<>%subtract(V[,1:(i-1)]%*%t(t(V[,i])%*%V[,1:(i-1)]))
      coord[i,i] %<>% multiply_by(V[,i] %>% matrix %>% norm('F'))
      V[,i]%<>%normalise
    }
  }
  return(V %>% structure(coord=coord))
}
# orth<-function(V){
#   V%<>%as.matrix
#   ncol(V)->n
#   V[,1]%<>%normalise
#   for(i in 2:n){
#     for(j in 1:(i-1)){
#       V[,i]%<>%subtract(V[,j]%*%t(t(V[,i])%*%V[,j]))
#     }
#     V[,i]%<>%normalise
#   }
#   return(V)
# }

proj<-function(u,V){
  V%<>%orth
  v<-0
  for(i in 1:ncol(V)){
    v%<>%add(inner(u,V[,i])*V[,i])
  }
  return(v)
}