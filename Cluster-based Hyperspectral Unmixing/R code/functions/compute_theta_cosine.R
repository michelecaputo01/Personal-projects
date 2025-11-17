compute_theta_cosine = function(data, labels){
  
  k = max(labels)
  
  mean_tot = colMeans(data)
  
  mean_clusters = matrix(0, k, dim(data)[2])
  
  for(i in 1:k){
    
    ind = which(labels == i)
    
    mean_clusters[i,] = colMeans(data[ind,])
    
  }
  
  mat_diff = c(proxy::dist(t(mean_tot), mean_clusters, method = "cosine"))^2
  
  B = sum(mat_diff * table(factor(labels, levels = seq(1,k))), na.rm = TRUE)
  
  W = 0
  
  for(i in 1:k){
    
    ind = which(labels == i)
    
    if(length(ind) > 0){
      
      dist = c(proxy::dist(t(mean_clusters[i,]), data[ind,], method = "cosine"))^2
      W = W + sum(dist)
      
    }
    
  }
  
  theta = B / (B + W)
  
  return(theta)
  
}