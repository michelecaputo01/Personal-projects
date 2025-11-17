compute_theta = function(data, labels){
  
  k = max(labels)
  
  mean_tot = colMeans(data)
  
  mean_clusters = matrix(0, k, dim(data)[2])
  
  for(i in 1:k){
    
    ind = which(labels == i)
    
    mean_clusters[i,] = colMeans(data[ind,])
    
  }
  
  mat_diff = mean_clusters - matrix(mean_tot, nrow = k, ncol = dim(data)[2], byrow = TRUE)
  
  B = sum(rowSums(mat_diff^2) * table(factor(labels, levels = seq(1,k))), na.rm = TRUE)
  
  W = 0
  
  for(i in 1:k){
    
    ind = which(labels == i)
    
    if(length(ind) > 0){
      
      w = data[ind,]
      diff_sq = (w - matrix(mean_clusters[i,], nrow = length(ind), ncol = dim(data)[2], byrow = TRUE))^2
      W = W + sum(diff_sq)
      
    }
    
  }
  
  theta = B / (B + W)
  
  return(theta)
  
}