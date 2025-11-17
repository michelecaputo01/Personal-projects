compute_local_representative = function(data, grid_points, seed_points){
  
  cross_dist = dist(seed_points[,c(1,2)])
  sigma2 = (max(cross_dist)/min(cross_dist))^2
  cov_mat = diag(rep(sigma2,2))
  
  w = numeric(dim(grid_points)[1])
  
  d = as.matrix(data)
  
  weighted_fun = NULL
  
  for(i in 1:dim(seed_points)[1]){
    
    ind = which(grid_points$cell == i)
    
    w = dmvnorm(as.matrix(grid_points[ind,c(1,2)]), as.matrix(seed_points[i,]), cov_mat)
    
    weighted_fun = rbind(weighted_fun, colWeightedMeans(d[ind, , drop = FALSE], w = w))
    
  }
  
  return(as.data.frame(weighted_fun))
}