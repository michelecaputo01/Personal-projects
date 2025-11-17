sum_distance_matrix = function(data, grid, contribution, alpha = 0.5){
  
  s_mean = colMeans(grid)
  sigma2_s = mean(rowSums((as.matrix(grid) - matrix(s_mean, nrow(grid), ncol(grid), byrow = TRUE))^2))
  space_contr = as.matrix(dist(grid, method = "euclidean", diag = TRUE, upper = TRUE))^2
  space_contr = space_contr/sigma2_s
  
  if(contribution == "l2"){
    
    f_mean = colMeans(data)
    sigma2_f = mean(rowSums((as.matrix(data) - matrix(f_mean, nrow(data), ncol(data), byrow = TRUE))^2))
    funct_contr = as.matrix(dist(data, method = "euclidean", diag = TRUE, upper = TRUE))^2
    funct_contr = funct_contr/sigma2_f
    
  } else if (contribution == "cosine"){
    
    f_mean = colSums(data)
    f_mean = f_mean/sqrt(sum(f_mean^2))
    sigma2_f = mean(proxy::dist(data, t(f_mean), method = "cosine")^2)
    funct_contr = as.matrix(proxy::dist(data, method = "cosine"))^2
    funct_contr = funct_contr/sigma2_f
    
  }
  
  dissimilarity = sqrt(space_contr * alpha + funct_contr * (1-alpha))
  
  plot(space_contr[-1,1]*alpha, dissimilarity[-1,1]^2, xlab = "Normalized spatial distance", 
       ylab = "Sum distance", main = "Comparison w.r.t. the first observation",
       pch = 19, col = "grey", lwd = 1, asp = 1)
  grid()
  lines(c(-10,50), c(-10,50), lwd = 3, col = "black", lty = 2)
  
  return(dissimilarity)
  
}