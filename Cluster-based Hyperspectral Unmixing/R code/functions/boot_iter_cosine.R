boot_iter_cosine = function(data, n, k, grid_points, print = FALSE){
  
  ############################ VORONOI TESSELATION ####################
  
  seed_points = grid_points[sample(nrow(grid_points), n),]
  
  distances = dist(grid_points[,c(1,2)], seed_points)
  grid_points$cell = apply(distances, 1, which.min)
  
  if(print){
    tess = deldir(seed_points$x, seed_points$y)
    plot(tess, wlines = "tess")
    title("Voronoi tesselation")
    points(seed_points, pch = 20, col = "red")
    points(grid_points$x, grid_points$y, pch = 1, col = as.factor(grid_points$cell))
  }
  
  ############################ COMPUTE LOCAL REPRESENTATIVES ####################
  
  local_representatives = compute_local_representative(data, grid_points, seed_points)
  
  if(print){
    random_cell = sample(1:n, size = 1)
    matplot(t(data[which(grid_points$cell == random_cell),]), type = 'l', col = 'black', lwd = 1, lty = 2)
    ind = as.numeric(rownames(seed_points[random_cell,]))
    matlines(t(data[ind,]), type = 'l', col = 'forestgreen', lwd = 3, lty = 1)
    matlines(t(local_representatives[random_cell,]), type = 'l', col = 'red', lwd = 3, lty = 1)
    legend("bottomright", legend = c("Observations", "Local center", "Local Representative"), col = c("black", "forestgreen", "red"), lty = c(2,1,1), lwd = c(1,3,3))
  }
  
  ######################### CLUSTERING #################
  
  dissimilarity = proxy::dist(local_representatives, method = "cosine")
  
  cluster_assignment = as.numeric(pam(dissimilarity, k, diss = TRUE, nstart = 30, cluster.only = TRUE))

  single_cluster = cluster_assignment[grid_points$cell]
  
  order = as.numeric(names(sort(table(single_cluster), decreasing = TRUE)))
  
  cl = numeric(length(single_cluster))
  
  for(i in order){
    ind = which(single_cluster == i)
    cl[ind] = which(order == i)
  }
  
  return(cl)
  
}
