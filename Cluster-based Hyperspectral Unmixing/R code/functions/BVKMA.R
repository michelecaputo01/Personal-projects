BVKMA = function(data, gridSizeX, gridSizeY, n, k, B, print = FALSE){
  
  grid_points = expand.grid(x = 1:gridSizeX, y = 1:gridSizeY)
  
  cluster_matrix = {}
  
  cat("\nBootstrap Simulations ...\n")
  
  for(b in 1:B){
    
    cat(sprintf("\rProgress: %d/%d", b, B))
    flush.console()
    
    local_cluster = boot_iter(data, n, k, grid_points, print)
    
    cluster_matrix = rbind(cluster_matrix, local_cluster)
    
  }  ### END BOOTSTRAP
  
  rownames(cluster_matrix) = 1:B
  
  ################################################################################
  ############################ CLUSTER MATCHING ##################################
  ################################################################################
  
  perm = ecr.iterative.1(cluster_matrix, k)$permutations
  
  for(i in 1:B){
    cluster_matrix[i,] = perm[i,][cluster_matrix[i,]]
  }

  ################################################################################
  ############################## FINAL RESULTS ###################################
  ################################################################################
  
  final_clusters = rep(0, dim(data)[1])
  
  for (i in 1:dim(data)[1]){
    final_clusters[i] = sort(unique(cluster_matrix[,i]))[which.max(table(cluster_matrix[,i]))]
  }
  
  ######################### PLOT CLUSTERS #################
  
  if(print){
    ind = sample(1:dim(data)[1], size = 100)
    matplot(t(data[ind,]), col = as.factor(final_clusters[ind]), type = 'l', lty = 1)
  }
  
  plot_grid(gridSizeX, gridSizeY, final_clusters, title = paste("Final clusterization with n =", n, "and k =", k), 
            scaled = 0, centroids = FALSE, original = FALSE)
  
  ######################### FREQUENCIES #################
  
  frequencies = rep(0, dim(data)[1])
  
  for(i in 1:dim(data)[1]){
    frequencies[i] = sum(as.numeric(cluster_matrix[,i] == final_clusters[i]))
  }
  
  ################################################################################
  ########################## EVALUATION OF ETA ###################################
  ################################################################################
  
  eta_values = compute_eta(cluster_matrix)
  
  eta_x = eta_values$eta_x
  
  eta = eta_values$eta   # lower --> better
  
  plot_grid(gridSizeX, gridSizeY, eta_x, title = paste("Entropy with", k, "clusters and n =", n),
            scaled = 100, centroids = FALSE, original = FALSE)
  
  ################################################################################
  ########################## EVALUATION OF THETA #################################
  ################################################################################
  
  theta = compute_theta(data, final_clusters)
  
  cat(paste("\nDone for n =", n, "and k =", k))
  cat("\n--------------------------\n")
  
  ################################################################################
  ############################## RESULTS #########################################
  ################################################################################
  
  results = list( cluster = final_clusters,
                  theta = theta,
                  eta = eta,
                  boot_clusters = cluster_matrix,
                  frequencies = frequencies/B )
  
  return(results)
}
