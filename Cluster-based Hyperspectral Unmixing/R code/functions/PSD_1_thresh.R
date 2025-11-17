PSD_1_thresh = function(data, gridSizeX, gridSizeY, t, method = "l2"){
  
  numNodes = gridSizeX*gridSizeY
  
  grid = expand.grid(x = 1:gridSizeX, y = 1:gridSizeY)
  
  if(method == "cosine"){
    
    distance = proxy::dist(data, method = "cosine")
    
  } else if (method == "arcosine"){
    
    distance = acos(-proxy::dist(data, method = "cosine") + 1)
    
  } else if (method == "l2"){
    
    distance = as.matrix(dist(data, method = "euclidean"))
    
  }
  
  geo_dist = as.matrix(dist(grid))
  
  dissimilarity = matrix(1e10, nrow = numNodes, ncol = numNodes)
  diag(dissimilarity) = 0
  
  for(i in 1:(numNodes-1)){
    
    cat(sprintf("\rProgress: Threshold = %f, progress %d/%d       ", t, i+1, numNodes))
    flush.console()
    
    neigh = as.numeric(which(geo_dist[i,] < t))
    neigh = neigh[neigh > i]
    
    for(j in neigh){
      
      path = direct_path(i, j, grid, gridSizeX)
      
      dissimilarity[i,j] = dissimilarity[j,i] = sum(distance[cbind(path[-length(path)], path[-1])])
      
    }
    
  }
  
  diss_with_threshold = as.dist(dissimilarity)
  
  return(diss_with_threshold)
  
}