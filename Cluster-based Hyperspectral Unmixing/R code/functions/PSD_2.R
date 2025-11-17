PSD_2 = function(data, gridSizeX, gridSizeY, method = "l2"){
  
  numNodes = gridSizeX*gridSizeY
  
  grid = expand.grid(x = 1:gridSizeX, y = 1:gridSizeY)
  
  if(method == "cosine"){
    
    distance = proxy::dist(data, method = "cosine")
    
  } else if (method == "arcosine"){
    
    distance = acos(-proxy::dist(data, method = "cosine") + 1)
    
  } else if (method == "l2"){
    
    distance = as.matrix(dist(data, method = "euclidean"))
    
  }
  
  dissimilarity = matrix(0, nrow = numNodes, ncol = numNodes)

  for(i in 1:(numNodes-1)){
    
    cat(sprintf("\rProgress: %d/%d            ", i+1, numNodes))
    flush.console()
    
    for(j in (i+1):numNodes){
      
      path = direct_path(i, j, grid, gridSizeX)
      
      dissimilarity[i,j] = dissimilarity[j,i] = 0.5*(sum(distance[cbind(path[1], path[-1])]) + sum(distance[cbind(path[length(path)], path[-length(path)])]))
      
    }
    
  }
  
  dissimilarity = as.dist(dissimilarity)
  
  return(dissimilarity)
  
}