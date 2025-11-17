generate_data = function(gridSizeX, gridSizeY, natt, sigma_f, scale, minshift, maxshift, print = TRUE){
  
  par(las = 1)
  
  grid = expand.grid(x = 1:gridSizeX, y = 1:gridSizeY)
  
  numNodes = gridSizeX * gridSizeY
  
  center_points = expand.grid(x = c(round(gridSizeX/4), round(gridSizeX*3/4)), 
                              y = c(round(gridSizeY/4), round(gridSizeY*3/4)))
  
  distances = as.matrix(dist(rbind(center_points, grid)))[-seq(1,4), seq(1,4)]
  rownames(distances) = 1:numNodes
  cluster = numeric(numNodes)
  
  for(i in 1:numNodes){
    cluster[i] = ifelse(sum(distances[i,] < min(gridSizeX, gridSizeY)/5) > 0, 
                        which(distances[i,] < min(gridSizeX, gridSizeY)/5), 5)
    if(cluster[i] == 1){ cluster[i] = 4 }
    cluster[i] = cluster[i] - 1
  }
  
  if(print){
    plot_grid(gridSizeX, gridSizeY, cluster, title = paste0("True labels in simulation - ", gridSizeX, "x", gridSizeY), 
              scaled = 0, centroids = FALSE, original = FALSE)
  }

  seq_of_points = seq(0, 5, length.out = natt)
  
  f = list()
  
  f[[1]] = function(x){return(-3*cos(x) + sin(x^2))}
  f[[2]] = function(x){return(3.5*cos(x-1)^3)}
  f[[3]] = function(x){return(cos(5*x)*2.5)}
  f[[4]] = function(x){return(rep(0.1, length(x)))}
  
  colors = c("red", "gold", "green", "skyblue")
  
  if(print){
    plot(1:natt, f[[1]](seq_of_points), type = 'l', col = colors[1], main = "Generator functions", ylab = '',
         xlab = '', lwd = 2, ylim = c(-4,4))
    lines(1:natt, f[[2]](seq_of_points), lwd = 2, col = colors[2])
    lines(1:natt, f[[3]](seq_of_points), lwd = 2, col = colors[3])
    lines(1:natt, f[[4]](seq_of_points), lwd = 2, col = colors[4])
    grid()
  }
  
  data = matrix(0, nrow = numNodes, ncol = natt)
  
  dist = outer(seq_of_points, seq_of_points, "-")^2
  
  cov = sigma_f^2 * exp(-0.5 / scale^2 * dist)
  
  for(i in 1:4){
    
    ind = which(cluster == i)
    
    for(j in ind){
      data[j,] = mvrnorm(1, mu = f[[i]](seq_of_points), Sigma = cov) + 
        runif(1, minshift, maxshift)
    }
    
  }
  
  if(sum(data <= 0) > 0){
    data = data + 0.5 + abs(min(data))
  }
  
  data = as.data.frame(data)
  
  if(print){
    
    for(k in 1:4){
      matplot(t(data[sample(numNodes, size = numNodes/4),]), col = "grey", ylim = c(min(data) - 0.5, max(data) + 0.5), main = paste("Cluster", k),
              type = 'l', lty = 1, xlab = '', ylab = '')
      matlines(t(data[which(cluster == k),]), type = 'l', lwd = 1, lty = 1, col = colors[k])
      grid()
    }
    
    ind = NULL
    for(k in 1:4){
      ind = c(ind, sample(which(cluster == k), size = 3))
    }
    
    matplot(1:natt, t(data[ind,]), type = 'l', lty = 1, col = "black", xlab = '', ylab = '', ylim = c(min(data) - 0.5, max(data) + 0.5),
            main = "Examples of generated functions (three for each cluster)", lwd = 1.5)
    grid()
    
    matplot(1:natt, t(data[ind,]), type = 'l', lty = 1, col = colors[cluster[ind]], xlab = '', ylab = '', ylim = c(min(data) - 0.5, max(data) + 0.5),
            main = "Examples of generated functions (three for each cluster)", lwd = 1.5)
    grid()
    
    matplot(1:natt, t(data), type = 'l', lty = 1, col = "grey", xlab = '', ylab = '', ylim = c(min(data) - 0.5, max(data) + 0.5),
            main = "Simulated data")
    grid()
    
    matplot(1:natt, t(data), type = 'l', lty = 1, col = colors[cluster], xlab = '', ylab = '', ylim = c(min(data) - 0.5, max(data) + 0.5),
            main = "Simulated data")
    grid()
    
  }
  
  return(list(data = data,
              real_labels = cluster))
  
}
