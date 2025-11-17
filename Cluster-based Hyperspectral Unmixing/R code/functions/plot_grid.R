plot_grid = function(gridSizeX, gridSizeY, labels, title = "", scaled = 0, 
                     centroids = FALSE, data, title_centr = "", 
                     original = FALSE, original_data){
  
  # labels as vector starting from bottomleft (one row at time)
  
  old_par = par(no.readonly = TRUE)
  
  k = max(labels)
  
  matrix_labels = matrix(labels, nrow = gridSizeX, ncol = gridSizeY, byrow = FALSE)
  
  if(scaled){
    
    scale = colorRampPalette(c("white", "red"))
    colors = scale(scaled)
    
  } else if(!scaled){
    
    colors = c("red", "gold", "green", "skyblue", "black", "violet", "blue", "orange", "grey", "forestgreen")
    colors = rep(colors, length.out = k)
    
  }
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(8, 1))
  
  par(oma = c(0,0,0,0) + 1.5)
  par(mar = c(1,1,2,1))
  par(las = 1)
  
  if(scaled){
    
    image(x = 1:gridSizeX, y = 1:gridSizeY, matrix_labels, 
          asp = 1, col = colors, xlab = "", ylab = "", zlim = c(0,1))
    
    image(1, seq(0, 1, length=100), t(matrix(1:100, ncol = 1)),
          col = colors, axes = FALSE, xlab = "", ylab = "")
    axis(4)
    box()
    
  } else if(!scaled){
    
    image(x = 1:gridSizeX, y = 1:gridSizeY, matrix_labels, 
          asp = 1, col = colors, xlab = "", ylab = "")
    
    image(1, seq(1, max(matrix_labels), length = k), t(1:k),
          col = colors, axes = FALSE, xlab = "", ylab = "")
    axis(4, at = 1:k, 1:k)
    box()
    
  }
  
  mtext(title, outer = TRUE, cex = 1.5, line = -1, font = 2)
  
  par(old_par)
  
  par(las = 1)
  
  if(centroids){
    
    repr = matrix(0, nrow = k, ncol = dim(data)[2])
    
    for(i in 1:k){
      
      ind = which(labels == i)
      
      if(length(ind) > 0){ repr[i,] = colMeans(data[ind,]) }
      
    }
    
    matplot(t(data[sample(1:nrow(data), size = 1000),]), type = 'l', lty = 1, ylab = "", col = "lightgrey")
    matlines(t(repr), type = 'l', lty = 1, col = colors, lwd = 3)
    grid()
    mtext(title_centr, cex = 1.5, line = 1, font = 2)
    
  }
  
  par(old_par)
  
  par(las = 1)
  
  if(original){
    
    repr = matrix(0, nrow = k, ncol = dim(original_data)[2])
    
    for(i in 1:k){
      
      ind = which(labels == i)
      
      if(length(ind) > 0){ repr[i,] = colMeans(original_data[ind,]) }
      
    }
    
    matplot(t(original_data[sample(1:nrow(original_data), size = 1000),]), type = 'l', lty = 1, ylab = "", col = "lightgrey")
    matlines(t(repr), type = 'l', lty = 1, col = colors, lwd = 3)
    grid()
    mtext(title_centr, cex = 1.5, line = 1, font = 2)
    
  }
  
  par(old_par)

}


