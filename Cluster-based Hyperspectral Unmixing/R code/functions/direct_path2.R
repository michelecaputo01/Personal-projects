direct_path2 = function(i, j, grid, gridSizeX){
  
  ## Points on the Segment
  
  y2 = grid[j,2]
  y1 = grid[i,2]
  
  x2 = grid[j,1]
  x1 = grid[i,1]
  
  if(y2 == y1){
    return(seq(i, j, by = 1))
  }
  
  if(x2 == x1){
    return(seq(i, j, by = gridSizeX))
  }
  
  ds = sign(x2 - x1)
  
  m = (y2 - y1)/(x2 - x1)
  
  delta_x = 1/sqrt(1 + m^2)
  
  seq = c(seq(x1, x2, by = ds*delta_x), x2)
  
  points = as.matrix(unique(round(cbind(seq, y1 + m*(seq-x1)))))
  
  p = points[,1] + (points[,2]-1)*gridSizeX
  
  return(p)
  
}
