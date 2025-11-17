direct_path1 = function(i, j, grid, gridSizeX){
  
  ## Bresenham Algo
  
  cells = NULL
  
  y1 = grid[j,2]
  y0 = grid[i,2]
  
  x1 = grid[j,1]
  x0 = grid[i,1]
  
  dx = abs(x1 - x0)
  dy = abs(y1 - y0)
  
  sx = ifelse(x0 < x1, 1, -1)
  sy = ifelse(y0 < y1, 1, -1)
  
  err = dx - dy
  
  while (TRUE) {
    
    cells = c(cells, x0 + gridSizeX*(y0-1))
    
    if (x0 == x1 && y0 == y1) break
    
    e2 = 2 * err
    
    if (e2 > -dy) {
      err = err - dy
      x0 = x0 + sx
    }
    
    if (e2 < dx) {
      err = err + dx
      y0 = y0 + sy
    }
    
  }
  
  return(cells)
  
}