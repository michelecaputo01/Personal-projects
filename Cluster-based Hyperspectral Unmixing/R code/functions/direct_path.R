direct_path = function(i, j, grid, gridSizeX){
  
  ## Angular Coefficient 
  
  y1 = grid[j,2]
  y0 = grid[i,2]
  
  x1 = grid[j,1]
  x0 = grid[i,1]
  
  dx = sign(x1 - x0)
  dy = sign(y1 - y0)
  
  if(y0 == y1){
    return(seq(i, j, by = dx))
  }
  
  if(x0 == x1){
    return(seq(i, j, by = dy*gridSizeX))
  }
  
  if(y1 - y0 == x1 - x0){
    return(seq(i, j, by = dy*(gridSizeX + 1)))
  }
  
  if(y1 - y0 == -x1 + x0){
    return(seq(i, j, by = dy*(gridSizeX - 1)))
  }
  
  cells = NULL 
  
  while (TRUE) {
    
    cells = c(cells, x0 + gridSizeX*(y0-1))
    
    if (x0 == x1 && y0 == y1) break
    
    m = (y1 - y0)/(x1 - x0)
    
    if(0 <= m && m <= 0.5){
      
      x0 = x0 + dx
      
    } else if(0.5 < m && m < 2){
      
      x0 = x0 + dx
      y0 = y0 + dy
      
    } else if(-0.5 <= m && m <= 0){
      
      x0 = x0 + dx
      
    } else if(-2 < m && m < -0.5){
      
      x0 = x0 + dx
      y0 = y0 + dy
      
    } else if(abs(m) >= 2){
      
      y0 = y0 + dy
      
    }
    
  }
  
  return(cells)
  
}