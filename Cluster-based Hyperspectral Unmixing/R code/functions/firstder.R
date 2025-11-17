firstder = function(data, to_remove = NULL){
  
  # to_remove = point indexes on the right of the discontinuities
  
  natt = ncol(data)
  
  d = as.matrix(data)
  
  diff = d[,-1] - d[,-natt]
  
  if(!is.null(to_remove)){
    diff = diff[,-c(to_remove-1)]
  }
  
  return(as.data.frame(diff))
  
}
