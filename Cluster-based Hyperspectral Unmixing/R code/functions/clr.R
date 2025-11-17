clr = function(data){
  
  g_mean = apply(data, 1, function(row) geometric.mean(row))
  
  g_mat = matrix(g_mean, nrow = nrow(data), ncol = ncol(data), byrow = FALSE)
  
  return(as.data.frame(log(data/g_mat)))
  
}