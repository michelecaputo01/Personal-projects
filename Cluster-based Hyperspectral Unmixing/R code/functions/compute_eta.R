compute_eta = function(cluster_matrix){
  
  k = max(cluster_matrix)
  
  numNodes = dim(cluster_matrix)[2]
  
  B = dim(cluster_matrix)[1]
  
  eta_matrix = apply(cluster_matrix, 2, FUN = function(i) table(factor(i, levels = seq(1,k))) / B)
  
  eta_x = -eta_matrix*log(eta_matrix)
  eta_x[is.na(eta_x)] = 0
  eta_x = apply(eta_x, 2, sum)/log(k)
  
  eta = mean(eta_x)
  
  results = list( eta_x = eta_x,
                  eta = eta )
  
  return(results)
  
}