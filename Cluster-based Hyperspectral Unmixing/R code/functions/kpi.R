kpi = function(real_labels, fitted){
  
  k = max(fitted)
  
  numNodes = length(fitted)
  
  perm = e1071::permutations(k)
  
  value = 0
  
  for(i in 1:nrow(perm)){
    
    copy_perm = perm[i,][fitted]
    
    new_val = sum(real_labels == copy_perm)/numNodes
    
    if(new_val > value){
      value = new_val
    }
    
  }
  
  return(value)
  
}