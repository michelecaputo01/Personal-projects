library(permute)
library(deldir)
library(Matrix)
library(MatchIt)
library(fpc)
library(data.table)
library(cluster)
library(proxy)
library(e1071)
library(progress)
library(ggplot2)
library(gridExtra)
library(e1071)
library(mvtnorm)
library(matrixStats)
library(clue)
library(label.switching)
library(ClusterR)
library(microbenchmark)

set.seed(1234)

source("functions/compute_local_representative.R")
source("functions/BVKMA.R")
source("functions/boot_iter.R")
source("functions/compute_eta.R")
source("functions/compute_theta.R")
source("functions/generate_data.R")

len = seq(30, 150, by = 10)

natt = 187

minshift = 0.75
maxshift = 1.25

sigma_f = 1
scale = 1

time = NULL

n = 50
k = 4

d = generate_data(min(len), min(len), natt, sigma_f, scale, minshift, maxshift, print = FALSE)$data

for(l in len){
  
  cat(sprintf("\rProgress: %d/%d", l, max(len)))
  flush.console()
  
  data = as.data.frame(d[sample(1:nrow(d), size = l^2, replace = TRUE),])
  grid_points = expand.grid(x = 1:l, y = 1:l)
  
  t = microbenchmark(boot_iter(data, n, k, grid_points), times = 30)$time
  
  time = c(time, mean(t)/1e6)
  
}

seq_of_points = seq(min(len)^2, max(len)^2, by = 5)
costant = 70
dim = seq(0, max(len^2))

plot(len^2, time, xlab = "Number of points", ylab = "Time", type = 'l', col = "black",
     main = "Time elapsed (ms) for an iteration of Bootstrap (n and k fixed)", lwd = 3)
points(len^2, time, pch = 19, lwd = 3)
lines(dim, dim/costant, lty = 2, lwd = 2, col = "forestgreen")
lines(dim, (dim/costant)^2, lty = 2, lwd = 2, col = "blue")
lines(dim, log2(dim/costant)*(dim/costant), lty = 2, lwd = 2, col = "violet")
grid()
legend("bottomright", legend = c("Real complexity", "O(N)", expression(O(N*log[2](N))), expression(O(N^2))),
       col = c("black", "forestgreen","violet", "blue"), lty = c(1,2,2,2), lwd = c(3,2,2,2), 
       pch = c(19, NA, NA, NA), cex = 0.8)

#####################################################

len = 50

natt = 187

noise = 0.3
minshift = 0.75
maxshift = 1.25

data = generate_data(len, len, natt, sigma_f, scale, minshift, maxshift, print = FALSE)$data

B = 50
k = 4

grid_points = expand.grid(x = 1:len, y = 1:len)

n_val = seq(50, 500, by = 50)

time = NULL

for(n in n_val){
  
  cat(sprintf("\rProgress: %d/%d", n, max(n_val)))
  flush.console()
  
  t = microbenchmark(boot_iter(data, n, k, grid_points), times = 30)$time
  
  time = c(time, mean(t)/1e6)
  
}

costant = 1.1
dim = seq(0, max(n_val))

plot(n_val, time, xlab = "n values", ylab = "Time", type = 'l', col = "black",
     main = "Time elapsed (ms) for an iteration of Bootstrap (data and k fixed)", lwd = 3)
points(n_val, time, pch = 19, lwd = 3)
lines(dim, dim/costant, lty = 2, lwd = 2, col = "forestgreen")
lines(dim, log2(dim/costant)*(dim/costant), lty = 2, lwd = 2, col = "violet")
grid()
legend("bottomright", legend = c("Real complexity", "O(n)", expression(O(n*log[2](n)))),
       col = c("black", "forestgreen","violet"), lty = c(1,2,2), lwd = c(3,2,2), 
       pch = c(19, NA, NA), cex = 0.8)


## Complexity of ECR1 (label switching) :
## - Hungarian Algo (used for single permutation) : O(K^3)
## - N permutations 
## So the total complexity of ECR1 is O(N * K^3)

