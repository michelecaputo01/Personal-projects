library(igraph)
library(ggplot2)
library(gridExtra)
library(akima)
library(fields)
library(MASS)
library(cluster)
library(ggplot2)
library(progress)
library(lattice)
library(progress)
library(mgcv)
library(data.table)
library(rgl)
library(psych)
library(R.matlab)
library(Polychrome)

source("functions/compute_theta.R")
source("functions/compute_theta_cosine.R")
source("functions/PSD_1.R")
source("functions/PSD_1_avg.R")
source("functions/direct_path.R")
source("functions/plot_grid.R")
source("functions/sum_distance_matrix.R")

set.seed(1234)

cuprite = readMat("Cuprite97_data.mat")

data_disorder = as.data.frame(t(cuprite$x))

gridSizeX = 190
gridSizeY = 250

indexes = NULL

original_grid = matrix(1:(gridSizeX * gridSizeY), nrow = gridSizeY, 
                       ncol = gridSizeX, byrow = FALSE)

for(i in gridSizeY:1){
  for(j in 1:gridSizeX){
    indexes = c(indexes, original_grid[i,j])
  }
}

original_data = data_disorder[indexes,]

data = original_data

n_attributes = dim(data)[2]

numNodes = gridSizeX * gridSizeY

grid = expand.grid(x = 1:gridSizeX, y = 1:gridSizeY)

########################## L2 contribution #######

dissimilarity = sum_distance_matrix(data, grid, "l2")

k_values = seq(2,20)
theta = NULL

for(k in k_values){
  
  cat(sprintf("\rProgress:  k = %d", k))
  flush.console()
  
  cluster = as.numeric(pam(dissimilarity, k, diss = TRUE, variant = "faster", cluster.only = TRUE))
  
  theta = c(theta, compute_theta(data, cluster))
  
}

plot(k_values, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Number of cluster", main = "Theta values after PAM - Mixed distance (L2 contr)", ylim = c(0,1))
points(k_values, theta, pch = 19)
axis(side = 1, at = k_values, labels = k_values)
grid()

k = 9
cl = as.numeric(pam(dissimilarity, k, diss = TRUE, variant = "faster", cluster.only = TRUE))

plot_grid(gridSizeX, gridSizeY, cl, title = "Example with PAM and 10 clusters (L2 contr)")

theta = NULL

cluster_matrix = matrix(0, nrow = k-1, ncol = numNodes)
rownames(cluster_matrix) = seq(2,k)

repr = matrix(0, nrow = k, ncol = dim(data)[2])

for(i in 1:k){
  ind = which(cl == i)
  repr[i,] = colMeans(data[ind,])
}

for(nclus in 2:(k-1)){
  
  cat(sprintf("\rProgress:  %d", nclus))
  flush.console()
  
  finalkm = kmeans(as.data.frame(repr), nclus, nstart = 20)
  final_cluster = finalkm$cluster[cl]
  
  cluster_matrix[as.character(nclus), ] = final_cluster
  
  theta = c(theta, compute_theta(data, final_cluster))
  
}

theta = c(theta, compute_theta(data, cl))

cluster_matrix[as.character(k), ] = cl

plot(2:k, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Number of cluster", main = "Theta values after second PAM - Mixed distance (L2 contr)", ylim = c(0,1))
points(2:k, theta, pch = 19)
axis(side = 1, at = 2:k, labels = 2:k)
grid()

k = 5

plot_grid(gridSizeX, gridSizeY, cluster_matrix[as.character(k),], title = paste("Example with", k, "clusters (L2 contr)"))

repr = matrix(0, nrow = k, ncol = dim(data)[2])

for(i in 1:k){
  ind = which(cluster_matrix[as.character(k), ] == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l', ylab = "", lty = 1, col = "lightgrey",
        main = paste("Centroids with", k, "clusters"))
matlines(t(repr), type = 'l', lty = 1, col = seq(1,k), lwd = 3)
grid()

data = fread("cuprite_50_50_pixels.csv", header = TRUE)

repr = matrix(0, nrow = k, ncol = dim(data)[2])

for(i in 1:k){
  ind = which(cluster_matrix[as.character(k), ] == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l', lty = 1, ylab = "", col = "lightgrey", 
        main = paste("Centroids with", k, "clusters"))
matlines(t(repr), type = 'l', lty = 1, col = seq(1,k), lwd = 3)
grid()

########################## Cosine contribution #######

dissimilarity = sum_distance_matrix(data, grid, "cosine")

k_values = seq(2,20)
theta = NULL

for(k in k_values){
  
  cat(sprintf("\rProgress:  k = %d", k))
  flush.console()
  
  cluster = as.numeric(pam(dissimilarity, k, diss = TRUE, variant = "faster", cluster.only = TRUE))
  
  theta = c(theta, compute_theta_cosine(data, cluster))
  
}

plot(k_values, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Number of cluster", main = "Theta values after PAM - Mixed distance (Cosine contr)", ylim = c(0,1))
points(k_values, theta, pch = 19)
axis(side = 1, at = k_values, labels = k_values)
grid()

k = 10
cl = as.numeric(pam(dissimilarity, k, diss = TRUE, variant = "faster", cluster.only = TRUE))

plot_grid(gridSizeX, gridSizeY, cl, title = "Example with PAM and 9 clusters (Cosine contr)")

theta = NULL

cluster_matrix = matrix(0, nrow = k-1, ncol = numNodes)
rownames(cluster_matrix) = seq(2,k)

repr = matrix(0, nrow = k, ncol = dim(data)[2])

for(i in 1:k){
  ind = which(cl == i)
  repr[i,] = colMeans(data[ind,])
}

repr_diss = proxy::dist(as.data.frame(repr), method = "cosine")

for(nclus in 2:(k-1)){
  
  cat(sprintf("\rProgress:  %d", nclus))
  flush.console()
  
  finalcl = as.numeric(pam(repr_diss, nclus, diss = TRUE, variant = "faster", cluster.only = TRUE))
  
  final_cluster = finalcl[cl]
  
  cluster_matrix[as.character(nclus), ] = final_cluster
  
  theta = c(theta, compute_theta_cosine(data, final_cluster))
  
}

theta = c(theta, compute_theta_cosine(data, cl))

cluster_matrix[as.character(k), ] = cl

plot(2:k, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Number of cluster", main = "Theta values after second PAM - Mixed distance (Cosine contr)", ylim = c(0,1))
points(2:k, theta, pch = 19)
axis(side = 1, at = 2:k, labels = 2:k)
grid()

k = 3

plot_grid(gridSizeX, gridSizeY, cluster_matrix[as.character(k),], title = "Example with 3 clusters (Cosine contr)")

repr = matrix(0, nrow = k, ncol = dim(data)[2])

for(i in 1:k){
  ind = which(cluster_matrix[as.character(k), ] == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l', ylab = "", lty = 1, col = "lightgrey",
        main = "Centroids with 4 clusters")
matlines(t(repr), type = 'l', lty = 1, col = seq(1,k), lwd = 3)
grid()