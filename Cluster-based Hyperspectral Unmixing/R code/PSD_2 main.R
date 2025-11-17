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
library(R.matlab)

source("functions/direct_path.R")
source("functions/compute_theta.R")
source("functions/PSD_2.R")
source("functions/compute_theta_cosine.R")
source("functions/plot_grid.R")

cuprite = readMat("Cuprite97_data.mat")

data = 

set.seed(1234)

gridSizeX = 50
gridSizeY = 50

n_attributes = dim(data)[2]

numNodes = gridSizeX * gridSizeY

grid = expand.grid(x = 1:gridSizeX, y = 1:gridSizeY)

################### Classic PSD_2 ####################

dissimilarity = PSD_2(data, gridSizeX, gridSizeY)

ngroups_first = c(10, 25, 50, 75, 100, 150, 200, 250, 300)

theta = NULL

cluster_matrix = matrix(0, nrow = length(ngroups_first), ncol = numNodes)
rownames(cluster_matrix) = ngroups_first

for(n in ngroups_first){
  
  cat(sprintf("\rProgress:  ngroups = %d", n))
  flush.console()
  
  cluster = as.numeric(pam(dissimilarity, n, diss = TRUE, variant = "faster", cluster.only = TRUE))
  
  theta = c(theta, compute_theta(data, cluster))
  
  cluster_matrix[as.character(n),] = cluster
  
}

plot(ngroups_first, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Subgroups", main = "Theta values after single k-means - PSD[2]", ylim = c(0,1))
points(ngroups_first, theta, pch = 19)
axis(side = 1, at = ngroups_first, labels = ngroups_first)

k_values = 2:8
ngroups = c(50, 75, 100)
theta_matrix = matrix(0, nrow = length(k_values), ncol = length(ngroups))
rownames(theta_matrix) = k_values
colnames(theta_matrix) = ngroups

for(j in ngroups){
  
  cluster = cluster_matrix[as.character(j),]
  
  repr = matrix(0, nrow = j, ncol = dim(data)[2])
  
  for(t in 1:j){
    
    ind = which(cluster == t)
    
    repr[t,] = colMeans(data[ind,])
    
  }
  
  for(i in k_values){
    
    cat(sprintf("\rProgress: k = %d, ngroups = %d", i, j))
    flush.console()
    
    finalkm = kmeans(as.data.frame(repr), i, nstart = 100, iter.max = 100)
    final_cluster = finalkm$cluster[cluster]
    
    theta_matrix[as.character(i), as.character(j)] = compute_theta(data, final_cluster)
    
  }
  
}

colors = colorRampPalette(c("lightblue", "darkblue"))
colors = colors(length(ngroups))

matplot(k_values, theta_matrix, type = 'p', xaxt = 'n', pch = 19, cex = 1.2, ylab = expression(theta), 
        col = colors, main = "PSD[2] - Theta values", ylim = c(0,1))
matlines(k_values, theta_matrix, lty = 1, lwd = 2, col = colors)
axis(side = 1, labels = k_values, at = k_values)
legend("bottomright", legend = c("n = 50", "n = 75", "n = 100"), 
       col = colors, lty = 1, lwd = 2, pch = 19, cex = 0.5)

############## Example with k = 4 and ngroups = 75 #################

ngroups = 75
k = 4

cluster = cluster_matrix[as.character(ngroups),]

plot_grid(gridSizeX, gridSizeY, cluster, 
          title = paste0("After single k-means (", ngroups, " subgroups)"))

repr = matrix(0, nrow = ngroups, ncol = dim(data)[2])

for(t in 1:ngroups){
  
  ind = which(cluster == t)
  
  repr[t,] = colMeans(data[ind,])
  
}

finalkm = kmeans(as.data.frame(repr), k, nstart = 100, iter.max = 100)
final_cluster = finalkm$cluster[cluster]

plot_grid(gridSizeX, gridSizeY, final_cluster, 
          title = paste("PSD[2] -", k, "clusters and", ngroups, "subgroups"))

repr = matrix(0, nrow = k, ncol = dim(data)[2])

for(i in 1:k){
  ind = which(final_cluster == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l', ylab = "", lty = 1, col = "lightgrey", main = paste("Centroids when k =", k, "and ngroups =", ngroups))
matlines(t(repr), type = 'l', lty = 1, col = seq(1,k), lwd = 3)

################### PSD_2 with cosine diss ####################

dissimilarity_cos = PSD_2(data, gridSizeX, gridSizeY, method = "cosine")

ngroups_first = c(10, 25, 50, 75, 100, 150, 200, 250, 300)

theta = NULL

cluster_matrix = matrix(0, nrow = length(ngroups_first), ncol = numNodes)
rownames(cluster_matrix) = ngroups_first

for(n in ngroups_first){
  
  cat(sprintf("\rProgress:  ngroups = %d", n))
  flush.console()
  
  cluster = as.numeric(pam(dissimilarity_cos, n, diss = TRUE, variant = "faster", cluster.only = TRUE))
  
  theta = c(theta, compute_theta_cosine(data, cluster))
  
  cluster_matrix[as.character(n),] = cluster
  
}

plot(ngroups_first, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Subgroups", main = "Theta values after single k-means - PSD[2] cos", ylim = c(0,1))
points(ngroups_first, theta, pch = 19)
axis(side = 1, at = ngroups_first, labels = ngroups_first)

k_values = 2:8
ngroups = c(50, 75, 100)
theta_matrix = matrix(0, nrow = length(k_values), ncol = length(ngroups))
rownames(theta_matrix) = k_values
colnames(theta_matrix) = ngroups

for(j in ngroups){
  
  cluster = cluster_matrix[as.character(j),]
  
  repr = matrix(0, nrow = j, ncol = dim(data)[2])
  
  for(t in 1:j){
    
    ind = which(cluster == t)
    
    repr[t,] = colMeans(data[ind,])
    
  }
  
  repr_dist = proxy::dist(repr, method = "cosine")
  
  for(i in k_values){
    
    cat(sprintf("\rProgress: k = %d, ngroups = %d", i, j))
    flush.console()
    
    cl = as.numeric(pam(repr_dist, i, diss = TRUE, variant = "faster", cluster.only = TRUE))
    final_cluster = cl[cluster]
    
    theta_matrix[as.character(i), as.character(j)] = compute_theta_cosine(data, final_cluster)
    
  }
  
}

colors = colorRampPalette(c("lightblue", "darkblue"))
colors = colors(length(ngroups))

matplot(k_values, theta_matrix, type = 'p', xaxt = 'n', pch = 19, cex = 1.2, ylab = expression(theta), 
        col = colors, main = "PSD[2] cos - Theta values", ylim = c(0,1))
matlines(k_values, theta_matrix, lty = 1, lwd = 2, col = colors)
axis(side = 1, labels = k_values, at = k_values)
legend("bottomright", legend = c("n = 50", "n = 75", "n = 100"), 
       col = colors, lty = 1, lwd = 2, pch = 19, cex = 0.5)

############## Example with k = 4 and ngroups = 75 #################

ngroups = 75
k = 4

cluster = cluster_matrix[as.character(ngroups),]

plot_grid(gridSizeX, gridSizeY, cluster, 
          title = paste0("After single k-means (", ngroups, " subgroups)"))

repr = matrix(0, nrow = ngroups, ncol = dim(data)[2])

for(t in 1:ngroups){
  
  ind = which(cluster == t)
  
  repr[t,] = colMeans(data[ind,])
  
}

repr_dist = proxy::dist(repr, method = "cosine")

cl = as.numeric(pam(repr_dist, k, diss = TRUE, variant = "faster", cluster.only = TRUE))
final_cluster = cl[cluster]

plot_grid(gridSizeX, gridSizeY, final_cluster, 
          title = paste("PSD[2] cos -", k, "clusters and", ngroups, "subgroups"))

repr = matrix(0, nrow = k, ncol = dim(data)[2])

for(i in 1:k){
  ind = which(final_cluster == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l', ylab = "", lty = 1, col = "lightgrey", main = paste("Centroids when k =", k, "and ngroups =", ngroups))
matlines(t(repr), type = 'l', lty = 1, col = seq(1,k), lwd = 3)

################### PSD_2 with arccosine diss ####################

dissimilarity_arccos = PSD_2(data, gridSizeX, gridSizeY, method = "arcos")

ngroups_first = c(10, 25, 50, 75, 100, 150, 200, 250, 300)

theta = NULL

cluster_matrix = matrix(0, nrow = length(ngroups_first), ncol = numNodes)
rownames(cluster_matrix) = ngroups_first

for(n in ngroups_first){
  
  cat(sprintf("\rProgress:  ngroups = %d", n))
  flush.console()
  
  cluster = as.numeric(pam(dissimilarity_arccos, n, diss = TRUE, variant = "faster", cluster.only = TRUE))
  
  theta = c(theta, compute_theta_cosine(data, cluster))
  
  cluster_matrix[as.character(n),] = cluster
  
}

plot(ngroups_first, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Subgroups", main = "Theta values after single k-means - PSD[2] arccos", ylim = c(0,1))
points(ngroups_first, theta, pch = 19)
axis(side = 1, at = ngroups_first, labels = ngroups_first)

k_values = 2:8
ngroups = c(50, 75, 100)
theta_matrix = matrix(0, nrow = length(k_values), ncol = length(ngroups))
rownames(theta_matrix) = k_values
colnames(theta_matrix) = ngroups

for(j in ngroups){
  
  cluster = cluster_matrix[as.character(j),]
  
  repr = matrix(0, nrow = j, ncol = dim(data)[2])
  
  for(t in 1:j){
    
    ind = which(cluster == t)
    
    repr[t,] = colMeans(data[ind,])
    
  }
  
  repr_dist = proxy::dist(repr, method = "cosine")
  
  for(i in k_values){
    
    cat(sprintf("\rProgress: k = %d, ngroups = %d", i, j))
    flush.console()
    
    cl = as.numeric(pam(repr_dist, i, diss = TRUE, variant = "faster", cluster.only = TRUE))
    final_cluster = cl[cluster]
    
    theta_matrix[as.character(i), as.character(j)] = compute_theta_cosine(data, final_cluster)
    
  }
  
}

colors = colorRampPalette(c("lightblue", "darkblue"))
colors = colors(length(ngroups))

matplot(k_values, theta_matrix, type = 'p', xaxt = 'n', pch = 19, cex = 1.2, ylab = expression(theta), 
        col = colors, main = "PSD[2] arccos - Theta values", ylim = c(0,1))
matlines(k_values, theta_matrix, lty = 1, lwd = 2, col = colors)
axis(side = 1, labels = k_values, at = k_values)
legend("bottomright", legend = c("n = 50", "n = 75", "n = 100"), 
       col = colors, lty = 1, lwd = 2, pch = 19, cex = 0.5)

############## Example with k = 4 and ngroups = 75 #################

ngroups = 75
k = 4

cluster = cluster_matrix[as.character(ngroups),]

plot_grid(gridSizeX, gridSizeY, cluster, 
          title = paste0("After single k-means (", ngroups, " subgroups)"))

repr = matrix(0, nrow = ngroups, ncol = dim(data)[2])

for(t in 1:ngroups){
  
  ind = which(cluster == t)
  
  repr[t,] = colMeans(data[ind,])
  
}

repr_dist = proxy::dist(repr, method = "cosine")

cl = as.numeric(pam(repr_dist, k, diss = TRUE, variant = "faster", cluster.only = TRUE))
final_cluster = cl[cluster]

plot_grid(gridSizeX, gridSizeY, final_cluster, 
          title = paste("PSD[2] arccos -", k, "clusters and", ngroups, "subgroups"))

repr = matrix(0, nrow = k, ncol = dim(data)[2])

for(i in 1:k){
  ind = which(final_cluster == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l', ylab = "", lty = 1, col = "lightgrey", main = paste("Centroids when k =", k, "and ngroups =", ngroups))
matlines(t(repr), type = 'l', lty = 1, col = seq(1,k), lwd = 3)