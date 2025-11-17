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
source("functions/PSD_1.R")
source("functions/compute_theta_cosine.R")
source("functions/cosine_firstdev.R")
source("functions/clr.R")
source("functions/plot_grid.R")

set.seed(1234)

cuprite = readMat("Cuprite97_data.mat")

data = 

# to_remove = c(101, 135)
# data = cosine_firstdev(data, to_remove)

# data = clr(data)

gridSizeX = 50
gridSizeY = 50

n_attributes = dim(data)[2]

numNodes = gridSizeX * gridSizeY

grid = expand.grid(x = 1:gridSizeX, y = 1:gridSizeY)

################## Classic PSD_1 #############################

dissimilarity = PSD_1(data, gridSizeX, gridSizeY)

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
     xlab = "Subgroups", main = "Theta values after single PAM - PSD[1]", ylim = c(0,1))
points(ngroups_first, theta, pch = 19)
axis(side = 1, at = ngroups_first, labels = ngroups_first)
grid()

k_values = 2:8
ngroups = c(75, 100, 150)
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
        col = colors, main = "PSD[1] - Theta values", ylim = c(0,1))
matlines(k_values, theta_matrix, lty = 1, lwd = 2, col = colors)
axis(side = 1, labels = k_values, at = k_values)
legend("bottomright", legend = c("n = 75", "n = 100", "n = 150"), 
       col = colors, lty = 1, lwd = 2, pch = 19, cex = 0.6)
grid()

############## Example with k = 4 and ngroups = 75 #################

k = 5
ngroups = 100

cluster = cluster_matrix[as.character(ngroups),]

plot_grid(gridSizeX, gridSizeY, cluster, 
          title = paste0("After single PAM (", ngroups, " subgroups)"))

repr = matrix(0, nrow = ngroups, ncol = dim(data)[2])

for(t in 1:ngroups){
  
  ind = which(cluster == t)
  
  repr[t,] = colMeans(data[ind,])
  
}

finalkm = kmeans(as.data.frame(repr), k, nstart = 100, iter.max = 100)
final_cluster = finalkm$cluster[cluster]

plot_grid(gridSizeX, gridSizeY, final_cluster, 
          title = paste("PSD[1] -", k, "clusters and", ngroups, "subgroups"))

repr = matrix(0, nrow = k, ncol = dim(data)[2])

for(i in 1:k){
  ind = which(final_cluster == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l',ylab = "", lty = 1, col = "lightgrey", main = paste("Centroids when k =", k, "and ngroups =", ngroups))
matlines(t(repr), type = 'l', lty = 1, col = seq(1,k), lwd = 3)
grid()

data = fread("cuprite_50_50_pixels.csv", header = TRUE)

repr = matrix(0, nrow = k, ncol = dim(data)[2])

for(i in 1:k){
  ind = which(cluster == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l', lty = 1, ylab = "", col = "lightgrey", main = paste("Centroids when k =", k, "and n =", ngroups))
matlines(t(repr), type = 'l', lty = 1, col = seq(1,k), lwd = 3)
grid()

################## PSD_1 with cosine diss ##########################

dissimilarity_cos = PSD_1(data, gridSizeX, gridSizeY, method = "cosine")

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
     xlab = "Subgroups", main = "Theta values after single PAM - PSD_1 cos", ylim = c(0,1))
points(ngroups_first, theta, pch = 19)
axis(side = 1, at = ngroups_first, labels = ngroups_first)
grid()

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
        col = colors, main = "PSD_1 cos - Theta values", ylim = c(0,1))
matlines(k_values, theta_matrix, lty = 1, lwd = 2, col = colors)
axis(side = 1, labels = k_values, at = k_values)
grid()
legend("bottomright", legend = c("n = 50", "n = 75", "n = 100"), 
       col = colors, lty = 1, lwd = 2, pch = 19, cex = 0.6)

############## Example with k = 4 and ngroups = 75 #################

k = 4
ngroups = 75

cluster = cluster_matrix[as.character(ngroups),]

plot_grid(gridSizeX, gridSizeY, cluster, 
          title = paste0("After single PAM (", ngroups, " subgroups)"))

repr = matrix(0, nrow = ngroups, ncol = dim(data)[2])

for(t in 1:ngroups){
  
  ind = which(cluster == t)
  
  repr[t,] = colMeans(data[ind,])
  
}

repr_dist = proxy::dist(repr, method = "cosine")

cl = as.numeric(pam(repr_dist, k, diss = TRUE, variant = "faster", cluster.only = TRUE))
final_cluster = cl[cluster]

plot_grid(gridSizeX, gridSizeY, final_cluster, 
          title = paste("PSD_1 cos -", k, "clusters and", ngroups, "subgroups"))

repr = matrix(0, nrow = k, ncol = dim(data)[2])

for(i in 1:k){
  ind = which(final_cluster == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l',ylab = "", lty = 1, col = "lightgrey", main = paste("Centroids when k =", k, "and ngroups =", ngroups))
matlines(t(repr), type = 'l', lty = 1, col = seq(1,k), lwd = 3)
grid()

################## PSD_1 with arcosine diss ########################

dissimilarity_arccos = PSD_1(data, gridSizeX, gridSizeY, method = "arcosine")

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
     xlab = "Subgroups", main = "Theta values after single k-means - PSD[1] arccos", ylim = c(0,1))
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
  
  repr_dist = acos(-proxy::dist(repr, method = "cosine") + 1)
  
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
        col = colors, main = "PSD[1] arccos - Theta values", ylim = c(0,1))
matlines(k_values, theta_matrix, lty = 1, lwd = 2, col = colors)
axis(side = 1, labels = k_values, at = k_values)
legend("bottomright", legend = c("n = 50", "n = 75", "n = 100"), 
       col = colors, lty = 1, lwd = 2, pch = 19, cex = 0.6)

############## Example with k = 4 and ngroups = 75 #################

k = 4
ngroups = 75

cluster = cluster_matrix[as.character(ngroups),]

plot_grid(gridSizeX, gridSizeY, cluster, 
          title = paste0("After single k-means (", ngroups, " subgroups)"))

repr = matrix(0, nrow = ngroups, ncol = dim(data)[2])

for(t in 1:ngroups){
  
  ind = which(cluster == t)
  
  repr[t,] = colMeans(data[ind,])
  
}

repr_dist = acos(-proxy::dist(repr, method = "cosine") + 1)

cl = as.numeric(pam(repr_dist, k, diss = TRUE, variant = "faster", cluster.only = TRUE))
final_cluster = cl[cluster]

plot_grid(gridSizeX, gridSizeY, final_cluster, 
          title = paste("PSD[1] arccos -", k, "clusters and", ngroups, "subgroups"))

repr = matrix(0, nrow = k, ncol = dim(data)[2])

for(i in 1:k){
  ind = which(final_cluster == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l',ylab = "", lty = 1, col = "lightgrey", main = paste("Centroids when k =", k, "and ngroups =", ngroups))
matlines(t(repr), type = 'l', lty = 1, col = seq(1,k), lwd = 3)