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

source("functions/compute_theta.R")
source("functions/compute_theta_cosine.R")
source("functions/plot_grid.R")

set.seed(1234)

cuprite = readMat("Cuprite97_data.mat")

data = 

gridSizeX = 50
gridSizeY = 50

n_attributes = dim(data)[2]

numNodes = gridSizeX * gridSizeY

grid = expand.grid(x = 1:gridSizeX, y = 1:gridSizeY)

###################### Raw kmeans ##############

ngroups = seq(2,15)

theta = NULL

cluster_matrix = matrix(0, nrow = length(ngroups), ncol = numNodes)
rownames(cluster_matrix) = ngroups

for(n in ngroups){
  
  cat(sprintf("\rProgress:  ngroups = %d", n))
  flush.console()
  
  cluster = kmeans(data, centers = n, nstart = 50, iter.max = 100)$cluster
  
  theta = c(theta, compute_theta(data, cluster))
  
  cluster_matrix[as.character(n),] = cluster
  
}

plot(ngroups, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Subgroups", main = "Theta values after classic k-means", ylim = c(0,1))
points(ngroups, theta, pch = 19)
axis(side = 1, at = ngroups, labels = ngroups)
grid()

n = 7

plot_grid(gridSizeX, gridSizeY, cluster_matrix[as.character(n),], 
          title = paste0("After classic k-means (", n, " clusters)"))

repr = matrix(0, nrow = n, ncol = n_attributes)

for(i in 1:n){
  ind = which(cluster_matrix[as.character(n),] == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l',ylab = "", lty = 1, col = "lightgrey", main = paste("Centroids when k =", n))
matlines(t(repr), type = 'l', lty = 1, col = seq(1,n), lwd = 3)
grid()

###################### Raw PAM with cosine diss ##############

ngroups = seq(2,15)

theta = NULL

cluster_matrix = matrix(0, nrow = length(ngroups), ncol = numNodes)
rownames(cluster_matrix) = ngroups

diss = proxy::dist(data, method = "cosine")

for(n in ngroups){
  
  cat(sprintf("\rProgress:  ngroups = %d", n))
  flush.console()
  
  cluster = as.numeric(pam(diss, n, diss = TRUE, variant = "faster", cluster.only = TRUE))
  
  theta = c(theta, compute_theta_cosine(data, cluster))
  
  cluster_matrix[as.character(n),] = cluster
  
}

plot(ngroups, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Subgroups", main = "Theta values after classic PAM with cosine dissimilarity", ylim = c(0,1))
points(ngroups, theta, pch = 19)
axis(side = 1, at = ngroups, labels = ngroups)
grid()

n = 4

plot_grid(gridSizeX, gridSizeY, cluster_matrix[as.character(n),], 
          title = paste0("After classic PAM with cosine dissimilarity (", n, " clusters)"))

repr = matrix(0, nrow = n, ncol = n_attributes)

for(i in 1:n){
  ind = which(cluster_matrix[as.character(n),] == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l',ylab = "", lty = 1, col = "lightgrey", main = paste("Centroids when k =", n))
matlines(t(repr), type = 'l', lty = 1, col = seq(1,n), lwd = 3)
grid()

###################### Raw PAM with arcosine diss ##############

ngroups = seq(2,15)

theta = NULL

cluster_matrix = matrix(0, nrow = length(ngroups), ncol = numNodes)
rownames(cluster_matrix) = ngroups

diss_arcos = acos(-proxy::dist(data, method = "cosine") + 1)

for(n in ngroups){
  
  cat(sprintf("\rProgress:  ngroups = %d", n))
  flush.console()
  
  cluster = as.numeric(pam(diss_arcos, n, diss = TRUE, variant = "faster", cluster.only = TRUE))
  
  theta = c(theta, compute_theta_cosine(data, cluster))
  
  cluster_matrix[as.character(n),] = cluster
  
}

plot(ngroups, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Subgroups", main = "Theta values after classic PAM with arcosine dissimilarity", ylim = c(0,1))
points(ngroups, theta, pch = 19)
axis(side = 1, at = ngroups, labels = ngroups)
grid()

n = 6

plot_grid(gridSizeX, gridSizeY, cluster_matrix[as.character(n),], 
          title = paste0("After classic PAM with arcosine dissimilarity (", n, " clusters)"))

repr = matrix(0, nrow = n, ncol = n_attributes)

for(i in 1:n){
  ind = which(cluster_matrix[as.character(n),] == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l',ylab = "", lty = 1, col = "lightgrey", main = paste("Centroids when k =", n))
matlines(t(repr), type = 'l', lty = 1, col = seq(1,n), lwd = 3)
grid()

###################### Hierarchical clustering L2 ###########

dist = dist(data, method = 'euclidean')

hc_single = hclust(dist, method = "single")
hc_average = hclust(dist, method = "average")
hc_complete = hclust(dist, method = "complete")

plot(hc_single, main = "Dendogram (HC single, L2 dist)", xlab = '', labels = FALSE, sub = "")
plot(hc_average, main = "Dendogram (HC average, L2 dist)", xlab = '', labels = FALSE, sub = "")
plot(hc_complete, main = "Dendogram (HC complete, L2 dist)", xlab = '', labels = FALSE, sub = "")

ngroups = seq(2,10)

theta = NULL

cluster_matrix = matrix(0, nrow = length(ngroups), ncol = numNodes)
rownames(cluster_matrix) = ngroups

for(n in ngroups){
  
  cat(sprintf("\rProgress:  ngroups = %d", n))
  flush.console()
  
  cluster = cutree(hc_complete, k = n)
  
  theta = c(theta, compute_theta(data, cluster))
  
  cluster_matrix[as.character(n),] = cluster
  
}

plot(ngroups, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Subgroups", main = "Theta values after HC complete, L2 dist", ylim = c(0,1))
points(ngroups, theta, pch = 19)
axis(side = 1, at = ngroups, labels = ngroups)
grid()

n = 5

plot_grid(gridSizeX, gridSizeY, cluster_matrix[as.character(n),], 
          title = paste0("After HC complete, L2 dist (", n, " clusters)"))

repr = matrix(0, nrow = n, ncol = n_attributes)

for(i in 1:n){
  ind = which(cluster_matrix[as.character(n),] == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l',ylab = "", lty = 1, col = "lightgrey", main = paste("Centroids when k =", n))
matlines(t(repr), type = 'l', lty = 1, col = seq(1,n), lwd = 3)
grid()

###################### Hierarchical clustering cosine ###########

dist = proxy::dist(data, method = 'cosine')

hc_single = hclust(dist, method = "single")
hc_average = hclust(dist, method = "average")
hc_complete = hclust(dist, method = "complete")

plot(hc_single, main = "Dendogram (HC single, cosine dist)", xlab = '', labels = FALSE, sub = "")
plot(hc_average, main = "Dendogram (HC average, cosine dist)", xlab = '', labels = FALSE, sub = "")
plot(hc_complete, main = "Dendogram (HC complete, cosine dist)", xlab = '', labels = FALSE, sub = "")

ngroups = seq(2,10)

theta = NULL

cluster_matrix = matrix(0, nrow = length(ngroups), ncol = numNodes)
rownames(cluster_matrix) = ngroups

for(n in ngroups){
  
  cat(sprintf("\rProgress:  ngroups = %d", n))
  flush.console()
  
  cluster = cutree(hc_complete, k = n)
  
  theta = c(theta, compute_theta_cosine(data, cluster))
  
  cluster_matrix[as.character(n),] = cluster
  
}

plot(ngroups, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Subgroups", main = "Theta values after HC complete, cosine dist", ylim = c(0,1))
points(ngroups, theta, pch = 19)
axis(side = 1, at = ngroups, labels = ngroups)
grid()

n = 5

plot_grid(gridSizeX, gridSizeY, cluster_matrix[as.character(n),], 
          title = paste0("After HC complete, cosine dist (", n, " clusters)"))

repr = matrix(0, nrow = n, ncol = n_attributes)

for(i in 1:n){
  ind = which(cluster_matrix[as.character(n),] == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l',ylab = "", lty = 1, col = "lightgrey", main = paste("Centroids when k =", n))
matlines(t(repr), type = 'l', lty = 1, col = seq(1,n), lwd = 3)
grid()

###################### Hierarchical clustering arcosine ####

dist = acos(-proxy::dist(data, method = 'cosine') + 1)

hc_single = hclust(dist, method = "single")
hc_average = hclust(dist, method = "average")
hc_complete = hclust(dist, method = "complete")

plot(hc_single, main = "Dendogram (HC single, arcosine dist)", xlab = '', labels = FALSE, sub = "")
plot(hc_average, main = "Dendogram (HC average, arcosine dist)", xlab = '', labels = FALSE, sub = "")
plot(hc_complete, main = "Dendogram (HC complete, arcosine dist)", xlab = '', labels = FALSE, sub = "")

ngroups = seq(2,10)

theta = NULL

cluster_matrix = matrix(0, nrow = length(ngroups), ncol = numNodes)
rownames(cluster_matrix) = ngroups

for(n in ngroups){
  
  cat(sprintf("\rProgress:  ngroups = %d", n))
  flush.console()
  
  cluster = cutree(hc_complete, k = n)
  
  theta = c(theta, compute_theta_cosine(data, cluster))
  
  cluster_matrix[as.character(n),] = cluster
  
}

plot(ngroups, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Subgroups", main = "Theta values after HC complete, arccosine dist", ylim = c(0,1))
points(ngroups, theta, pch = 19)
axis(side = 1, at = ngroups, labels = ngroups)
grid()

n = 5

plot_grid(gridSizeX, gridSizeY, cluster_matrix[as.character(n),], 
          title = paste0("After HC complete, arcosine dist (", n, " clusters)"))

repr = matrix(0, nrow = n, ncol = n_attributes)

for(i in 1:n){
  ind = which(cluster_matrix[as.character(n),] == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l',ylab = "", lty = 1, col = "lightgrey", main = paste("Centroids when k =", n))
matlines(t(repr), type = 'l', lty = 1, col = seq(1,n), lwd = 3)
grid()