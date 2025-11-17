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
library(gtools)
library(psych)
library(colorspace)

source("functions/compute_theta.R")
source("functions/compute_theta_cosine.R")
source("functions/PSD_1.R")
source("functions/PSD_1_avg.R")
source("functions/direct_path.R")
source("functions/generate_data.R")
source("functions/kpi.R")
source("functions/plot_grid.R")
source("functions/sum_distance_matrix.R")
source("functions/clr.R")
source("functions/firstder.R")

set.seed(1234)

gridSizeX = 50
gridSizeY = 50

numNodes = gridSizeX * gridSizeY

natt = 187

grid = expand.grid(x = 1:gridSizeX, y = 1:gridSizeY)

minshift = -10
maxshift = 10

sigma_f = 3      ## amplitude
scale = 0.3     ## bigger --> smooth

d = generate_data(gridSizeX, gridSizeY, natt, sigma_f, scale, minshift, maxshift)
original_data = d$data
real_labels = d$real_labels

clr_data = clr(original_data)

first_der_data = firstder(original_data)

data = clr_data

########################## PSD_1 contribution #######

s_mean = colMeans(grid)
sigma2_s = mean(rowSums((as.matrix(grid) - matrix(s_mean, nrow(grid), ncol(grid), byrow = TRUE))^2))
space_contr = as.matrix(dist(grid, method = "euclidean", diag = TRUE, upper = TRUE))^2

psd_matrix = as.matrix(PSD_1(data, gridSizeX, gridSizeY, method = "l2"))

min_dist_point = as.numeric(which.min(rowSums(psd_matrix^2)))
sigma2_f = mean(psd_matrix[,min_dist_point]^2)
funct_contr = psd_matrix^2

dissimilarity = sqrt(space_contr/sigma2_s + funct_contr/sigma2_f)

plot(space_contr[-1,1]/sigma2_s, dissimilarity[-1,1]^2, xlab = "Normalized spatial distance", 
     ylab = "Our distance", main = "Comparison w.r.t. the first observation (PSD_1 contr)",
     pch = 19, col = "grey", lwd = 1, asp = 1)
grid()
lines(c(-10,50), c(-10,50), lwd = 3, col = "black", lty = 2)

k_values = seq(2,20)
theta = NULL

for(k in k_values){
  
  cat(sprintf("\rProgress:  k = %d", k))
  flush.console()
  
  cluster = as.numeric(pam(dissimilarity, k, diss = TRUE, variant = "faster", cluster.only = TRUE))
  
  theta = c(theta, compute_theta(data, cluster))
  
}

plot(k_values, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Number of cluster", main = "Theta values after PAM - Mixed distance (PSD_1 contr)", ylim = c(0,1))
points(k_values, theta, pch = 19)
grid()
axis(side = 1, at = k_values, labels = k_values)

plot_grid(gridSizeX, gridSizeY, as.numeric(pam(dissimilarity, 13, diss = TRUE, variant = "faster", cluster.only = TRUE)),
          title = main = "Example with PAM and 13 clusters (PSD_1 contr)", scaled = 0, centroids = FALSE, original = FALSE)

########################## Cosine contribution #########

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

k = 6
cl = as.numeric(pam(dissimilarity, k, diss = TRUE, variant = "faster", cluster.only = TRUE))

plot_grid(gridSizeX, gridSizeY, cl, title = "Example with PAM and 6 clusters (Cosine contr)",
          scaled = 0, centroid = FALSE, original = FALSE)

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
  
  repr_diss = proxy::dist(as.data.frame(repr), method = "cosine")
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

kpi_value = kpi(real_labels, cluster_matrix[as.character(k), ])

plot_grid(gridSizeX, gridSizeY, cluster_matrix[as.character(k), ], 
          title = paste("FD + Sum distance (Cosine), hard - KPI =", kpi_value), 
          scaled = 0, centroids = TRUE, data = data, title_centr = "Centroids with 3 clusters", 
          original = TRUE, original_data = original_data)

########################## L2 contribution ###########

dissimilarity = sum_distance_matrix(data, grid, "l2")

k_values = seq(5,50, by = 5)
theta = NULL

for(k in k_values){
  
  cat(sprintf("\rProgress:  k = %d", k))
  flush.console()
  
  cluster = as.numeric(pam(dissimilarity, k, diss = TRUE, variant = "faster", cluster.only = TRUE))
  
  theta = c(theta, compute_theta(data, cluster))
  
}

plot(k_values, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Number of clusters", main = "Theta values after PAM - Mixed distance (L2 contr)", ylim = c(0,1))
points(k_values, theta, pch = 19)
axis(side = 1, at = k_values, labels = k_values)
grid()

k = 30
cl = as.numeric(pam(dissimilarity, k, diss = TRUE, variant = "faster", cluster.only = TRUE))

plot_grid(gridSizeX, gridSizeY, cl, title = "Example with PAM and 30 clusters (L2 contr)", 
          scaled = 0, centroid = FALSE, original = FALSE)

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

k = 8

kpi_value = kpi(real_labels, cluster_matrix[as.character(k), ])

plot_grid(gridSizeX, gridSizeY, cluster_matrix[as.character(k),], 
          title = paste("CLR + Sum distance (L2), hard - KPI =", kpi_value), 
          scaled = 0, centroids = TRUE, data = data, title_centr = "Centroids with 8 clusters", 
          original = TRUE, original_data = original_data)

########################## L2 contribution (only functional) #######

dissimilarity = as.matrix(dist(data, method = "euclidean", diag = TRUE, upper = TRUE))

k_values = seq(2,20)
theta = NULL

for(k in k_values){
  
  cat(sprintf("\rProgress:  k = %d", k))
  flush.console()
  
  cluster = as.numeric(pam(dissimilarity, k, diss = TRUE, variant = "faster", cluster.only = TRUE))
  
  theta = c(theta, compute_theta(data, cluster))
  
}

plot(k_values, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Number of cluster", main = "Theta values after PAM - Only functional (L2) distance", ylim = c(0,1))
points(k_values, theta, pch = 19)
axis(side = 1, at = k_values, labels = k_values)
grid()

k = 4
cl = as.numeric(pam(dissimilarity, k, diss = TRUE, variant = "faster", cluster.only = TRUE))

kpi_value = kpi(real_labels, cl)

plot_grid(gridSizeX, gridSizeY, cl, title = paste("Only functional distance (L2), hard - KPI =", kpi_value), 
          scaled = 0, centroids = TRUE, data = data, title_centr = "Centroids with 4 clusters", 
          original = TRUE, original_data = original_data)