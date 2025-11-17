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
source("functions/PSD_1_thresh.R")
source("functions/plot_grid.R")

cuprite = readMat("Cuprite97_data.mat")

data =

set.seed(1234)

gridSizeX = 50
gridSizeY = 50

n_attributes = dim(data)[2]

numNodes = gridSizeX * gridSizeY

grid = expand.grid(x = 1:gridSizeX, y = 1:gridSizeY)

diss_with_threshold = list()

thresh_values = c(10, 15, 20, 25)

circle = function(center_x, center_y, r, col, lwd) {
  theta = seq(0, 2 * pi, length.out = 500)
  x = center_x + r * cos(theta)
  y = center_y + r * sin(theta)
  lines(x, y, col = col, lwd = lwd)
}

plot(grid, pch = 19, asp = 1, xlim = c(1, gridSizeX), ylim = c(1,gridSizeY))
for(t in thresh_values){
  circle(10, 10, t, 'red', 2)
}

for(t in thresh_values){
  
  diss_with_threshold[[as.character(t)]] = PSD_1_thresh(data, gridSizeX, gridSizeY, t)
  
}

ngroups = c(50, 75, 100, 125, 150, 200)
thresh_ngroups_matrix = matrix(0, nrow = length(thresh_values), ncol = length(ngroups))
rownames(thresh_ngroups_matrix) = thresh_values
colnames(thresh_ngroups_matrix) = ngroups

for(t in thresh_values){
  
  for(j in ngroups){
    
    cat(sprintf("\rProgress: Threshold = %f, ngroups = %d       ", t, j))
    flush.console()
    
    cluster = as.numeric(pam(diss_with_threshold[[as.character(t)]], j, diss = TRUE, variant = "faster", cluster.only = TRUE))
    
    thresh_ngroups_matrix[as.character(t), as.character(j)] = compute_theta(data, cluster)
    
  }
  
}

colors = colorRampPalette(c("lightblue", "darkblue"))
colors = colors(length(thresh_values))

matplot(ngroups, t(thresh_ngroups_matrix), type = 'p', ylab = '', col = colors,
        pch = 19, ylim = c(0,1), main = "Theta as function of ngroups with different thresholds")
matlines(ngroups, t(thresh_ngroups_matrix), lty = 1, lwd = 2, col = colors)
axis(side = 1, labels = ngroups, at = ngroups)
legend("bottomright", legend = c("t = 10", "t = 15", "t = 20", "t = 25"), 
       col = colors, lty = 1, lwd = 2, pch = 19, cex = 0.6)

k_values = 2:8
ngroups = c(100, 125, 150)
theta_matrix_list = rep(list(matrix(0, nrow = length(k_values), ncol = length(ngroups))), length(thresh_values))
names(theta_matrix_list) = thresh_values

for(t in thresh_values){
  
  rownames(theta_matrix_list[[as.character(t)]]) = k_values
  colnames(theta_matrix_list[[as.character(t)]]) = ngroups
  
  for(j in ngroups){
    
    cluster = as.numeric(pam(diss_with_threshold[[as.character(t)]], j, diss = TRUE, variant = "faster", cluster.only = TRUE))
    
    repr = matrix(0, nrow = j, ncol = dim(data)[2])
    
    for(r in 1:j){
      
      ind = which(cluster == r)
      
      repr[r,] = colMeans(data[ind,])
      
    }
    
    for(i in k_values){
      
      cat(sprintf("\rProgress: Threshold = %f, k = %d, ngroups = %d       ", t, i, j))
      flush.console()
      
      finalkm = kmeans(as.data.frame(repr), i, nstart = 100, iter.max = 100)
      final_cluster = finalkm$cluster[cluster]
      
      theta_matrix_list[[as.character(t)]][as.character(i), as.character(j)] = compute_theta(data, final_cluster)
      
    }
    
  }
  
}

colors = colorRampPalette(c("lightblue", "darkblue"))
colors = colors(length(ngroups))

for(t in thresh_values){
  
  matplot(k_values, theta_matrix_list[[as.character(t)]], type = 'p', lty = 1, xaxt = 'n', ylab = expression(theta),
          pch = 19, main = paste("Theta for threshold =", t), ylim = c(0,1), col = colors)
  matlines(k_values, theta_matrix_list[[as.character(t)]], lty = 1, lwd = 2, col = colors)
  axis(side = 1, labels = k_values, at = k_values)
  legend("bottomleft", legend = c("n = 100", "n = 125", "n = 150"), 
         col = colors, lty = 1, lwd = 2, pch = 19, cex = 0.6)
  
}

############## Example with k = 4 and ngroups = 150, threshold = 20 #################

ngroups = 150
k = 4
t = 20

cluster = as.numeric(pam(diss_with_threshold[[as.character(t)]], j, diss = TRUE, variant = "faster", cluster.only = TRUE))

plot_grid(gridSizeX, gridSizeY, cluster)

repr = matrix(0, nrow = ngroups, ncol = dim(data)[2])

for(r in 1:ngroups){
  
  ind = which(cluster == r)
  
  repr[r,] = colMeans(data[ind,])
  
}

finalkm = kmeans(as.data.frame(repr), k, nstart = 100, iter.max = 100)
final_cluster = finalkm$cluster[cluster]

plot_grid(gridSizeX, gridSizeY, final_cluster, 
          title = paste("PSD[1] THRESHOLD -", k, "clusters and", ngroups, "subgroups"))

repr = matrix(0, nrow = k, ncol = dim(data)[2])

for(i in 1:k){
  ind = which(final_cluster == i)
  repr[i,] = colMeans(data[ind,])
}

matplot(t(data[sample(numNodes, size = 1000),]), type = 'l', lty = 1, col = "lightgrey", main = paste("Centroids when k =", k, "and ngroups =", ngroups), ylab = "")
matlines(t(repr), type = 'l', lty = 1, col = seq(1,k), lwd = 3)
