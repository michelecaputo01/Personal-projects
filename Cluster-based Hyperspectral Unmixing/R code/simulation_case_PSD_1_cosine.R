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
library(igraph)
library(gtools)
library(psych)

set.seed(1234)

source("functions/direct_path.R")
source("functions/compute_theta_cosine.R")
source("functions/PSD_1.R")
source("functions/generate_data.R")
source("functions/kpi.R")
source("functions/plot_grid.R")
source("functions/firstder.R")

gridSizeX = 50
gridSizeY = 50

numNodes = gridSizeX * gridSizeY

natt = 187

grid = expand.grid(x = 1:gridSizeX, y = 1:gridSizeY)

minshift = -2
maxshift = 2

sigma_f = 1    ## amplitude
scale = 0.75      ## bigger --> smooth

d = generate_data(gridSizeX, gridSizeY, natt, sigma_f, scale, minshift, maxshift)
original_data = d$data
real_labels = d$real_labels

first_der_data = firstder(original_data)

data = first_der_data

###################################################################

dissimilarity = PSD_1(data, gridSizeX, gridSizeY, method = "cosine")

ngroups_first = c(10, 25, 50, 75, 100, 150, 200, 250, 300)

theta = NULL

cluster_matrix = matrix(0, nrow = length(ngroups_first), ncol = numNodes)
rownames(cluster_matrix) = ngroups_first

for(n in ngroups_first){
  
  cat(sprintf("\rProgress:  ngroups = %d", n))
  flush.console()
  
  cluster = as.numeric(pam(dissimilarity, n, diss = TRUE, variant = "faster", cluster.only = TRUE))
  
  theta = c(theta, compute_theta_cosine(data, cluster))
  
  cluster_matrix[as.character(n),] = cluster
  
}

plot(ngroups_first, theta, type = 'l', col = "black", xaxt = 'n', ylab = expression(theta), 
     xlab = "Subgroups", main = "Theta values after single PAM", ylim = c(0,1))
points(ngroups_first, theta, pch = 19)
axis(side = 1, at = ngroups_first, labels = ngroups_first)
grid()

k_values = 2:8
ngroups = c(150, 200, 250)
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
    
    cl = as.numeric(pam(repr_dist, i, diss = TRUE, nstart = 30, cluster.only = TRUE))
    final_cluster = cl[cluster]
    
    theta_matrix[as.character(i), as.character(j)] = compute_theta_cosine(data, final_cluster)
    
  }
  
}

colors = colorRampPalette(c("lightblue", "darkblue"))
colors = colors(length(ngroups))

matplot(k_values, theta_matrix, type = 'p', xaxt = 'n', pch = 19, cex = 1.2, ylab = expression(theta), 
        col = colors, main = "PSD (cosine) (sim) - Theta values", ylim = c(0,1), xlab = "Number of clusters")
matlines(k_values, theta_matrix, lty = 1, lwd = 2, col = colors)
axis(side = 1, labels = k_values, at = k_values)
grid()
legend("bottomright", legend = c("n = 150", "n = 200", "n = 250"), 
       col = colors, lty = 1, lwd = 2, pch = 19, cex = 0.6)

############## Example with k = 4 and ngroups = 200 #################

k = 4
ngroups = 200

cluster = cluster_matrix[as.character(ngroups),]

plot_grid(gridSizeX, gridSizeY, cluster, title = paste0("After single PAM (", ngroups, " subgroups)"),
          scaled = 0, centroids = FALSE, original = FALSE)

repr = matrix(0, nrow = ngroups, ncol = dim(data)[2])

for(t in 1:ngroups){
  
  ind = which(cluster == t)
  
  repr[t,] = colMeans(data[ind,])
  
}

repr_dist = proxy::dist(repr, method = "cosine")

cl = as.numeric(pam(repr_dist, k, diss = TRUE, nstart = 30, cluster.only = TRUE))
final_cluster = cl[cluster]

kpi_value = kpi(real_labels, final_cluster)

plot_grid(gridSizeX, gridSizeY, final_cluster, title = paste("First der. + PSD cosine, hard - KPI =", kpi_value),
          scaled = 0, centroids = TRUE, data = data, 
          title_centr = paste("Centroids when k =", k, "and ngroups =", ngroups), 
          original = TRUE, original_data = original_data)
