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
library(R.matlab)

set.seed(1234)

source("functions/compute_local_representative.R")
source("functions/BVKMA_cosine.R")
source("functions/boot_iter_cosine.R")
source("functions/compute_eta.R")
source("functions/compute_theta_cosine.R")
source("functions/plot_grid.R")

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
  
to_remove = c(101, 135)
first_der_data = cosine_firstdev(original_data, to_remove)

data = original_data

gridSizeX = 190
gridSizeY = 250

numNodes = gridSizeX * gridSizeY

B = 50
print = FALSE

n = c(50, 100, 150, 200, 250, 300, 400, 500, 600)
k = seq(2,8)

eta_matrix = matrix(0, nrow = length(n), ncol = length(k))
rownames(eta_matrix) = as.character(n)
colnames(eta_matrix) = as.character(k)
theta_matrix = eta_matrix

for(i in 1:length(n)){
  for(j in 1:length(k)){
    a = BVKMA_cosine(data, gridSizeX, gridSizeY, n[i], k[j], B, print = FALSE)
    eta_matrix[i,j] = a$eta
    theta_matrix[i,j] = a$theta
  }
}

colors = colorRampPalette(c("lightblue", "darkblue"))
colors = colors(length(k))

matplot(n, eta_matrix, pch = 19, type = 'p', xaxt = 'n', xlab = "n", ylab = "Normalized Entropy",
        col = colors, cex = 1.2, main = "BVKMA (Cosine) - Normalized Entropy", ylim = c(0,1))
matlines(n, eta_matrix, pch = 19, lty = 1, col = colors, lwd = 2)
axis(side = 1, at = n, labels = n)
grid()
legend("bottomright", legend = c("k = 2", "k = 3", "k = 4", "k = 5", "k = 6", "k = 7", "k = 8"), 
       col = colors, lty = 1, lwd = 2, pch = 19, cex = 0.75)


colors = colorRampPalette(c("lightblue", "darkblue"))
colors = colors(4)

matplot(k, t(theta_matrix[c("150", "200", "250", "300"),]), type = 'p', pch = 19, xlab = "Number of clusters", ylab = expression(theta * "-cosine"), cex = 1.2, 
        col = colors, main = "BVKMA (Cosine) - Theta values", ylim = c(0,1))
matlines(k, t(theta_matrix[c("150", "200", "250", "300"),]), lty = 1, lwd = 2, col = colors)
grid()
legend("topleft", legend = c("n = 150", "n = 200", "n = 250", "n = 300"), 
       col = colors, lty = 1, lwd = 2, pch = 19, cex = 0.8)

############################## Example with n = 200 and k = 3 #######################

n = 200
k = 4

cluster = BVKMA_cosine(data, gridSizeX, gridSizeY, n, k, B, print = TRUE)$cluster

plot_grid(gridSizeX, gridSizeY, cluster, title = paste("BVKMA (cosine) -", k, "clusters and n =", n),
          scaled = 0, centroids = TRUE, data = data, 
          title_centr = paste("Centroids when k =", k, "and n =", n), 
          original = TRUE, original_data = original_data)
