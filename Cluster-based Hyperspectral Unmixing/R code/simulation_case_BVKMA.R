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
library(gtools)
library(psych)
library(MASS)
library(Polychrome)

set.seed(1234)

source("functions/compute_local_representative.R")
source("functions/BVKMA.R")
source("functions/boot_iter.R")
source("functions/compute_eta.R")
source("functions/compute_theta.R")
source("functions/generate_data.R")
source("functions/kpi.R")
source("functions/plot_grid.R")
source("functions/clr.R")
source("functions/firstder.R")

gridSizeX = 50
gridSizeY = 50

numNodes = gridSizeX * gridSizeY

natt = 187

B = 100

grid = expand.grid(x = 1:gridSizeX, y = 1:gridSizeY)

# easy : -2, 2, 1, 0.75
# hard : -10, 10, 3, 0.3

minshift = -10
maxshift = 10

sigma_f = 3  ## amplitude
scale = 0.3   ## bigger --> smooth

d = generate_data(gridSizeX, gridSizeY, natt, sigma_f, scale, minshift, maxshift)
original_data = d$data
real_labels  = d$real_labels

clr_data = clr(original_data)

data = original_data

#################### BVKMA #################

n = c(10, 25, 50, 75, 100, 150, 200)
k = seq(2,8)

eta_matrix = matrix(0, nrow = length(n), ncol = length(k))
rownames(eta_matrix) = as.character(n)
colnames(eta_matrix) = as.character(k)
theta_matrix = eta_matrix

for(i in 1:length(n)){
  for(j in 1:length(k)){
    a = BVKMA(data, gridSizeX, gridSizeY, n[i], k[j], B, print = FALSE)
    eta_matrix[i,j] = a$eta
    theta_matrix[i,j] = a$theta
  }
}

colors = colorRampPalette(c("lightblue", "darkblue"))
colors = colors(length(k))

matplot(n, eta_matrix, pch = 19, type = 'p', xaxt = 'n', xlab = "n ", ylab = "Normalized Entropy",
        col = colors, cex = 1.2, main = "BVKMA (sim) - Normalized Entropy", ylim = c(0,1))
matlines(n, eta_matrix, pch = 19, lty = 1, col = colors, lwd = 2)
axis(side = 1, at = n, labels = n)
grid()
legend("bottomleft", legend = c("k = 2", "k = 3", "k = 4", "k = 5", "k = 6", "k = 7", "k = 8"), 
       col = colors, lty = 1, lwd = 2, pch = 19, cex = 0.75)

colors = colorRampPalette(c("lightblue", "darkblue"))
colors = colors(3)

matplot(k, t(theta_matrix[c("50", "75", "100"),]), type = 'p', pch = 19, xlab = "Number of clusters", ylab = expression(theta), cex = 1.2, 
        col = colors, main = "BVKMA (sim) - Theta values", ylim = c(0,1))
grid()
matlines(k, t(theta_matrix[c("50", "75", "100"),]), lty = 1, lwd = 2, col = colors)
legend("topleft", legend = c("n = 50", "n = 75", "n = 100"), 
       col = colors, lty = 1, lwd = 2, pch = 19, cex = 0.6)

#################### Example with n = 100 and k = 4 #############

n = 100
k = 4

cluster = BVKMA(data, gridSizeX, gridSizeY, n, k, B, print = TRUE)$cluster

kpi_value = kpi(real_labels, cluster)

plot_grid(gridSizeX, gridSizeY, cluster, title = paste("CLR + BVKMA, hard - KPI =", kpi_value),
          scaled = 0, centroids = TRUE, data = data, 
          title_centr = paste("Centroids when k =", k, "and n =", n), 
          original = TRUE, original_data = original_data)
