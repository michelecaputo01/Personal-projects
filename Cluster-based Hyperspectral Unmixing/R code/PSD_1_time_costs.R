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
library(microbenchmark)

set.seed(1234)

source("functions/direct_path.R")
source("functions/PSD_1.R")
source("functions/generate_data.R")
source("functions/PSD_1_no_tracing.R")

len = seq(5, 30, by = 5)

natt = 187

minshift = -1
maxshift = 1

sigma_f = 1    ## amplitude
scale = 1      ## bigger --> smooth

time = NULL

for(l in len){
  
  data = generate_data(l, l, natt, sigma_f, scale, minshift, maxshift, print = FALSE)$data
  
  t = microbenchmark(PSD_1(data, l, l), times = ceiling((max(len)^2)/(4*l^2)))$time
  
  time = c(time, mean(t)/1e9)
  
}

seq_of_points = seq(min(len)^2, max(len)^2, by = 5)
costant = 130
n = seq(0, max(len^2))

plot(len^2, time, xlab = "Number of points", ylab = "Time", type = 'l', col = "black",
     main = expression(bold("Time elapsed (in seconds) for PSD"[1])), lwd = 3)
points(len^2, time, pch = 19, lwd = 3)
lines(n, n/costant, lty = 2, lwd = 2, col = "forestgreen")
lines(n, (n/costant)^2, lty = 2, lwd = 2, col = "blue")
lines(n, (n/costant)^3, lty = 2, lwd = 2, col = "red")
grid()
legend("topleft", legend = c("Real complexity", "O(N)", expression(O(N^2)), expression(O(N^3))),
       col = c("black", "forestgreen", "blue", "red"), lty = c(1,2,2,2), lwd = c(3,2,2,2), 
       pch = c(19, NA, NA, NA), cex = 0.8)

######################### With tracing vs without tracing #########

time_with = NULL
time_without = NULL

for(l in len){
  
  data = generate_data(l, l, natt, sigma_f, scale, minshift, maxshift, print = FALSE)$data
  
  t = microbenchmark(PSD_1(data, l, l), times = ceiling((max(len)^2)/(4*l^2)))$time
  time_with = c(time_with, mean(t)/1e9)
  
  t = microbenchmark(PSD_1_no_tracing(data, l, l), times = ceiling((max(len)^2)/(4*l^2)))$time
  time_without = c(time_without, mean(t)/1e9)
  
}

plot(len^2, time_with, xlab = "Number of points", ylab = "Time", type = 'l', col = "black",
     main = expression(bold("Time elapsed (in seconds) for PSD"[1])), lwd = 3)
points(len^2, time_with, pch = 19, lwd = 3, col = "black")
lines(len^2, time_without, lwd = 3, col = 'blue')
points(len^2, time_without, lwd = 3, pch = 19, col = "blue")
grid()
legend("bottomright", legend = c("With tracing", "Without tracing"),
       col = c("black", "blue"), lty = 1, lwd = 3, pch = 19, cex = 0.8)