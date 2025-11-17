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

source("functions/direct_path.R")
source("functions/direct_path1.R")
source("functions/direct_path2.R")

gridSizeX = gridSizeY = 30
numNodes = gridSizeX * gridSizeY
grid = expand.grid(x = 1:gridSizeX, y = 1:gridSizeY)

################### Direct_path ##############

microbenchmark({
  for(i in 1:(numNodes-1)){
    cat(sprintf("\rProgress: %d/%d               ", i+1, numNodes))
    flush.console()
    for(j in (i+1):numNodes){
      path = direct_path(i, j, grid, gridSizeX)
    }
    
  }
}, times = 1)

plot(grid, asp = 1, main = "From 30 to 700, Angular Coefficient Method", col = "lightgrey")
points(grid[direct_path(30,700, grid, gridSizeX),], col = "red", pch = 19)
lines(grid[c(30,700),], lwd = 2, lty = 2)

################### Direct_path1 ##############

microbenchmark({
  for(i in 1:(numNodes-1)){
    cat(sprintf("\rProgress: %d/%d               ", i+1, numNodes))
    flush.console()
    for(j in (i+1):numNodes){
      path = direct_path1(i, j, grid, gridSizeX)
    }
    
  }
}, times = 1)

plot(grid, asp = 1, main = "From 30 to 700, Bresenham Algorithm", col = "lightgrey")
points(grid[direct_path1(30, 700, grid, gridSizeX),], col = "red", pch = 19)
lines(grid[c(30,700),], lwd = 2, lty = 2)

################### Direct_path2 ##########

microbenchmark({
  for(i in 1:(numNodes-1)){
    cat(sprintf("\rProgress: %d/%d               ", i+1, numNodes))
    flush.console()
    for(j in (i+1):numNodes){
      path = direct_path2(i, j, grid, gridSizeX)
    }
    
  }
}, times = 1)

plot(grid, asp = 1, main = "From 30 to 700, points on the Segment", col = "lightgrey")
points(grid[direct_path2(30, 700, grid, gridSizeX),], col = "red", pch = 19)
lines(grid[c(30,700),], lwd = 2, lty = 2)
