# Cluster-based Hyperspectral Unmixing
In this folder the code developed as part of my university thesis is reported. In particular, the folder is divided into two subfolders that contain the MATLAB code, used for hyperspectral unmixing, and the R code, used to explore the possibilities related to clustering.

## Matlab code
The aim of this work is to improve the performance of an existing hyperspectral unmixing algorithm by using clustering to better separate the data points. After the cluster unmixing, which is carried out independently, a reunification step is performed. 
The folder contains the existing hyperspectral unmixing functions along with some auxiliary scripts, as well as comparisons between the two approaches (with and without clustering). The reported clusterization is performed using the functions in the "R code" folder.

## R code
In this folder the problem of "spatial" clustering is explored. The distinguishing feature of these algorithms is that they perform clustering of spatial functional data taking into account both the non-spatial variables (in this case, a function of the reflectance) and the coordinates of the data points. The algorithms explored are: 
- Bagging-Voronoi distance (P. Secchi, S. Vantini, and V. Vitelli. Bagging voronoi classifiers for clustering spatial
functional data. International Journal of Applied Earth Observation and Geoinformation, 22:53–64, 2013);
- Penalized Spatial Distance (B. Zhang, W. J. Yin, M. Xie, and J. Dong. Geo-spatial clustering with non-spatial
attributes and geographic non-overlapping constraint: A penalized spatial distance
measure. In Advances in Knowledge Discovery and Data Mining, pages 1072–1079,
Berlin, Heidelberg, 2007. Springer Berlin Heidelberg);
- Classic K-means clustering, where the distance is computed as combination of functional and spatial attributes).
