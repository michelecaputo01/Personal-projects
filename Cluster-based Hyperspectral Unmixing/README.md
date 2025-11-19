# Cluster-based Hyperspectral Unmixing
In this folder the code developed as part of my university thesis is reported. The aim of this work is to accurately estimate the physical composition of a given terrain—typically the soil of celestial bodies—and the percentage presence of each material contained within it. 
The folder is divided into two subfolders that contain the MATLAB code, used for hyperspectral unmixing, and the R code, used to explore the possibilities related to clustering.

## Matlab code
The aim of this work is to improve the performance of an existing hyperspectral unmixing algorithms (\textit(J. M. Bioucas-Dias and M. A. T. Figueiredo. Alternating direction algorithms for constrained sparse regression: Application to hyperspectral unmixing. In 2010 2nd Workshop on Hyperspectral Image and Signal Processing: Evolution in Remote Sensing, pages 1–4, 2010)) by using clustering to better separate the data points. After the cluster unmixing, which is carried out independently, a reunification step is performed. 
The folder contains the existing hyperspectral unmixing functions along with some auxiliary scripts, as well as comparisons between the two approaches (with and without clustering). The reported clusterization is performed using the functions in the "R code" folder.

## R code
In this folder the problem of "spatial" clustering is explored. The distinguishing feature of these algorithms is that they perform clustering of spatial functional data taking into account both the non-spatial variables (in this case, a function of the reflectance) and the coordinates of the data points. The algorithms explored are: 
- Bagging-Voronoi classifier (P. Secchi, S. Vantini, and V. Vitelli. Bagging voronoi classifiers for clustering spatial
functional data. International Journal of Applied Earth Observation and Geoinformation, 22:53–64, 2013);
- Classic K-means clustering using Penalized Spatial Distance (B. Zhang, W. J. Yin, M. Xie, and J. Dong. Geo-spatial clustering with non-spatial
attributes and geographic non-overlapping constraint: A penalized spatial distance
measure. In Advances in Knowledge Discovery and Data Mining, pages 1072–1079,
Berlin, Heidelberg, 2007. Springer Berlin Heidelberg);
- Classic K-means clustering, where the distance is computed as combination of functional and spatial attributes).
