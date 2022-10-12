# R Clustering Algorithms
Clustering is an unsupervised machine learning technique that groups data points into clusters by similar features. This repository contains a <a href="https://github.com/TomMakesThings/R-Clustering-Algorithms/blob/main/ClusteringAlgorithms.R">script</a> in which k-means, DBSCAN and agglomerative hierarchical clustering algorithms are implemented from scratch in R. These are tested using the well known Iris dataset and scRNA-seq data by <a href="https://github.com/LuyiTian/sc_mixology">Tian et al</a>.

## Clustering Algorithms
#### K-Means
K-means is a non-deterministic, centroid based clustering method that attempts to group similar data points into a set number of clusters k. The algorithm begins by setting k randomly placed centroids, in this case through randomly selecting k points, then iteratively optimising them until the positions converge or a max number of iterations is reached. During optimisation, the Euclidean distance between each point and each centroid is calculated to find its closest cluster. At the end of each iteration, the centroid positions are updated as the mean of their assigned points. Once convergence is reached, the cluster assignment for each point is returned.

#### Agglomerative Hierarchical Clustering
Hierarchical clustering is a type of deterministic clustering algorithm that repeatedly merges or splits of clusters to form a dendrogram. The two approaches referred to as top-down (divisive) and bottom-up (agglomerative). Here, agglomerative was selected as it is more frequently used and is more computationally efficient. To implement the algorithm, each point is initially set as its own cluster. The Euclidean distance between all clusters is calculated using a distance matrix and the closest
two clusters are combined. This process is average group linkage, a popular distance metric that avoids creating either large or compact clusters. The linkage process repeats until the desired number of clusters is obtained and the cluster assignments returned.

#### DBSCAN
Density-Based Spatial Clustering of Applications with Noise is an algorithm better suited to clustering arbitrary shapes with varying densities than k-means or hierarchical clustering. Unlike hierarchical clustering, it is not completely deterministic as border points can be clustered differently depending on the order data is input. To begin, a point is selected at random and the distance between it and all other points calculated. Like with the other two algorithms, the distance metric is set as Euclidean as this is suitable for relatively low dimensionality, though in theory other measures such as Manhattan could be used. Points within a distance ε of the selected point are considered to be neighbours and are added to form a candidate cluster. The neighbourhood for each newly added point is located and added to the candidate cluster, then the process repeated until no new neighbouring points are found. If a minimum number of points were located, the points are considered a cluster, otherwise they are marked as noise. While there remains points that have not been clustered or set as noise, a new point is again picked at random and the process of finding neighbouring points to form a new candidate cluster repeated. Finally once all points have been marked, the labels are returned.

## Data

#### Iris
The Iris dataset contains three species of flower with 50 instances of each. Samples have four features representing their sepal and petal length and width. One species is linearly separable from the others, but the other two are not linearly separable from one another.

<p><img width = 600 src="https://github.com/TomMakesThings/R-Clustering-Algorithms/blob/assets/Images/Iris-UMAP-Results.png"></p>

#### Cell Line
The scRNA-seq data consists of counts for 11,786 genes for 3,918 cell across five cancer cell lines (H838, H2228, HCC827, A549, H1975).

<p><img width = 600 src="https://github.com/TomMakesThings/R-Clustering-Algorithms/blob/assets/Images/Cell-tSNE-Results.png"></p>
