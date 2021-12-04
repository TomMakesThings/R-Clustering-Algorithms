# Euclidean distance: https://www.analyticsvidhya.com/blog/2020/02/4-types-of-distance-metrics-in-machine-learning/

library(clue)
library(fossil)
library(ggplot2)

# Calculate Euclidean distance for n-dimensions
calculateEuclidean <- function(p1, p2) {
  # Add the squared difference of each dimension in both points
  total <- 0
  
  for (i in 1:length(names(p1))) {
    total <- total + ((p1[i] - p2[i]) ** 2)
  }
  
  # Square root the summed difference
  euclidean_dist <- unname(sqrt(total))
  
  return(euclidean_dist)
}

formatClusterResults <- function(centroid_assignments, points) {
  # Add a columns to the results of cluster assignment
  for (i in 1:length(centroid_assignments)) {
    cluster_points_df <- data.frame(centroid_assignments[[i]])
    cluster_points_df$cluster <- i
    centroid_assignments[[i]] <- cluster_points_df
  }
  
  results <- do.call("rbind", centroid_assignments)
  
  if (is.null(row.names(points))) {
    row.names(results) <- 1:nrow(results)
  } else {
    row.names(results) <- row.names(points)
  }
  
  return(results)
}

# K-means clustering algorithm
k_means <- function(points, k = 3, max_iterations = 100) {
  # Randomly initialize k centroids
  centroids <- points[sample(1:nrow(points), k), ]
  
  # Optimize centers up to n times
  for (i in 1:max_iterations) {
    # Initialize list to record points assigned to each centroid
    centroid_assignments <- vector("list", length = k)
    point_clusters <- c()
    
    # For each point, find its closest centroid
    for (p in 1:nrow(points)) {
      point <- points[p,]
      
      min_dist <- Inf
      best_centroid <- NA
      
      # Test each centroid
      for (c in 1:k) {
        centroid <- centroids[c,]
        # Calculate Euclidean distance between point and centroid
        euclidean <- calculateEuclidean(point, centroid)
        
        # If point is closest to this centroid than previous attempts, record the centroid
        if (euclidean < min_dist) {
          min_dist <- euclidean
          best_centroid <- c
        }
      }
      
      # Add the point's coordinates to its best centroid
      if (is.null(centroid_assignments[best_centroid][[1]])) {
        centroid_assignments[[best_centroid]] <- point
      } else {
        centroid_assignments[[best_centroid]] <- rbind(centroid_assignments[[best_centroid]], point)
      }
      
      point_clusters <- c(point_clusters, best_centroid)
    }
    
    # Record current centroid assignment
    centroids_old <- centroids
    
    # Calculate new centroids
    for (c in 1:k) {
      centroids[c,] <- colMeans(centroid_assignments[[c]])
    }
    
    # If centroids are the same, algorithm has converged so stop
    if (identical(centroids, centroids_old)) {
      return(as.factor(point_clusters))
    }
  }
  
  return(as.factor(point_clusters))
}

# Bottom up hierarchical clustering
agglomerativeClustering <- function(points, n_clusters = 3) {
  # Initialise list to store clusters
  clusters <- list()
  
  # Set each point as a new cluster
  for (i in 1:nrow(points)) {
    clusters[[i]] <- data.frame(t(points[i,]))
    clusters[[i]]$id <- i
  }
  
  # Set the cluster index names
  names(clusters) <- paste("cluster", 1:length(clusters), sep = "_")
  
  # Create dataframe of average cluster positions
  cluster_averages <- do.call(rbind.data.frame, clusters)
  cluster_averages$id <- NULL
  
  #return(cluster_averages)
  
  # Iteratively use average group linkage to merge clusters
  for (j in 1:(nrow(points) - n_clusters)) {
    # Calculate Euclidean distance between every cluster
    average_cluster_dist <- as.matrix(dist(cluster_averages))
    # Set distance matrix diagonal to infinity to stop any cluster being compare to itself
    diag(average_cluster_dist) <- Inf
    # Find the two clusters with the smallest distance
    min_idx <- arrayInd(which.min(average_cluster_dist), dim(average_cluster_dist))
    cluster_1 <- rownames(average_cluster_dist)[min_idx[1]]
    cluster_2 <- rownames(average_cluster_dist)[min_idx[2]]

    # Merge both into the first cluster
    clusters[[cluster_1]] <- rbind(clusters[[cluster_1]], clusters[[cluster_2]])
    # Remove the second cluster
    clusters[[cluster_2]] <- NULL
    cluster_averages <- cluster_averages[-c(min_idx[2]),]

    print("now")
    print(clusters[cluster_1])
    print(clusters[cluster_2])
    
    # Update the center point of the combined cluster
    cluster_averages[cluster_1,] <- colMeans(clusters[[cluster_1]][1:(ncol(clusters[[cluster_1]])-1)])
  }
  
  return(clusters)
}


# Extract features and metadata
iris_features <- iris[, !names(iris) %in% c("Species")]
iris_metadata <- iris[, "Species"]

# Run principal component analysis with normalization
iris_pca <- prcomp(iris_features, center = TRUE, scale. = TRUE)$x
row.names(iris_pca) <- row.names(iris)

iris_pca <- iris_pca[,c(1,2)]

# Run k-means clustering
kmean_clusters <- k_means(iris_pca, max_iterations = 50)

# Run agglomerative hierarchical clustering
hierarchical_clusters <- agglomerativeClustering(iris_pca)


clusters <- lapply(clusters, function(x) split(x, 2, 1:length(clusters)))
clusters <- lapply(clusters, function(x) setNames(x, "points"))

for (c in 1:length(hierarchical_clusters)) {
  # Calculate average inter-cluster distance
  hierarchical_clusters[[c]]$average_dist <- calculateEuclidean(hierarchical_clusters[[c]]$point,
                                                                hierarchical_clusters[[c]]$point)
}

iris_pca_clusters <- data.frame(iris_pca)
iris_pca_clusters$label <- iris_metadata
iris_pca_clusters$kmeans <- kmean_clusters

ggplot(iris_pca_clusters, aes(x = PC1, y = PC2, color = kmeans, shape = label)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#fc8021", "#462cc7", "#3ab03a")) +
  labs(title = "K-Means Clustering")

# Find the best mapping between labels and predicted clusters using the Hungarian matching algorithm
clusterLabelMatch <- function(labels, predictions) {
  n_clusters <- length(levels(predictions))
  n_samples <- length(predictions)
  
  # Construct a bipartite graph via an adjacency matrix
  bipartite_graph <- matrix(0, n_clusters, n_clusters)
  
  for (i in 1:n_samples) {
    # Add 1 for every intersection between rows and columns
    bipartite_graph[predictions[i], labels[i]] <- bipartite_graph[predictions[i], labels[i]] + 1
  }
  # Solve the linear sum assignment problem through the Hungarian method
  best_assignment <- solve_LSAP(max(bipartite_graph) - bipartite_graph, maximum = FALSE)
  
  # Create a list mapping labels to cluster predictions
  assignment_map <- list() 
  
  for (j in 1:length(best_assignment)) {
    assignment_map[levels(labels)[best_assignment[j]]] <- j
  }
  
  return(assignment_map)
}

kmean_map <- clusterLabelMatch(iris_metadata, kmean_clusters)

convertTruth <- function(labels, label_to_cluster_map) {
  unlist(lapply(labels, function(x) label_to_cluster_map[[x]]))
}

kmeans_truth <- convertTruth(iris_metadata, kmean_map)

kmeans_adj_rand <- adj.rand.index(kmeans_truth, kmean_clusters)
kmeans_adj_rand