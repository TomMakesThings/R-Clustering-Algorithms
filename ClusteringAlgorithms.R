# Euclidean distance: https://www.analyticsvidhya.com/blog/2020/02/4-types-of-distance-metrics-in-machine-learning/
# https://towardsdatascience.com/the-5-clustering-algorithms-data-scientists-need-to-know-a36d136ef68

library(dplyr)
library(clue)
library(fossil)
library(ggplot2)
library(clValid)

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

# Bottom-up hierarchical clustering
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
  
  # Create list of cluster assignment per points
  point_clusters <- paste("cluster", 1:length(clusters), sep = "_")
  names(point_clusters) <- rownames(points)
  
  distances <- c()
  
  # Create dataframe of average cluster positions
  cluster_averages <- do.call(rbind.data.frame, clusters)
  cluster_averages$id <- NULL
  
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
    
    distances <- c(distances, min(average_cluster_dist))

    # Merge both into the first cluster
    clusters[[cluster_1]] <- rbind(clusters[[cluster_1]], clusters[[cluster_2]])
    
    # Update cluster assignment of each point
    for (p in 1:nrow(clusters[[cluster_2]])) {
      point_clusters[clusters[[cluster_2]][p,]$id] <- cluster_1
    }
    
    # Remove the second cluster
    clusters[[cluster_2]] <- NULL
    cluster_averages <- cluster_averages[-c(min_idx[2]),]
    
    # Update the center point of the combined cluster
    cluster_averages[cluster_1,] <- colMeans(clusters[[cluster_1]][1:(ncol(clusters[[cluster_1]])-1)])
  }
  
  results <- list(clusters = factor(as.numeric(factor(unname(point_clusters)))),
                  distances = distances)
  
  return(results)
}

# Density-Based spatial clustering of applications with noise (DBSCAN)
DBSCAN <- function(points, epsilon, min_points) {
  # Initialise list to store all clusters and points assigned as noise
  clusters <- list()
  noise <- NULL
  
  # Randomly select a point to form the first cluster
  shuffled_points <- sample_n(data.frame(points), nrow(points))
  current_cluster <- NULL
  new_points <- NULL
  previous_points <- shuffled_points[1,]
  remaining_points <- shuffled_points[2:nrow(shuffled_points),]

  # Iterate until all points are classified into clusters or noise
  while (nrow(remaining_points) > 0) {
    
    # Compare newly clustered points to those unclustered
    for (x in 1:nrow(previous_points)) {
      point_1 <- previous_points[x,]
      for (y in 1:nrow(remaining_points)) {
        point_2 <- remaining_points[y,]
        
        # Check the point has not yet been added
        if (!(rownames(point_2) %in% rownames(new_points))) {
          # Calculate Euclidean distance between the two points
          euclidean <- calculateEuclidean(point_1, point_2)[,]
          
          # Check if the distance is less than the threshold epsilon
          if (euclidean < epsilon) {
            # If so, record the point to be added to the cluster
            if (is.null(new_points)) {
              new_points <- point_2
            } else {
              new_points <- rbind(new_points, point_2)
            }
          }
        }
      }
    }
    
    # Add previously discovered points to the current cluster
    if (is.null(current_cluster)) {
      current_cluster <- previous_points
    } else {
      current_cluster <- rbind(current_cluster, previous_points)
    }

    # Check if any new points were clustered
    if (!is.null(new_points)) {

      # Update remaining unclassified points by removing newly clustered points
      remaining_rows <- rownames(remaining_points)[!rownames(remaining_points) %in% rownames(new_points)]
      remaining_points <- remaining_points[remaining_rows,]
      
      # Reset point placeholders
      previous_points <- new_points
      new_points <- NULL
      
    } else {
      
      # If no new points clustered, check if current cluster has n points or more
      if (nrow(current_cluster) >= min_points) {
        # If so, record the cluster
        clusters <- append(clusters, list(current_cluster))
      } else {
        # Otherwise, label the points as noise
        if (is.null(noise)) {
          noise <- current_cluster
        } else {
          noise <- rbind(noise, current_cluster)
        }
      }
      
      # Reset variables
      current_cluster <- NULL
      new_points <- NULL
      previous_points <- remaining_points[1,]
      remaining_points <- remaining_points[!(row.names(remaining_points) %in% rownames(previous_points)),]
    }
  }
  
  return(list(clusters, noise))
}

# Extract features and metadata
iris_features <- iris[, !names(iris) %in% c("Species")]
iris_metadata <- iris[, "Species"]

# Run principal component analysis with normalization
iris_pca <- prcomp(iris_features, center = TRUE, scale. = TRUE)$x
row.names(iris_pca) <- row.names(iris)

# Run k-means clustering
kmean_clusters <- k_means(iris_pca, max_iterations = 50)

# Run agglomerative hierarchical clustering
hierarchical_results <- agglomerativeClustering(iris_pca)
hierarchical_clusters <- hierarchical_results$clusters

# Run DBSCAN clustering
dbscan <- DBSCAN(iris_pca, 0.6, 3)

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
hierarchical_map <- clusterLabelMatch(iris_metadata, hierarchical_clusters)

convertTruth <- function(labels, label_to_cluster_map) {
  unlist(lapply(labels, function(x) label_to_cluster_map[[x]]))
}

kmeans_truth <- convertTruth(iris_metadata, kmean_map)
hierarchical_map_truth <- convertTruth(iris_metadata, hierarchical_map)

# Calculate adjusted Rand index
kmeans_adj_rand <- adj.rand.index(kmeans_truth, kmean_clusters)
hierarchical_adj_rand <- adj.rand.index(hierarchical_map_truth, hierarchical_clusters)

# Create graphs
iris_pca_clusters <- data.frame(iris_pca)
iris_pca_clusters$label <- iris_metadata
iris_pca_clusters$kmeans <- kmean_clusters
iris_pca_clusters$hierarchical <- hierarchical_clusters

ggplot(iris_pca_clusters, aes(x = PC1, y = PC2, color = kmeans, shape = label)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#fc8021", "#462cc7", "#3ab03a")) +
  labs(title = paste("K-Means Clustering - ARI", signif(kmeans_adj_rand, 3)))

ggplot(iris_pca_clusters, aes(x = PC1, y = PC2, color = hierarchical, shape = label)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#fc8021", "#462cc7", "#3ab03a")) +
  labs(title = paste("Agglomerative Hierarchical Clustering - ARI",
                     signif(hierarchical_adj_rand, 3)))