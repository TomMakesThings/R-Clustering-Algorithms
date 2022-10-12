library(dplyr)
library(clue)
library(fossil)
library(ggplot2)
library(clValid)
library(umap)
library(Rtsne)
library(gridExtra)
library(cluster)

# Set directory to be repository containing code
setwd("C:/Users/redds/Documents/GitHub/R-Clustering-Algorithms")

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

### Clustering algorithms ###

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
        
        # If point is closest to this centroid than previous attempts,
        # record the centroid
        if (euclidean < min_dist) {
          min_dist <- euclidean
          best_centroid <- c
        }
      }
      
      # Add the point's coordinates to its best centroid
      if (is.null(centroid_assignments[best_centroid][[1]])) {
        centroid_assignments[[best_centroid]] <- point
      } else {
        centroid_assignments[[best_centroid]] <- rbind(centroid_assignments[[best_centroid]],
                                                       point)
      }
      
      point_clusters <- c(point_clusters, best_centroid)
    }
    
    # Record current centroid assignment
    centroids_old <- centroids
    
    # Calculate new centroids
    for (c in 1:k) {
      if (!is.null(centroid_assignments[[c]])) {
        centroids[c,] <- colMeans(centroid_assignments[[c]])
      }
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
    # Set distance matrix diagonal to infinity to stop any cluster
    # being compare to itself
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
  
  # Assign cluster numbers to the input points
  clustered_points <- data.frame(points)
  clustered_points$cluster <- NA
  
  for (i in 1:length(clusters)) {
    for (j in 1:nrow(clusters[[i]]))
      clustered_points[rownames(clusters[[i]][j,]),]$cluster <- i
  }
  
  return(factor(clustered_points$cluster))
}

### Data pre-processing ###

# Apply dimensionality reduction algorithm
reduceDimensions <- function(data, reduction = "pca", n_pcs = 0, perplexity = 30) {
  if (reduction == "umap") {
    # UMAP dimensionality reduction
    features = umap(data)$layout
    colnames(features) <- c("UMAP1", "UMAP2")
    
  } else if (reduction == "tsne") {
    # t-SNE dimensionality reduction
    features <- Rtsne(distinct(data), dims = 2, perplexity = perplexity,
                      max_iter = 500)$Y
    colnames(features) <- c("tSNE1", "tSNE2")
    
  } else {
    # PCA with normalization
    features <- prcomp(data, center = TRUE, scale. = TRUE)$x
    if (n_pcs > 0) {
      # Extract n principal components
      features <- features[,1:n_pcs]
    }
  }
  
  row.names(features) <- row.names(data)
  
  return(features)
}

### Running clustering ###

# Find the best mapping between labels and predicted clusters using the
# Hungarian matching algorithm
clusterLabelMatch <- function(labels, predictions) {
  n_clusters <- length(levels(predictions))
  n_samples <- length(predictions)
  
  # Construct a bipartite graph via an adjacency matrix
  bipartite_graph <- matrix(0, n_clusters, n_clusters)
  
  for (i in 1:n_samples) {
    # Add 1 for every intersection between rows and columns
    bipartite_graph[predictions[i], labels[i]] <- bipartite_graph[predictions[i],
                                                                  labels[i]] + 1
  }
  # Solve the linear sum assignment problem through the Hungarian method
  best_assignment <- solve_LSAP(max(bipartite_graph) - bipartite_graph,
                                maximum = FALSE)
  
  # Create a list mapping labels to cluster predictions
  assignment_map <- list() 
  
  for (j in 1:length(best_assignment)) {
    assignment_map[levels(labels)[best_assignment[j]]] <- j
  }
  
  # Use the assignment map to convert labels into cluster numbers
  truth_clusters <- unlist(lapply(labels, function(x) assignment_map[[x]]))
  
  return(list(assignment_map = assignment_map, truth_clusters = truth_clusters))
}

runClustering <- function(data, metadata, algorithm, epsilon = 0.8, 
                          min_points = 3) {
  # Record time
  start_time <- proc.time()
  
  if (algorithm == "dbscan") {
    # Run density based clustering
    clusters <- DBSCAN(data, epsilon, min_points)
  } else if (algorithm == "hierarchical") {
    # Run hierarchical clustering
    hierarchical_results <- agglomerativeClustering(data,
                                                    n_clusters = length(levels(metadata)))
    clusters <- hierarchical_results$clusters
  } else {
    # Run k-means clustering
    clusters <- k_means(data, max_iterations = 50,
                        k = length(levels(metadata)))
  }
  
  # Calculate run time
  run_time <- proc.time() - start_time
  
  # # Use the ground truth to find the best matching labels to the clusters
  # expected_clusters <- clusterLabelMatch(metadata, clusters)$truth_clusters
  
  # Calculate ARI between labels and predictions
  adjusted_rand <- adj.rand.index(as.numeric(metadata), clusters)
  # Calculate average silhouette coefficient
  if (algorithm == "dbscan") {
    # Remove points classed as noise
    silhouette_scores <- silhouette(as.numeric(clusters[which(!is.na(clusters))]),
                                    dist(data[-which(is.na(as.numeric(clusters))),]))
  } else {
    silhouette_scores <- silhouette(as.numeric(clusters), dist(data))
  }
  
  suppressWarnings(if (!is.na(silhouette_scores)) {
    average_silhouette = mean(silhouette_scores[, 3])
  } else {
    average_silhouette = NA
  })
  
  return(list(clusters = clusters, run_time = run_time,
              ari = adjusted_rand, silhouette = average_silhouette))
}

### Plotting ###

plotClusters <- function(data, algorithm, colours, title = NA,
                         silhouette_score = NA, ari = NA) {
  # Set the algorithm name for the title
  if (is.na(algorithm)) {
    plot <- ggplot(data,
                   aes_string(x = names(data)[1],
                              y = names(data)[2],
                              color = "label"))
  } else {
    plot <- ggplot(data,
                   aes_string(x = names(data)[1],
                              y = names(data)[2],
                              color = algorithm,
                              shape = "label"))
    
    if (algorithm == "kmeans") {
      algorithm_title = "K-Means"
    } else if (algorithm == "hierarchical") {
      algorithm_title = "Agglomerative Hierarchical"
    } else {
      algorithm_title = "DBSCAN"
    }
  }
  
  # Create scatter graph
  plot <- plot + geom_point(size = 2) +
    scale_color_manual(values = colours)
  
  
  if (is.na(algorithm)) {
    plot <- plot + labs(title = title, color = 'Label',
                        x = names(data[1]), y = names(data[2]))
  } else {
    plot <- plot + labs(title = paste(algorithm_title, " \nSilhouette: ",
                                      signif(silhouette_score, 3),
                                      ", ARI: ", signif(ari, 3), sep = ""),
                        x = names(data[1]), y = names(data[2]),
                        color = 'Cluster', shape = 'Ground Truth')
  }
  
  return(plot)
}

plot_colours <- c("#fc8021", "#462cc7", "#3ab03a", "#ff2181", "#385df2", 
                  "#db43fa", "#5ce1ed")

### Iris clustering ###

# Number of samples use
n_iris_samples <- 150
unique_iris <- distinct(iris)

if (n_iris_samples < 150) {
  # Randomly select n samples
  unique_iris <- sample_n(unique_iris, n_iris_samples)
}

# Extract features and metadata
iris_features <- unique_iris[, !names(unique_iris) %in% c("Species")]
iris_metadata <- unique_iris[, "Species"]

plotIris <- function(colours) {
  # Plot data features
  pairs(iris_features, col = colours, lower.panel = NULL,
        cex.labels = 1, pch = 19, cex = 0.5)
  # Add species legend
  par(xpd = TRUE)
  legend(x = 0.05, y = 0.4, cex = 0.6,
         legend = as.character(levels(iris_metadata)),
         fill = unique(colours))
  par(xpd = NA)
}

# Create plot of original iris features
pdf("iris_features.pdf", width = 5, height = 4)
plotIris(colours = rev(plot_colours)[1:3][as.numeric(iris_metadata)])
dev.off()

iris_pca <- reduceDimensions(iris_features, reduction = "pca", n_pcs = 2)
iris_tsne <- reduceDimensions(iris_features, reduction = "tsne", perplexity = 30)
iris_umap <- reduceDimensions(iris_features, reduction = "umap")

# Set the dimensionality reduction to use
iris_reduction <- as.matrix(iris_features[c(1,2)])
iris_reduction <- iris_umap

iris_kmeans <- runClustering(iris_reduction, iris_metadata, algorithm = "kmeans")
iris_hierarchical <- runClustering(iris_reduction, iris_metadata, algorithm = "hierarchical")
iris_dbscan <- runClustering(iris_reduction, iris_metadata, algorithm = "dbscan",
                             epsilon = 0.8, min_points = 3)

# Create graph data
iris_plot_data <- data.frame(iris_reduction)
iris_plot_data$label <- iris_metadata
iris_plot_data$kmeans <- iris_kmeans$clusters
iris_plot_data$hierarchical <- iris_hierarchical$clusters
iris_plot_data$dbscan <- iris_dbscan$clusters

# Plot ground truth
iris_label_plot <- plotClusters(data = iris_plot_data,
                                algorithm = NA,
                                colours = rev(plot_colours)[1:length(levels(iris_plot_data$kmeans))],
                                title = "Iris Ground Truth")

# Plot clustering results
iris_kmeans_plot <- plotClusters(data = iris_plot_data,
                                 algorithm = "kmeans",
                                 colours = plot_colours[1:length(levels(iris_plot_data$kmeans))],
                                 silhouette_score = iris_kmeans$silhouette,
                                 ari = iris_kmeans$ari)

iris_hierarchical_plot <- plotClusters(data = iris_plot_data,
                                       algorithm = "hierarchical",
                                       colours = plot_colours[1:length(levels(iris_plot_data$hierarchical))],
                                       silhouette_score = iris_hierarchical$silhouette,
                                       ari = iris_hierarchical$ari)

iris_dbscan_plot <- plotClusters(data = iris_plot_data,
                                 algorithm = "dbscan",
                                 colours = plot_colours[1:length(levels(iris_plot_data$dbscan))],
                                 silhouette_score = iris_dbscan$silhouette,
                                 ari = iris_dbscan$ari)

# Add graphs as subplots
grid.arrange(iris_label_plot, iris_kmeans_plot,
             iris_hierarchical_plot, iris_dbscan_plot, nrow = 2)

# Save as PDF
pdf("iris_clusters_umap.pdf", width = 10, height = 8)
grid.arrange(iris_label_plot, iris_kmeans_plot,
             iris_hierarchical_plot, iris_dbscan_plot, nrow = 2)
dev.off()

# Plot pre-recorded algorithm run times
time_data <- data.frame(n_samples = rep(c(25, 50, 75, 100, 125, 150), 3),
                        time = c(0, 0.01, 0.01, 0.02, 0.02, 0.03,
                                 0.01, 0.02, 0.04, 0.05, 0.07, 0.08,
                                 0.18, 0.64, 1.34, 2.42, 3.57, 5.21),
                        algorithm = c(rep(c("k-means"), 6),
                                      rep(c("hierarchical"), 6),
                                      rep(c("DBSCAN"), 6)))

pdf("iris_timings.pdf", width = 5, height = 4)
ggplot(data = time_data, aes(x = n_samples, y = time, color = algorithm,
                             group = algorithm)) +
  geom_line() +
  geom_point() + 
  labs(title = "Algorithm Run Time on Iris Data", x = "Number of samples",
       y = "Time in seconds", color = "Algorithm")
dev.off()


### Cell line clustering ###

cell_counts <- read.csv("Data/sc_10x_5cl.count.csv.gz")
cell_metadata <- read.csv("Data/sc_10x_5cl.metadata.csv.gz")

# Flip the data so rows are cells (samples) and columns are genes (features)
cell_data <- data.frame(t(cell_counts))
# Extract the cell line labels
cell_labels <- factor(cell_metadata$cell_line)

n_cell_samples <- 2000

if (n_cell_samples < 3918) {
  # Join labels and data
  cell_data$cell_line <- cell_labels
  # Randomly select n samples
  cell_data <- sample_n(cell_data, n_cell_samples)
  cell_labels <- cell_data$cell_line
  # Remove labels
  cell_data <- cell_data[ , -which(names(cell_data) %in% c("cell_line"))]
}

print(paste("Number of cells:", nrow(cell_data)))
print(paste("Number of genes:", ncol(cell_data)))
print(paste("Cell lines:", toString(unique(cell_labels))))

cell_pca <- reduceDimensions(cell_data, reduction = "pca", n_pcs = 2)
cell_tsne <- reduceDimensions(cell_data, reduction = "tsne")
cell_umap <- reduceDimensions(cell_data, reduction = "umap")

cell_reduction <- cell_tsne

cell_kmeans <- runClustering(cell_reduction, cell_labels, algorithm = "kmeans")
cell_hierarchical <- runClustering(cell_reduction, cell_labels, algorithm = "hierarchical")
cell_dbscan <- runClustering(cell_reduction, cell_labels, algorithm = "dbscan",
                             epsilon = 1.3, min_points = 30)

cell_plot_data <- data.frame(cell_reduction)
cell_plot_data$label <- cell_labels
cell_plot_data$kmeans <- cell_kmeans$clusters
cell_plot_data$hierarchical <- cell_hierarchical$clusters
cell_plot_data$dbscan <- cell_dbscan$clusters

# Plot clustering results
cell_kmeans_plot <- plotClusters(data = cell_plot_data,
                                 algorithm = "kmeans",
                                 colours = plot_colours[1:length(levels(cell_plot_data$kmeans))],
                                 silhouette_score = cell_kmeans$silhouette,
                                 ari = cell_kmeans$ari)

cell_hierarchical_plot <- plotClusters(data = cell_plot_data,
                                       algorithm = "hierarchical",
                                       colours = plot_colours[1:length(levels(cell_plot_data$hierarchical))],
                                       silhouette_score = cell_hierarchical$silhouette,
                                       ari = cell_hierarchical$ari)

cell_dbscan_plot <- plotClusters(data = cell_plot_data,
                                 algorithm = "dbscan",
                                 colours = plot_colours[1:length(levels(cell_plot_data$dbscan))],
                                 silhouette_score = cell_dbscan$silhouette,
                                 ari = cell_dbscan$ari)

pdf("cell_clusters_tsne.pdf", width = 10, height = 8)
grid.arrange(cell_label_plot_tsne, cell_kmeans_plot,
             cell_hierarchical_plot, cell_dbscan_plot, nrow = 2)
dev.off()

# Ground truth plot data
cell_plot_pca <- data.frame(cell_pca)
cell_plot_tsne <- data.frame(cell_tsne)
cell_plot_umap <- data.frame(cell_umap)
cell_plot_pca$label <- cell_labels
cell_plot_tsne$label <- cell_labels
cell_plot_umap$label <- cell_labels

# Plot ground truth
cell_label_plot_pca <- plotClusters(data = cell_plot_pca,
                                    algorithm = NA,
                                    colours = rev(plot_colours)[1:length(unique(cell_labels))],
                                    title = "Cell Line Ground Truth")
cell_label_plot_tsne <- plotClusters(data = cell_plot_tsne,
                                     algorithm = NA,
                                     colours = rev(plot_colours)[1:length(unique(cell_labels))],
                                     title = "Cell Line Ground Truth")
cell_label_plot_umap <- plotClusters(data = cell_plot_umap,
                                     algorithm = NA,
                                     colours = rev(plot_colours)[1:length(unique(cell_labels))],
                                     title = "Cell Line Ground Truth")

pdf("cell_labels.pdf", width = 10, height = 4)
grid.arrange(cell_label_plot_pca, cell_label_plot_umap, nrow = 1)
dev.off()

# Plot algorithm run times
time_data <- data.frame(n_samples = rep(c(125, 250, 500, 1000, 2000), 3),
                        time = c(0.01, 0.08, 0.14, 0.25, 0.95,
                                 0.1, 0.17, 1.07, 6.42, 51.59,
                                 4.45, 18.95, 79.56, 311.33, 1262.86),
                        algorithm = c(rep(c("k-means"), 5),
                                      rep(c("hierarchical"), 5),
                                      rep(c("DBSCAN"), 5)))

kmeans_time_data <- data.frame(n_samples = c(125, 250, 500, 1000, 2000),
                               time = c(0.01, 0.04, 0.12, 0.25, 0.95))
hierarchical_time_data <- data.frame(n_samples = c(125, 250, 500, 1000, 2000),
                               time = c(0.1, 0.17, 1.07, 6.42, 51.59))
dbscan_time_data <- data.frame(n_samples = c(125, 250, 500, 1000, 2000),
                               time = c(4.45, 18.95, 79.56, 311.33, 1262.86))

kmeans_time_plot <- ggplot(data = kmeans_time_data,
                           aes(x = n_samples,
                               y = time)) +
  geom_line(color = "red") +
  geom_point() + 
  labs(title = "K-Means Time", x = "Number of samples",
       y = "Time in seconds")
hierarchical_time_plot <- ggplot(data = hierarchical_time_data,
                                 aes(x = n_samples,
                                     y = time)) +
  geom_line(color = "green") +
  geom_point() + 
  labs(title = "Agglomerative Hierarchical Time", x = "Number of samples",
       y = "Time in seconds")
dbscan_time_plot <- ggplot(data = dbscan_time_data,
                           aes(x = n_samples,
                               y = time)) +
  geom_line(color = "blue") +
  geom_point() + 
  labs(title = "DBSCAN Time", x = "Number of samples",
       y = "Time in seconds")

pdf("cell_timings.pdf", width = 10, height = 4)
grid.arrange(kmeans_time_plot, hierarchical_time_plot,
             dbscan_time_plot, nrow = 1)
dev.off()