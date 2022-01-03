# Euclidean distance: https://www.analyticsvidhya.com/blog/2020/02/4-types-of-distance-metrics-in-machine-learning/
# https://towardsdatascience.com/the-5-clustering-algorithms-data-scientists-need-to-know-a36d136ef68

library(dplyr)
library(clue)
library(fossil)
library(ggplot2)
library(clValid)
library(umap)
library(Rtsne)
library(gridExtra)
library(cluster)

setwd("C:/Users/redds/Documents/GitHub/Iris-Clustering")

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
  
  # Assign cluster numbers to the input points
  clustered_points <- data.frame(points)
  clustered_points$cluster <- NA
  
  for (i in 1:length(clusters)) {
    for (j in 1:nrow(clusters[[i]]))
      clustered_points[rownames(clusters[[i]][j,]),]$cluster <- i
  }
  
  return(factor(clustered_points$cluster))
}

# Apply dimensionality reduction algorithm
reduceDimensions <- function(data, reduction = "pca", n_pcs = 0) {
  if (reduction == "umap") {
    # UMAP dimensionality reduction
    features = umap(data)$layout
    colnames(features) <- c("UMAP1", "UMAP2")
    
  } else if (reduction == "tsne") {
    # t-SNE dimensionality reduction
    features <- Rtsne(distinct(data), dims = 2, perplexity = 30, max_iter = 500)$Y
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

# Find the best mapping between labels and predicted clusters using the
# Hungarian matching algorithm
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
    hierarchical_results <- agglomerativeClustering(data)
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
  average_silhouette = mean(silhouette_scores[, 3])
  
  return(list(clusters = clusters, run_time = run_time,
              ari = adjusted_rand, silhouette = average_silhouette))
}

### Iris clustering ###

# Extract features and metadata
unique_iris <- distinct(iris)
iris_features <- unique_iris[, !names(unique_iris) %in% c("Species")]
iris_metadata <- unique_iris[, "Species"]

iris_pca <- reduceDimensions(iris_features, reduction = "pca", n_pcs = 2)
iris_tsne <- reduceDimensions(iris_features, reduction = "tsne")
iris_umap <- reduceDimensions(iris_features, reduction = "umap")

# Set the dimensionality reduction to use
iris_reduction <- as.matrix(iris_features[c(1,2)])
iris_reduction <- iris_umap

iris_kmeans <- runClustering(iris_reduction, iris_metadata, algorithm = "kmeans")
iris_hierarchical <- runClustering(iris_reduction, iris_metadata, algorithm = "hierarchical")
iris_dbscan <- runClustering(iris_reduction, iris_metadata, algorithm = "dbscan",
                             epsilon = 0.8, min_points = 3)

# Create graphs
iris_plot_data <- data.frame(iris_reduction)
iris_plot_data$label <- iris_metadata
iris_plot_data$kmeans <- iris_kmeans$clusters
iris_plot_data$hierarchical <- iris_hierarchical$clusters
iris_plot_data$dbscan <- iris_dbscan$clusters

cluster_colours <- c("#fc8021", "#462cc7", "#3ab03a", "#ff2181", "#385df2", 
                     "#db43fa", "#5ce1ed")

iris_label_plot <- ggplot(iris_plot_data,
                          aes(x = get(names(iris_plot_data)[1]),
                              y = get(names(iris_plot_data)[2]),
                              color = label)) +
  geom_point(size = 2) +
  scale_color_manual(values = rev(cluster_colours)[1:length(levels(iris_plot_data$kmeans))]) +
  labs(title = "Iris Ground Truth", color = 'Label',
       x = names(iris_plot_data[1]), y = names(iris_plot_data[2]))

iris_kmeans_plot <- ggplot(iris_plot_data,
                           aes(x = get(names(iris_plot_data)[1]),
                               y = get(names(iris_plot_data)[2]),
                               color = kmeans, shape = label)) +
  geom_point(size = 2) +
  scale_color_manual(values = cluster_colours[1:length(levels(iris_plot_data$kmeans))]) +
  labs(title = paste("K-Means \nSilhouette: ", signif(iris_kmeans$silhouette, 3),
                     ", ARI: ", signif(iris_kmeans$ari, 3), sep = ""),
       x = names(iris_plot_data[1]), y = names(iris_plot_data[2]),
       color = 'Cluster', shape = 'Ground Truth')

iris_hierarchical_plot <- ggplot(iris_plot_data, 
                                 aes(x = get(names(iris_plot_data)[1]),
                                     y = get(names(iris_plot_data)[2]),
                                     color = hierarchical, shape = label)) +
  geom_point(size = 2) +
  scale_color_manual(values = cluster_colours[1:length(levels(iris_plot_data$hierarchical))]) +
  labs(title = paste("Agglomerative Hierarchical \n",
                     "Silhouette: ", signif(iris_hierarchical$silhouette, 3),
                     ", ARI: ", signif(iris_hierarchical$ari, 3), sep = ""),
       x = names(iris_plot_data[1]), y = names(iris_plot_data[2]),
       color = 'Cluster', shape = 'Ground Truth')

iris_dbscan_plot <- ggplot(iris_plot_data,
                           aes(x = get(names(iris_plot_data)[1]),
                               y = get(names(iris_plot_data)[2]),
                               color = dbscan, shape = label)) +
  geom_point(size = 2) +
  scale_color_manual(values = cluster_colours[1:length(levels(iris_plot_data$dbscan))]) +
  labs(title = paste("DBSCAN \nSilhouette: ", signif(iris_dbscan$silhouette, 3),
       ", ARI: ", signif(iris_dbscan$ari, 3), sep = ""),
       x = names(iris_plot_data[1]), y = names(iris_plot_data[2]),
       color = 'Cluster', shape = 'Ground Truth')

grid.arrange(iris_label_plot, iris_kmeans_plot,
             iris_hierarchical_plot, iris_dbscan_plot, nrow = 2)

pdf("iris_labels_raw.pdf", width = 5, height = 4)
iris_label_plot
dev.off()

### Mouse clustering ###

# Read sci-RNA-seq3 dataset
gene_counts <- readRDS("gene_count_cleaned_sampled_100k.RDS", refhook = NULL)
gene_counts <- as.data.frame(summary(gene_counts))
names(gene_counts) <- c("gene_id", "cell_id", "count")

# Read the metadata
gene_annotation <- read.csv("gene_annotate.csv")
cell_annotation <- read.csv("GSE119945_cell_annotate.csv.gz")
cell_annotation <- cell_annotation[c("id", "sex", "day", "Main_Cluster")]

# Extract 1000 cells
n_cells = 1000
gene_counts = gene_counts[gene_counts$cell_id <= n_cells,]
cell_annotation = cell_annotation[1:n_cells,]

# Count the number of unique genes
n_genes = length(unique(gene_counts$gene_id))

# Create empty count matrix
count_df <- as.data.frame(matrix(0, ncol = n_cells, nrow = n_genes))
colnames(count_df) <- unique(gene_counts$cell_id)
rownames(count_df) <- unique(gene_counts$gene_id)

# Fill the count matrix 
for (cell_id in 1:n_cells) {
  cell_rows <- gene_counts[gene_counts$cell_id == cell_id,]
  
  for (row_idx in 1:length(cell_rows)) {
    # Record each gene count for the cell
    row <- cell_rows[row_idx,]
    count_df[row$gene_id, cell_id] <- row$count
  }
}

# Set the rows of the count data as Ensembl IDs
gene_names = c()

for (gene_id in rownames(count_df)) {
  gene_names = c(gene_names, gene_annotation[as.numeric(gene_id),]$gene_id)
}

rownames(count_df) <- gene_names
# Remove genes with all zero counts
count_df <- count_df[rowSums(count_df != 0) > 0, ]

mouse_pca <- reduceDimensions(t(count_df), reduction = "pca", n_pcs = 2)
mouse_tsne <- reduceDimensions(t(count_df), reduction = "tsne")
mouse_umap <- reduceDimensions(t(count_df), reduction = "umap")

s = data.frame(standardize(mouse_reduction))

mouse_reduction <- mouse_umap
mouse_metadata <- factor(cell_annotation$Main_Cluster)

mouse_kmeans <- runClustering(mouse_reduction, mouse_metadata, algorithm = "kmeans")

mouse_plot_data <- data.frame(mouse_reduction)
mouse_plot_data$label <- mouse_metadata
mouse_plot_data$kmeans <- mouse_kmeans$clusters
# mouse_plot_data$hierarchical <- iris_hierarchical$clusters
# mouse_plot_data$dbscan <- iris_dbscan$clusters

#mouse_plot_data <- mouse_plot_data[which(mouse_plot_data$PC2 < 1),]

mouse_label_plot <- ggplot(mouse_plot_data,
                          aes(x = get(names(mouse_plot_data)[1]),
                              y = get(names(mouse_plot_data)[2]),
                              color = label)) +
  geom_point(size = 2) +
  scale_color_manual(values = rev(cluster_colours)[1:length(levels(mouse_plot_data$kmeans))]) +
  labs(title = "Iris Ground Truth", color = 'Label',
       x = names(mouse_plot_data[1]), y = names(mouse_plot_data[2]))