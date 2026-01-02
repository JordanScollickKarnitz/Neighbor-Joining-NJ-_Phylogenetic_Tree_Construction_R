############################################################
## Neighbor-Joining (NJ) Phylogenetic Tree Construction
## Distance metric: Jukes–Cantor
############################################################

## Example DNA sequences
seq1 <- "ACTGGGCT"
seq2 <- "CGTGAGCT"

############################################################
## Jukes–Cantor distance function
############################################################
jukes_cantor_distance <- function(seq1, seq2) {
  s1 <- strsplit(seq1, "")[[1]]
  s2 <- strsplit(seq2, "")[[1]]
  
  # Proportion of differing sites
  p <- sum(s1 != s2) / length(s1)
  
  # Jukes–Cantor correction
  d <- (-3/4) * log(1 - (4 * p / 3))
  return(d)
}

############################################################
## Build initial distance matrix
############################################################
allseq <- list(seq1 = seq1, seq2 = seq2)
n <- length(allseq)
seq_name <- names(allseq)

distance_matrix <- matrix(0, nrow = n, ncol = n)
rownames(distance_matrix) <- seq_name
colnames(distance_matrix) <- seq_name

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      distance_matrix[i, j] <- jukes_cantor_distance(allseq[[i]], allseq[[j]])
    }
  }
}

cat("Initial Distance Matrix:\n")
print(distance_matrix)

############################################################
## Neighbor-Joining setup
############################################################
active_clusters <- seq_name

tree_edges <- data.frame(
  parent = character(),
  child1 = character(),
  child2 = character(),
  branch1 = numeric(),
  branch2 = numeric(),
  stringsAsFactors = FALSE
)

next_node_id <- n + 1
current_matrix <- distance_matrix

############################################################
## Compute the Q-matrix
############################################################
calculate_Q_matrix <- function(dist_matrix) {
  n <- nrow(dist_matrix)
  Q <- matrix(0, n, n)
  rownames(Q) <- rownames(dist_matrix)
  colnames(Q) <- colnames(dist_matrix)
  
  row_sums <- rowSums(dist_matrix)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      Q[i, j] <- (n - 2) * dist_matrix[i, j] - row_sums[i] - row_sums[j]
      Q[j, i] <- Q[i, j]
    }
  }
  
  diag(Q) <- Inf
  return(Q)
}

############################################################
## Identify minimum entry in Q-matrix
############################################################
find_min_Q <- function(Q_matrix) {
  Q_matrix[lower.tri(Q_matrix)] <- Inf
  min_value <- min(Q_matrix)
  min_location <- which(Q_matrix == min_value, arr.ind = TRUE)[1, ]
  
  list(
    min_value = min_value,
    row_idx = min_location[1],
    col_idx = min_location[2]
  )
}

############################################################
## Compute branch lengths for joined clusters
############################################################
calculate_branch_lengths <- function(dist_matrix, cluster1, cluster2) {
  n <- nrow(dist_matrix)
  d_ij <- dist_matrix[cluster1, cluster2]
  
  sum_i <- sum(dist_matrix[cluster1, ]) - d_ij
  sum_j <- sum(dist_matrix[cluster2, ]) - d_ij
  
  branch_i <- (d_ij / 2) + ((sum_i - sum_j) / (2 * (n - 2)))
  branch_j <- d_ij - branch_i
  
  list(branch1 = branch_i, branch2 = branch_j)
}

############################################################
## Update distance matrix after joining two clusters
############################################################
update_matrix_nj <- function(matrix, cluster1, cluster2, new_cluster_name) {
  idx1 <- match(cluster1, colnames(matrix))
  idx2 <- match(cluster2, colnames(matrix))
  
  other_clusters <- colnames(matrix)[-c(idx1, idx2)]
  d_ij <- matrix[cluster1, cluster2]
  
  new_distances <- numeric(length(other_clusters))
  for (i in seq_along(other_clusters)) {
    k <- other_clusters[i]
    new_distances[i] <- (matrix[cluster1, k] + matrix[cluster2, k] - d_ij) / 2
  }
  
  updated_matrix <- matrix[other_clusters, other_clusters, drop = FALSE]
  updated_matrix <- rbind(updated_matrix, new_distances)
  rownames(updated_matrix)[nrow(updated_matrix)] <- new_cluster_name
  
  new_col <- c(new_distances, 0)
  updated_matrix <- cbind(updated_matrix, new_col)
  colnames(updated_matrix)[ncol(updated_matrix)] <- new_cluster_name
  
  return(updated_matrix)
}

############################################################
## Main Neighbor-Joining loop
############################################################
while (length(active_clusters) > 2) {
  Q_matrix <- calculate_Q_matrix(current_matrix)
  min_result <- find_min_Q(Q_matrix)
  
  cluster1 <- rownames(current_matrix)[min_result$row_idx]
  cluster2 <- colnames(current_matrix)[min_result$col_idx]
  
  branches <- calculate_branch_lengths(current_matrix, cluster1, cluster2)
  
  new_node <- paste0("Node", next_node_id)
  next_node_id <- next_node_id + 1
  
  tree_edges <- rbind(
    t
    