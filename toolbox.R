MAX_ATTEMPTS <- 100

from_incidence_matrix_to_bipartite <- function(tmp){
  g <- graph_from_incidence_matrix(tmp)
  cl <- clusters(g, "strong")
  g2 <- induced.subgraph(g, which(cl$membership == which.max(cl$csize)))
  return(g2)
}

## a more consistent sampling function
my.sample <- function(x, size=1) x[sample(length(x), size)]

# create a random walk through the complete graph
# add a link every time a new node is discovered
# this will provide a skeleton for the graph: 
# a (random) spanning tree will connect all of the nodes
random_spanning_tree <- function(num_rows, num_cols) {
  tree_B <- matrix(0, num_rows, num_cols)
  discovered_rows <- rep(0, num_rows)
  discovered_cols <- rep(0, num_cols)
  my_row <- sample(1:num_rows, 1)
  my_col <- sample(1:num_cols, 1)
  discovered_rows[my_row] <- 1
  discovered_cols[my_col] <- 1
  tree_B[my_row, my_col] <- 1
  while(sum(c(discovered_cols, discovered_rows) == 0) > 0){
    my_row <- sample(1:num_rows, 1)
    if (discovered_rows[my_row] == 0){
      discovered_rows[my_row] <- 1
      tree_B[my_row, my_col] <- 1 
    }
    my_col <- sample(1:num_cols, 1)
    if (discovered_cols[my_col] == 0){
      discovered_cols[my_col] <- 1
      tree_B[my_row, my_col] <- 1 
    }
  }
  return(tree_B)
}

sparse_erdos_renyi <- function(B) {
  ## create a random spanning tree to ensure connectedness
  ER_B <- random_spanning_tree(nrow(B), ncol(B))
  ## now add remaining links to make number of edges match original network
  ER_B <- as.vector(ER_B)
  ER_B[sample(which(ER_B == 0), sum(B) - sum(ER_B), replace=FALSE)] <- 1
  ER_B <- matrix(ER_B, nrow(B), ncol(B))
  dimnames(ER_B) <- dimnames(B)
  ## convert back to igraph object
  return(graph_from_incidence_matrix(ER_B))
}

erdos_renyi <- function(g) {
  ## get B from the igraph object
  B <- get.incidence(g)
  ## take some measures to ensure connectedness:
  connected <- FALSE
  n_attempts <- 0
  while (!connected) {
    if (n_attempts > MAX_ATTEMPTS) {
      cat("\nReached", MAX_ATTEMPTS, "ER randomizations without finding a connected one")
      return(NULL)
    }
    n_attempts <- n_attempts + 1
    if (ecount(g) > max(dim(B)) * log(vcount(g))) {
      ## then the ER graph is almost surely connected
      ER_g <- sample_bipartite(nrow(B), ncol(B), m=ecount(g), type="gnm")
    } else {
      ## ensure no empty rows before assigning remaining links
      ER_g <- sparse_erdos_renyi(B)
    }
    ## double check that it is connected
    connected <- is.connected(ER_g)
  }
  #cat("     ER")
  return(ER_g)
}

compute_stats <- function(B){
  g_connected <- from_incidence_matrix_to_bipartite(B)
  B <- get.incidence(g_connected)
  
  n_rows <- nrow(B)
  n_cols <- ncol(B)
  n_links <- sum(B)
  size <- n_rows * n_cols
  if (n_rows > n_cols) B <- t(B)
  # compute eigenvalues
  M <- B %*% t(B)
  eM <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  l1 <- sqrt(eM[1])
  l2 <- sqrt(eM[2])
  # compute expectations
  d_row <- rowSums(B)
  d_col <- colSums(B)
  pr <- (n_links - n_rows) /  ((n_cols - 1) * n_rows)
  pc <- (n_links - n_cols) / ((n_rows - 1) * n_cols)
  # expected l1 according to ER model
  l1_er <- sqrt((1 + (n_cols - 1) * pr * (3 + (n_cols - 2) * pr)) / (1 + (n_cols - 1) * pr) *
                  (1 + (n_rows - 1) * pc * (3 + (n_rows - 2) * pc)) / (1 + (n_rows - 1) * pc))
  # expected l1 according to CM model
  l1_cm <- sqrt((mean(d_col^2) / mean(d_col)) * (mean(d_row^2) / mean(d_row)))
  # expected l2 according to ER model (Marchenko-Pastur approximation)
  C <- (sum(d_row) - l1_cm^2) / size
  l2_mp <- sqrt(n_rows * C * (1 + sqrt(n_cols/n_rows))^2)
  return(tibble(l1 = l1, l2 = l2, l1_er = l1_er, l1_cm = l1_cm, l2_mp = l2_mp, 
                v1 = 1 - l1_er / l1,
                v2 = 1 - l1_cm / l1,
                v3 = 1 - l2_mp / l2))
}

compute_stats_all <- function(ID, Nrand = 1) {
  B <- load_web(here("network_data", ID))
  g_connected <- from_incidence_matrix_to_bipartite(B)
  1:Nrand %>% 
    map_dfr(~compute_stats(get.incidence(erdos_renyi(g_connected)))) %>% 
    mutate(data_type = 'er') %>% 
    bind_rows(
      mutate(compute_stats(get.incidence(g_connected)), data_type = 'ori')
    )
}