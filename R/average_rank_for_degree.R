#' Average rank of edges corresponding to known gene-gene interactions of a specific degree.
#'
#' This function computes the average rank of edges in a reconstructed network where these 
#' edges correspond to a specific degree in a known gene-gene interaction network. Degree 1 
#' corresponds to first degree edges, while higher degrees correspond to indirect interactions
#' in the known interaction network.
#'
#' @param net matrix or data.frame. A gene x gene matrix representing edge weights
#' between genes in a co-expression network. See details.
#'
#' @param known matrix or data.frame. A gene x gene matrix representing the probability
#' that edges between genes true. See details.
#'
#' @param degrees list. A list of positive-valued degrees in the STRING network.
#'
#' @details
#' Each value in \code{known} must be in the range \[0, 1\], where
#' a nonzero value indicates the presence of an edge.
#' While the values in \code{net} are not limited to any range,
#' each value should represent the relative probability that
#' the corresponding edge is true. In other words, larger values should
#' represent higher confidence in corresponding edges.
#' Values in \code{net} must be unsigned.
#'
#' Both \code{net} and \code{known} must be square matrices of same dimension.
#' Either both or none of the matrices should have row and column names.
#' If available, the set of row names and column names must be unique and same in each matrix.
#' The set of row and columns names of both matrices should also be same.
#' Both matrices must be symmetric when rows and columns are in the same order.
#' Diagonal entries in the matrices are ignored.
#'
#' @return A list object with the following items.
#' \item{degrees}{The degrees input as \code{degrees}}.
#' \item{average_ranks}{A list object containing the average ranks corresponding to the degree values 
#' in \code{degrees}.}
#'
#' @export
#' @examples
#' genes = sprintf("G%d", 1:10)
#' dummy_net = matrix(rnorm(length(genes)^2), nrow = length(genes), dimnames = list(genes, genes))
#' dummy_net = abs((dummy_net + t(dummy_net))/2)                  # symmetric network
#' dummy_ppi = matrix(rbinom(length(genes)^2, 1 , 0.5), nrow = length(genes), dimnames = list(genes, genes))
#' dummy_ppi[lower.tri(dummy_ppi)] = t(dummy_ppi)[lower.tri(dummy_ppi)]  # symmetric ppi
#' degrees = c(1, 2, 3)
#' average_ranks = average_rank_for_degree(
#'   net = dummy_net,
#'   known = dummy_ppi,
#'   degrees = degrees
#' )
#' plot(average_ranks$degrees, average_ranks$average_ranks)

requireNamespace('igraph', quietly = T)

find_shortest_paths = function(network){
  g = igraph::graph_from_adjacency_matrix(network, mode = c("undirected"), diag = FALSE)
  paths = igraph::shortest.paths(g)
  return(paths)
}

find_motif = function(paths, dist){
  return(paths == dist)
}

get_edge_ranks = function(network){ # Note: highest value corresponds to highest rank
  network[upper.tri(network, diag = TRUE)] = NA
  
  output = matrix(rank(-network, na.last = "keep", ties.method = "average"), nrow = nrow(network), ncol = ncol(network))
  rownames(output) = rownames(network)
  colnames(output) = colnames(network)
  output[is.na(output)] = 0
  return(output)
}

determine_avg_rank = function(ranked_network, edges){
  edges_to_use = edges * ranked_network 
  num_edges = sum(edges)/2
  edge_sum = sum(edges_to_use)
  return(edge_sum/num_edges)
}

average_rank_for_degree <- function(net, known, degrees){
  
  ### check arguments
  check_row_col_compatibility_of_net_and_known(net = net, known = known)
  
  neg.treat = "error"
  check_value_range_of_net_and_known(net = net, known = known, neg.treat = neg.treat)

  ### ensure known is binary
  known[known != 0] = 1
  
  ### convert net and known to a matrix
  if(!is.matrix(net))
    net = as.matrix(net)
  if(!is.matrix(known))
    known = as.matrix(known)

  ### ensure identical gene ordering in both net and known
  if(length(rownames(known)) > 0){
    genes = rownames(known)
    if(any(colnames(known) != genes ))  # check to avoid using extra memory
      known = known[genes, genes]
    if(any(rownames(net) != genes) || any(colnames(net) != genes))
      net = net[genes, genes]
  }

  ### check matrix symmetry
  check_matrix_symmetry_of_net_and_known(net = net, known = known)

  ### compute average ranks for degrees
  average_ranks = c()
  ranked_net = get_edge_ranks(net)
  
  paths = find_shortest_paths(known)
  for (degree in degrees){
    value = find_motif(paths, degree)
    average_rank = determine_avg_rank(ranked_net, value)
    average_ranks = append(average_ranks, average_rank)
  }

  return(list(degrees = degrees, average_ranks = average_ranks))
}


