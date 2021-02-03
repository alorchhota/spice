#' Mean rank of edges corresponding to known gene-gene interactions of a specific distance.
#'
#' This function computes the mean rank of edges in a reconstructed network, 
#' where each edge corresponds to two genes with a given geodesic distance 
#' in the known gene-gene interaction network. 
#' A distance of 1 corresponds to direct interactions, 
#' while higher distances correspond to indirect interactions 
#' in the known interaction network.
#'
#' @param net matrix or data.frame. A gene x gene matrix representing edge weights
#' between genes in a co-expression network. See details.
#'
#' @param known matrix or data.frame. A gene x gene matrix representing the probability
#' that edges between genes true. See details.
#'
#' @param d integer vector. A list of positive-valued geodesic distances in the known interaction network.
#'
#' @param score.threshold numeric. The value must be in the rage \[0,1\].
#' If known interaction score is greater than or equal to \code{score.threshold},
#' the corresponding edge is considered true.
#'
#' @details
#' Each value in \code{known} must be in the range \[0, 1\], where
#' a value greater than or equal to \code{score.threshold} indicates 
#' the presence of an edge. 
#' A \code{NA} value will be treated as an edge with a zero value.
#' 
#' Values in \code{net} must be unsigned.
#' Each value should represent the relative probability that
#' the corresponding edge is true. In other words, larger values should
#' represent higher confidence in corresponding edges. 
#' A \code{NA} value will be treated as an edge with a \code{-Inf} value, 
#' i.e. an edge with very low confidence.
#'
#' Both \code{net} and \code{known} must be square matrices of same dimension.
#' Either both or none of the matrices should have row and column names.
#' If available, the set of row names and column names must be unique and same in each matrix.
#' The set of row and columns names of both matrices should also be same.
#' Both matrices must be symmetric when rows and columns are in the same order.
#' Diagonal entries in the matrices are ignored.
#'
#' @return A numeric vector of the same length as \code{d} containing 
#' the mean ranks corresponding to the distance values in \code{d}.
#'
#' @export
#' @examples
#' genes = sprintf("G%d", 1:10)
#' dummy_net = matrix(rnorm(length(genes)^2), nrow = length(genes), dimnames = list(genes, genes))
#' dummy_net = abs((dummy_net + t(dummy_net))/2)                  # symmetric network
#' dummy_ppi = matrix(rbinom(length(genes)^2, 1 , 0.5), nrow = length(genes), dimnames = list(genes, genes))
#' dummy_ppi[lower.tri(dummy_ppi)] = t(dummy_ppi)[lower.tri(dummy_ppi)]  # symmetric ppi
#' d = 1:3
#' mean_ranks = mean_rank_for_known_geodesic_distance(net = dummy_net,
#'                                                       known = dummy_ppi,
#'                                                       d = d)
#' barplot(
#'   mean_ranks,
#'   main = "Mean ranks of edges in learned network",
#'   xlab = "Known interaction distance",
#'   ylab = "Mean Rank",
#'   names.arg = d
#' )
#' barplot(
#'   c(
#'     mean_ranks[1] - mean_ranks[2],
#'     mean_ranks[2] - mean_ranks[3]
#'   ),
#'   main = "Differences between mean ranks of two known distances",
#'   xlab = "Known geodesic distances",
#'   ylab = "Difference in mean rank",
#'   names.arg = c("1-2", "2-3")
#' )

mean_rank_for_known_geodesic_distance <- function(net,
                                                  known,
                                                  d = 1:3,
                                                  score.threshold = 1e-10) {
  requireNamespace('igraph', quietly = T)
  requireNamespace('data.table', quietly = T)
  
  ### check arguments
  check_row_col_compatibility_of_net_and_known(net = net, known = known)
  neg.treat = "error"
  check_value_range_of_net_and_known(net = net, known = known, neg.treat = neg.treat)
  
  if (!is.numeric(d) || d %% 1 != 0 || any(d < 1)) {
    stop("d must be a vector of positive integers.")
  }
  
  if(!is.numeric(score.threshold) ||
     length(score.threshold) != 1 ||
     score.threshold < 0 ||
     score.threshold > 1) {
    stop("score.threshold must be a number between 0 and 1.")
  } 
  
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
  
  ### handle NA values
  net[is.na(net)] = -Inf
  known[is.na(known)] = 0
  
  ### check matrix symmetry
  check_matrix_symmetry_of_net_and_known(net = net, known = known)
  
  ### ensure known is binary
  known[known >= score.threshold] = 1
  known[known < score.threshold] = 0
  
  ### compute mean ranks for distances
  ranked_net = get_edge_ranks(net)
  paths = find_shortest_paths(known)
  mean_ranks = sapply(d, function(distance){
    motif = find_motif(paths, distance)
    mean_rank = determine_mean_rank(ranked_net, motif)
  })
  
  return(mean_ranks)
}

find_shortest_paths = function(net){
  g = igraph::graph_from_adjacency_matrix(net, mode = c("undirected"), diag = FALSE)
  paths = igraph::distances(g)
  return(paths)
}

find_motif = function(paths, distance){
  return(paths == distance)
}

get_edge_ranks = function(net){
  net[upper.tri(net, diag = TRUE)] = NA
  
  output = matrix(
    data.table::frank(as.numeric(-net), na.last = "keep", ties.method = "average"),
    nrow = nrow(net),
    ncol = ncol(net)
  )
  output = pmax(output, t(output), na.rm = T)
  rownames(output) = rownames(net)
  colnames(output) = colnames(net)
  return(output)
}

determine_mean_rank = function(ranked_net, edges){
  return(mean(ranked_net[edges]))
}
