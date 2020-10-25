#' Precision of a co-expression network.
#'
#' This function computes the precision of a co-expression network at top edges
#' i.e., the fraction of true edges in a given number of top weighted edges.
#' It has options to set the number of top edges and the known score thresholds
#' to define which edges are true.
#'
#' @param net matrix or data.frame. A gene x gene matrix representing edge weights
#' between genes in a co-expression network. See details.
#'
#' @param known matrix or data.frame. A gene x gene matrix representing the probability
#' that edges between genes true. See details.
#'
#' @param na.ignore character representing how \code{NA}'s should be handled.
#' Accepted values are \code{'net'}, \code{'known'} and  \code{'any'}.
#' If \code{'net'}, edges with \code{NA} weight in \code{net} are ignored.
#' If \code{'known'}, edges with \code{NA} weight in \code{known} are ignored.
#' If \code{'any'}, edges with \code{NA} weight in either \code{net} or \code{known} are ignored.
#'
#' @param neg.treat character representing how negative values in \code{net} should be treated.
#' Accepted values are \code{'none'}, \code{'warn'} and  \code{'error'}.
#' If \code{'allow'}, negative values are allowed.
#' If \code{'warn'}, a warning is generated.
#' If \code{'error'}, an error is generated.
#'
#' @param score.thresholds numeric vector. Each value must be in the rage \[0,1\].
#' If known interaction score is equal to or greater than a score threshold,
#' the corresponding edge is considered true.
#'
#' @param n.top.edges vector with numeric values or \code{NA}. The number of
#' top weighted edges in \code{net} to use to compute the precision.
#' If \code{NA}, the number of edges with a known interaction score
#' equal to or greater than the given score threshold is used.
#' If the number of edges is greater than the number of available edges,
#' all available edges are used.
#' If multiple edges have the same weight as the last top-weighted edge,
#' all of those edges are used.
#'
#' @details
#' Each value in \code{known} must be in the range \[0, 1\] representing
#' the probability that the corresponding edge (interaction) is true.
#' While the values in \code{net} are not limited to any range,
#' each value should represent the relative probability that
#' the corresponding edge is true. In other words, larger values should
#' represent higher confidence in corresponding edges.
#' If the sign of values in \code{net} represents positive or negative
#' associations between genes, you probably should provide absolute values.
#' If you still want to allow negative values in \code{net},
#' you may set \code{neg.treat = "allow"}.
#' In this case, any negative value will represent lower confidence than
#' any non-negative value.
#'
#' Both \code{net} and \code{known} must be square matrices of same dimension.
#' Either both or none of the matrices should have row and column names.
#' If available, the set of row names and column names must be unique and same in each matrix.
#' The set of row and columns names of both matrices should also be same.
#' Both matrices must be symmetric when rows and columns are in the same order.
#' Diagonal entries in the matrices are ignored.
#'
#'
#' @return A matrix with \code{length(score.thresholds)} rows and
#' \code{length(n.top.edges)} columns.
#' Each value in the matrix represent the precision at a given number of top edges (column)
#' using a given score threshold (row).
#'
#' @export
#' @examples
#' genes = sprintf("G%d", 1:10)
#' dummy_net = matrix(rnorm(length(genes)^2), nrow = length(genes), dimnames = list(genes, genes))
#' dummy_net = abs((dummy_net + t(dummy_net))/2)                    # symmetric network
#' dummy_ppi = abs(dummy_net + rnorm(length(dummy_net)))
#' dummy_ppi = (dummy_ppi + t(dummy_ppi)) / (2 * max(dummy_ppi))    # symmetric ppi
#' net_precision = coexpression_known_interactions_precision_at_top(
#'   net = dummy_net,
#'   known = dummy_ppi,
#'   n.top.edges = c(10, 30, NA),
#'   score.thresholds = c(0.25, 0.5)
#' )
#' print(net_precision)

coexpression_known_interactions_precision_at_top <- function (net, known,
                                                              score.thresholds,
                                                              n.top.edges = NA,
                                                              na.ignore = "known",
                                                              neg.treat = "error")
{
  ### check arguments
  check_row_col_compatibility_of_net_and_known(net = net, known = known)
  check_value_range_of_net_and_known(net = net, known = known, neg.treat = neg.treat)
  if(!is.numeric(score.thresholds) || length(score.thresholds)==0)
    stop("score.thresholds must be numeric.")
  if(sum(score.thresholds<0) > 0 || sum(score.thresholds>1) > 0)
    stop("score.thresholds must be in the rage [0,1]")
  if(!all(is.na(n.top.edges) | (is.numeric(n.top.edges) & n.top.edges>=1)))
    stop("n.top.edges must contain numeric values >=1 or NA")

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

  ### check na.ignore compatibility
  check_na.ignore_compatibility_of_net_and_known(net = net, known = known, na.ignore = na.ignore)
  check_matrix_symmetry_of_net_and_known(net = net, known = known)

  # known counts by different thresholds
  ppi_score_values = get_lower_triangle_vector(known)
  net_weight_values = get_lower_triangle_vector(net)

  ### exclude na
  if(na.ignore == "known"){
    non_na_idx = which(!is.na(ppi_score_values))
  } else if(na.ignore == "net"){
    non_na_idx = which(!is.na(net_weight_values))
  } else if(na.ignore == "any"){
    non_na_idx = which((!is.na(ppi_score_values)) & (!is.na(net_weight_values)))
  }

  ppi_score_values = ppi_score_values[non_na_idx]
  net_weight_values = net_weight_values[non_na_idx]

  ### compute precision at top
  top_known_factions = sapply(n.top.edges, function(top){
    if(is.na(top)){
      sapply(score.thresholds, function(sth){
        # count the number of available true edges
        n_known = sum(ppi_score_values >= sth)
        top_weight_threshold = stats::quantile(net_weight_values, probs =  max(1 - n_known / length(net_weight_values), 0))
        top_weight_index = net_weight_values >= top_weight_threshold
        n_known_at_top = sum(ppi_score_values[top_weight_index] >= sth)
        frac_known_at_top = n_known_at_top / sum(top_weight_index)
        return(frac_known_at_top)
      })
    } else{
      top_weight_threshold = stats::quantile(net_weight_values, probs =  max(1 - top / length(net_weight_values), 0))
      top_weight_index = net_weight_values >= top_weight_threshold
      sapply(score.thresholds, function(sth){
        n_known_at_top = sum(ppi_score_values[top_weight_index] >= sth)
        frac_known_at_top = n_known_at_top / sum(top_weight_index)
        return(frac_known_at_top)
      })
    }
  })

  if(!is.matrix(top_known_factions)){
    top_known_factions = matrix(top_known_factions,
                                nrow = length(score.thresholds))
  }
  n.top.edges[is.na(n.top.edges)] = "NA"
  rownames(top_known_factions) = as.character(score.thresholds)
  colnames(top_known_factions) = as.character(n.top.edges)

  return(top_known_factions)
}
