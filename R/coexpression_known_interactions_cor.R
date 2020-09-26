#' Correlation between a co-expression network and a known interaction network.
#'
#' This function computes correlations between edge weights of a co-expression network
#' and corresponding known gene-gene interaction scores.
#'
#' @param net matrix or data.frame. A gene x gene matrix representing edge weights
#' between genes in a co-expression network. See details.
#'
#' @param known matrix or data.frame. A gene x gene matrix representing the probability
#' that edges between genes true. See details.
#'
#' @param method a character string indicating which correlation coefficient is to be
#' used for the test. One of \code{"spearman"}, \code{"pearson"}, or \code{"kendall"},
#' can be abbreviated.
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
#' @return A list with class \code{"htest"}, same as output from \code{\link[stats]{cor.test}}.
#'
#' @export
#' @examples
#' genes = sprintf("G%d", 1:10)
#' dummy_net = matrix(rnorm(length(genes)^2), nrow = length(genes), dimnames = list(genes, genes))
#' dummy_net = abs((dummy_net + t(dummy_net))/2)                    # symmetric undirected nework
#' dummy_ppi = abs(dummy_net + rnorm(length(dummy_net)))
#' dummy_ppi = (dummy_ppi + t(dummy_ppi)) / (2 * max(dummy_ppi))    # symmetric known interaction probability
#' cor_res = coexpression_known_interactions_cor(net = dummy_net, known = dummy_ppi, method = "spearman")
#' print(sprintf('Spearman roh: %g, p.value: %s', cor_res$estimate, cor_res$p.value))

coexpression_known_interactions_cor <- function (net, known,
                                                 method = "spearman",
                                                 na.ignore = "known",
                                                 neg.treat = "error")
{
  ### check arguments
  check_row_col_compatibility_of_net_and_known(net = net, known = known)
  check_value_range_of_net_and_known(net = net, known = known, neg.treat = neg.treat)

  ### convert net and known to a matrix
  if(!is.matrix(net))
    net = as.matrix(net)
  if(!is.matrix(known))
    known = as.matrix(known)

  ### ensure identical gene ordering in net and known
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

  ### correlation test
  known_score_values = get_lower_triangle_vector(known)
  net_weight_values = get_lower_triangle_vector(net)
  rm(list = c("net", "known"))

  ### exclude NA
  if(na.ignore == "known"){
    non_na_idx = which(!is.na(known_score_values))
  } else if(na.ignore == "net"){
    non_na_idx = which(!is.na(net_weight_values))
  } else if(na.ignore == "any"){
    non_na_idx = which((!is.na(known_score_values)) & (!is.na(net_weight_values)))
  }

  known_score_values = known_score_values[non_na_idx]
  net_weight_values = net_weight_values[non_na_idx]

  ### correlation test
  cor_test_res = cor.test(net_weight_values, known_score_values, method = method)
  return(cor_test_res)
}
