#' Area under the curve of a co-expression network to find known gene-gene interactions.
#'
#' This function computes the area under the precision-recall curve and the ROC curve.
#' The area is computed by the \code{PRROC} package using known gene-gene interaction scores.
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
#' @param curve logical. Should the curves be returned?
#' @param max.compute logical. Should the maximum area under the curve be computed?
#' @param min.compute logical. Should the minimum area under the curve be computed?
#' @param rand.compute logical. Should the are under the curve for a random classifier be computed?
#' @param dg.compute logical. Should the area under the precision-recall curve
#' according to the interpolation of Davis and Goadrich be computed?
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
#' @return A list object with the following items.
#' \item{pr}{Precision-recall curve object. See \code{\link[PRROC]{pr.curve}} for details.}
#' \item{roc}{ROC curve object. See \code{\link[PRROC]{roc.curve}} for details.}
#'
#' @export
#' @examples
#' genes = sprintf("G%d", 1:10)
#' dummy_net = matrix(rnorm(length(genes)^2), nrow = length(genes), dimnames = list(genes, genes))
#' dummy_net = abs((dummy_net + t(dummy_net))/2)                    # symmetric undirected nework
#' dummy_ppi = abs(dummy_net + rnorm(length(dummy_net)))
#' dummy_ppi = (dummy_ppi + t(dummy_ppi)) / (2 * max(dummy_ppi))    # symmetric known interaction probability
#' auc_res = coexpression_known_interactions_auc(net = dummy_net, known = dummy_ppi, curve = T, max.compute = T, min.compute = T)
#' print(sprintf('AUC under the precision-recall curve: %s', auc_res$pr$auc.integral))
#' print(sprintf('AUC under the ROC curve: %s', auc_res$roc$auc))
#' plot(auc_res$pr, max.plot = TRUE, min.plot = TRUE, fill.area = TRUE)
#' plot(auc_res$roc, max.plot = TRUE, min.plot = TRUE, fill.area = TRUE)

coexpression_known_interactions_auc <- function(net, known,
                                                curve = F,
                                                max.compute=F,
                                                min.compute=F,
                                                rand.compute=F,
                                                dg.compute=F,
                                                na.ignore = "known",
                                                neg.treat = "error"){
  suppressPackageStartupMessages(require('PRROC'))

  ### check arguments
  check_row_col_compatibility_of_net_and_known(net = net, known = known)
  check_value_range_of_net_and_known(net = net, known = known, neg.treat = neg.treat)

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

  ### exclude NA
  if(na.ignore == "known"){
    na_idx = which(is.na(known))
  } else if(na.ignore == "net"){
    na_idx = which(is.na(net))
  } else if(na.ignore == "any"){
    na_idx = which(is.na(known) | is.na(net))
  }

  net[na_idx] = NA
  known[na_idx] = NA

  subset_idx = lower.tri(net, diag = F) & !is.na(known)
  probj = pr.curve(scores.class0 = net[subset_idx],
                   weights.class0 = known[subset_idx],
                   curve = curve,
                   max.compute = max.compute,
                   min.compute = min.compute,
                   rand.compute = rand.compute,
                   dg.compute = dg.compute)
  rocobj = roc.curve(scores.class0 = net[subset_idx],
                     weights.class0 = known[subset_idx],
                     curve = curve,
                     max.compute = max.compute,
                     min.compute = min.compute,
                     rand.compute = rand.compute)

  return(list(pr = probj, roc = rocobj))
}
