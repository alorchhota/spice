#' Evaluate co-expression network's ability to find edges where both genes of an edge share a pathway.
#'
#' If both genes of an edge is available in at least one pathway (gene set),
#' the edge is considered positive.
#' Using the edge weights as the score for the positive class,
#' this function computes both the area under the precision-recall curve
#' and the area under the ROC curve.
#'
#' @param net matrix or data.frame. A gene x gene matrix representing edge weights
#' between genes in a co-expression network.
#' Gene names must be available as row and column names. See details.
#'
#' @param pathways list. List of pathways where each entry contains the genes in each pathway.
#'
#' @param neg.treat character representing how negative values in \code{net} should be treated.
#' Accepted values are \code{'none'}, \code{'warn'} and  \code{'error'}.
#' If \code{'allow'}, negative values are allowed.
#' If \code{'warn'}, a warning is generated.
#' If \code{'error'}, an error is generated.
#'
#' @param curve logical. Should the curves, required for plotting, be generated?
#' @param max.compute logical. Should the maximum area under the curve be computed?
#' @param min.compute logical. Should the minimum area under the curve be computed?
#' @param rand.compute logical. Should the are under the curve for a random classifier be computed?
#' @param dg.compute logical. Should the area under the precision-recall curve
#' according to the interpolation of Davis and Goadrich be computed?
#'
#' @details
#' Each value in \code{net} should represent the relative probability that
#' the corresponding edge is true. In other words, larger values should
#' represent higher confidence in corresponding edges.
#' If the sign of values in \code{net} represents positive or negative
#' associations between genes, you probably should provide absolute values.
#' If you still want to allow negative values in \code{net},
#' you may set \code{neg.treat = "allow"}.
#' In this case, any negative value will represent lower confidence than
#' any non-negative value.
#'
#' \code{net} must be a square matrix.
#' Gene names must be available as row and column names.
#' Gene names must be unique.
#' \code{net} must be symmetric when rows and columns are identically ordered.
#' Diagonal entries are ignored.
#'
#' @return A list object with the following items.
#' \item{pr}{Precision-recall curve object. See \code{\link[PRROC]{pr.curve}} for details.}
#' \item{roc}{ROC curve object. See \code{\link[PRROC]{roc.curve}} for details.}
#'
#' @export
#' @examples
#' require(msigdbr)
#' msigdb_df = as.data.frame(msigdbr(species = "Homo sapiens", category = "H"))
#' pathways = tapply(msigdb_df$human_gene_symbol, msigdb_df$gs_id, FUN = c)
#' genes = c('TP53', 'RBM3', 'SF3', 'LIM12', 'ATM', 'TMEM160', 'BCL2L1', 'MDM2', 'PDR', 'MEG3', 'EGFR', 'CD96', 'KEAP1', 'SRSF1', 'TSEN2')
#' dummy_net = matrix(rnorm(length(genes)^2), nrow = length(genes), dimnames = list(genes, genes))
#' dummy_net = abs((dummy_net + t(dummy_net))/2)                    # symmetric undirected nework
#' auc_res = coexpression_shared_pathway_auc(net = dummy_net,
#'                                           pathways = pathways,
#'                                           curve = T)
#' print(sprintf('AUC under the precision-recall curve: %s', auc_res$pr$auc.integral))
#' print(sprintf('AUC under the ROC curve: %s', auc_res$roc$auc))
#' plot(auc_res$pr)
#' plot(auc_res$roc)

coexpression_shared_pathway_auc <- function(net, pathways,
                                            curve = F,
                                            max.compute=F,
                                            min.compute=F,
                                            rand.compute=F,
                                            dg.compute=F,
                                            na.rm = F,
                                            neg.treat = "error"){

  suppressPackageStartupMessages(require('PRROC'))

  ### check arguments
  net = get_checked_coexpression_network(net = net, varname = "net", check.names = T, check.symmetry = T, check.na = !na.rm)
  check_negative_values_of_net(net = net, varname = "net", neg.treat = neg.treat)
  if(!is.list(pathways))
    stop("pathways must be a list.")

  ### create shared pathway matrix
  shared_pathway_mat = matrix(0, nrow = nrow(net), ncol = ncol(net), dimnames = list(rownames(net), colnames(net)))
  bg = rownames(net)
  tmp <- lapply(pathways, function(genes){
    common_genes = intersect(genes, bg)
    if(length(common_genes) < 2)
      return(NULL)
    shared_pathway_mat[common_genes, common_genes] <<- 1
    return(NULL)
  })

  ### get vectors
  net_vec = get_lower_triangle_vector(net)
  rm("net")
  shared_pathway_vec = get_lower_triangle_vector(shared_pathway_mat)
  rm("shared_pathway_mat")

  ### remove negative values
  if(na.rm){
    non_na_idx = !is.na(net_vec)
    net_vec = net_vec[non_na_idx]
    shared_pathway_vec = shared_pathway_vec[non_na_idx]
    rm("non_na_idx")
  }

  ### compute AUC
  probj = pr.curve(scores.class0 = net_vec,
                   weights.class0 = shared_pathway_vec,
                   curve = curve,
                   max.compute = max.compute,
                   min.compute = min.compute,
                   rand.compute = rand.compute,
                   dg.compute = dg.compute)
  rocobj = roc.curve(scores.class0 = net_vec,
                     weights.class0 = shared_pathway_vec,
                     curve = curve,
                     max.compute = max.compute,
                     min.compute = min.compute,
                     rand.compute = rand.compute)

  return(list(pr = probj, roc = rocobj))
}
