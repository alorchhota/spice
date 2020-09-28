#' Enrichment of pathways (gene sets) in a co-expression network.
#'
#' This function computes enrichment p-values of pathways in a co-expression network.
#' See details.
#'
#' @param net matrix or data.frame. A gene x gene matrix representing edge weights
#' between genes in a co-expression network.
#' Gene names must be available as row and column names. See details.
#'
#' @param pathways list. List of pathways where each entry contains the genes in each pathway.
#' Pathway names may be provided as \code{names(pathways)}.
#' If provided, pathway names must be unique.
#'
#' @param neg.treat character representing how negative values in \code{net} should be treated.
#' Accepted values are \code{'none'}, \code{'warn'} and  \code{'error'}.
#' If \code{'allow'}, negative values are allowed.
#' If \code{'warn'}, a warning is generated.
#' If \code{'error'}, an error is generated.
#'
#' @param na.rm logical. Should edges with \code{NA} weights be excluded?
#' If FALSE, \code{net} cannot have any edge with \code{NA} weight.
#'
#' @param min.gene integer. Each pathway must have at least \code{min.gene} genes
#' in \code{net}. Otherwise enrichment is not computed for the pathway.
#'
#' @param max.gene integer. Each pathway must have at most \code{max.gene} genes
#' in \code{net}. Otherwise enrichment is not computed for the pathway.
#'
#' @param iter integer. The number of random iterations or the number of random
#' gene sets to compute the null distribution.
#'
#' @param seed integer or NULL. Random number generator seed.
#'
#' @details To compute the enrichment p-value of a pathway in a co-expression network,
#' we define the score of a pathway as the the sum of weights in \code{net}
#' of all possible edges between the genes in the pathway. The enrichment p-value is
#' then defined as the probability that the score of the pathway is at least as big as
#' a random gene set with the same number of genes.
#'
#' To get the null distribution, we generate a number of (\code{iter}) random gene sets
#' where each gene set consists of the same number of randomly selected genes, compute
#' their scores, and fit a normal distribution.
#'
#' Enrichment of a pathway is computed only if at least \code{min.gene} and
#' at most \code{max.gene} genes from the pathway are available in \code{net}.
#' This criteria helps to avoid too small or to large pathways.
#'
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
#' @return A \code{data.frame} with the following columns.
#' \item{pathway}{Pathway name taken from \code{names(pathways)}.
#'                If pathway names are not provided, pathway index is used.}
#' \item{n.gene}{Number of genes from the pathway available in \code{net}.}
#' \item{p}{p-value for the pathway (computed using a fitted normal distribution as null).}
#' \item{p.empirical}{Empirical p-value for the pathway.}
#' \item{fdr}{False discovery rate computed using Benjamini-Hochberg method.}
#' \item{fdr.empirical}{Empirical false discovery rate computed using Benjamini-Hochberg method.}
#'
#' @export
#' @examples
#' genes = c('TP53', 'RBM3', 'SF3', 'LIM12', 'ATM', 'TMEM160', 'BCL2L1', 'MDM2', 'PDR', 'MEG3', 'EGFR', 'CD96', 'KEAP1', 'SRSF1', 'TSEN2')
#' dummy_net = matrix(rnorm(length(genes)^2), nrow = length(genes), dimnames = list(genes, genes))
#' dummy_net = abs((dummy_net + t(dummy_net))/2)                    # symmetric undirected nework
#' dummy_pathways = list(pathway1=c('TP53', 'RBM3', 'SF1', 'SF5'),
#'                        pathway2=c('LIM12', 'MDM2', 'BCL2L1', 'TMEM160', 'ATM'),
#'                        pathway3=c('EGFR', 'TP53', 'CD96', 'SRSF1', 'RBM14'))
#' enrich_res = coexpression_pathway_enrichment(net = dummy_net,
#'                                                 pathways = dummy_pathways,
#'                                                 min.gene = 3)
#' print(enrich_res)
#' n_sig = sum(enrich_res$fdr <= 0.05)
#' print(sprintf('Number of significantly enriched pathways: %d', n_sig))

coexpression_pathway_enrichment <- function(net, pathways,
                                            min.gene=5,
                                            max.gene=100,
                                            iter=10000,
                                            seed = NULL,
                                            na.rm = F,
                                            neg.treat = "error"){
  suppressPackageStartupMessages(require(MASS))
  net = get_checked_coexpression_network(net = net, varname = "net", check.names = T, check.symmetry = T, check.na = !na.rm)
  check_negative_values_of_net(net = net, varname = "net", neg.treat = neg.treat)
  check_pathways(pathways = pathways, varname = "pathways")
  if(length(names(pathways))==0)
    names(pathways) = as.character(seq_len(length(pathways)))

  if(sum(!is.finite(net)) > 0)
    stop('every edge weight in the network must be finite.')

  ### diagonal weights should be 0 (needed to compute sum-of-weights)
  diag(net) = 0

  ### set seed
  if(is.numeric(seed))
    set.seed(seed)

  ### process given gene sets to contain a valid number of background genes only
  bg_genes = rownames(net)
  pathways = lapply(pathways, intersect, y=bg_genes)
  pathways_lengths = sapply(pathways, length)
  pathways= pathways[pathways_lengths>=min.gene & pathways_lengths <= max.gene]

  ### compute null distribution of sum-of-weights
  pathways_lengths = sapply(pathways, length)
  all_pvalues = lapply(unique(pathways_lengths), function(plen){
    ### compute null distribution for given pathway length
    bg_idx = seq_len(length(bg_genes))
    null_scores = sapply(1:iter, function(it){
      iter_genes = sample(bg_idx, size = plen)
      weight_sum = sum(net[iter_genes, iter_genes], na.rm = na.rm)/2
      return(weight_sum)
    })

    null_dist =  fitdistr(null_scores, densfun = 'normal')
    null_mean = null_dist$estimate['mean']
    null_sd = null_dist$estimate['sd']

    ### compute pvalue for each gene set given length
    pathway_names = names(pathways_lengths[pathways_lengths==plen])
    pathway_pvalues = sapply(pathway_names, function(pname){
      p_genes = pathways[[pname]]
      observed_score = sum(net[p_genes, p_genes], na.rm = na.rm)/2
      empirical_p = (sum(null_scores >= observed_score)+1)/(length(null_scores)+1)
      normal_p = 1 - pnorm(observed_score, mean = null_mean, sd = null_sd)
      return(list(pathway = pname,
                  n.gene = plen,
                  null.mean = null_mean,
                  null.sd = null_sd,
                  p = normal_p,
                  p.empirical = empirical_p))

    })

    ret_df = as.data.frame(t(pathway_pvalues), stringsAsFactors = F)
    for(ci in seq_len(length(ret_df)))
      ret_df[,ci] = unlist(ret_df[,ci])  # convert list to atomic vector
    return(ret_df)
  })

  all_pvalues_df = as.data.frame(do.call(rbind, args = all_pvalues), stringsAsFactors = F)
  if(nrow(all_pvalues_df) > 0){
    all_pvalues_df$fdr = p.adjust(all_pvalues_df$p, method = 'BH')
    all_pvalues_df$fdr.empirical = p.adjust(all_pvalues_df$p.empirical, method = 'BH')
  }
  return(all_pvalues_df)
}
