#' SPICE: Spanning tree based inference of co-expression networks

#'
#' This function reconstructs a co-expression network from gene expression profiles.
#' @param expr matrix or data.frame. Processed gene expression data (gene x sample).
#' The rows must have unique gene names.
#' @param method a character string indicating how to compute association between a pair of genes.
#' This must be one of "pearson" (default), "spearman", "mi.empirical", "mi.spearman", and "mi.pearson". See details.
#' @param iter numeric. Number of iterations.
#' @param frac.gene numeric. Fraction of genes used in each iteration.
#' @param frac.sample numeric. Fraction of samples used in each iteration.
#' @param n.cores numeric. The number of cores to use.
#' @param rank.ties a character string indicating how ties are treated in ranking.
#' Accepted values include \code{"average", "first", "last", "random", "max", "min", "dense"}.
#' For details, see \code{ties.method} parameter of \code{\link[data.table]{frank}}.
#' @param filter.mat  NULL or a gene x gene matrix indicating which edges to exclude.
#' An edge is excluded if the corresponding entry is \code{TRUE} or greater than 0.
#' Row and column names of the matrix must be same as \code{rownames(expr)}.
#' @param weight.method a character string indicating how the weights are assigned.
#' It must be either \code{"qnorm"} (default) or \code{"inverse.rank"}.
#' @param adjust.weight logical. Should the weights be raised to a power to exhibit a scale-free topology.
#' @param adjust.weight.powers numeric vector. Powers for which the scale free topology fit indices are to be calculated.
#' @param adjust.weight.rsquared numeric. Desired minimum scale free topology fitting index R-squared.
#' If no power achieves \code{adjust.weight.rsquared}, the power with the maximum R-squared is chosen.
#' @param adjust.weight.bins numeric. Number of bins in connectivity histograms used to compute scale-free fit indices.
#' @param adjust.clr logical. Should the network weights be adjusted using context likelihood or relatedness?
#' See \code{\link[minet]{clr}} for details.
#' @param seed integer or NULL. Random number generator seed.
#' @param verbose logical or numeric. Print messages if verbose > 0.
#'
#' @importFrom foreach %dopar%
#' @importFrom Rcpp evalCpp
#' @details
#' The expression matrix \code{expr} is expected to be properly normalized, processed and free of batch effects.
#' This function may run slow when there are missing values in \code{expr}.
#'
#'
#'
#' @return Returns a matrix representing a co-expression network.
#' @export
#' @rawNamespace useDynLib(spice)
#' @examples
#' n_gene = 10
#' n_sample = 100
#' expr = matrix(rnorm(n_gene * n_sample),
#'               nrow = n_gene,
#'               ncol = n_sample,
#'               dimnames = list(sprintf("Gene%s", seq_len(n_gene)),
#'                               sprintf("Sample%s", seq_len(n_sample))))
#' spice_net = spice(expr, iter = 10, verbose = TRUE)

spice <- function(expr,
                  method='pearson',
                  iter=100,
                  frac.gene = 0.8,
                  frac.sample= 0.8,
                  n.cores=1,
                  rank.ties = "average",
                  filter.mat = NULL,
                  weight.method='qnorm',
                  adjust.weight=T,
                  adjust.weight.powers = c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
                  adjust.weight.rsquared = 0.8,
                  adjust.weight.bins = 10,
                  adjust.clr = F,
                  seed = NULL,
                  verbose = F){

  ### check parameters
  stopifnot(is.matrix(expr) || is.data.frame(expr))
  stopifnot(nrow(expr) > 1 && ncol(expr) > 1)
  stopifnot(length(unique(rownames(expr))) == nrow(expr))
  stopifnot(method %in% c('spearman', 'pearson', 'mi.empirical', 'mi.spearman', 'mi.pearson'))
  stopifnot(is.numeric(iter) && iter >= 1)
  stopifnot(is.numeric(frac.gene) && frac.gene > 0 && frac.gene <= 1)
  stopifnot(is.numeric(frac.sample) && frac.sample > 0 && frac.sample <= 1)
  stopifnot(is.numeric(n.cores) && n.cores >= 1)
  stopifnot(rank.ties %in% c("average", "first", "last", "random", "max", "min", "dense"))
  stopifnot(is.null(filter.mat) || (is.matrix(filter.mat) && setequal(rownames(filter.mat), rownames(expr)) && setequal(colnames(filter.mat), rownames(expr))))
  stopifnot(weight.method %in% c('inverse.rank', 'qnorm'))
  if(adjust.weight == T){
    stopifnot(is.numeric(adjust.weight.powers))
    stopifnot(is.numeric(adjust.weight.rsquared))
    stopifnot(is.numeric(adjust.weight.bins) && adjust.weight.bins>=2)
  }

  ### load required packages
  requireNamespace('parallel', quietly = T)
  requireNamespace('minet', quietly = T)
  requireNamespace('data.table', quietly = T)
  requireNamespace('foreach', quietly = T)
  requireNamespace('igraph', quietly = T)
  requireNamespace('SharedObject', quietly = T)
  requireNamespace('doParallel', quietly = T)
  requireNamespace('flock', quietly = T)

  ### initialize variables
  rankprod_matrix = matrix(data = 0, nrow = choose(nrow(expr), 2), ncol = 2)
  rankprod_matrix = SharedObject::share(rankprod_matrix, copyOnWrite=F)

  rankprod_lock = tempfile()
  expr_df = SharedObject::share(as.matrix(expr))
  rm("expr")
  tmp = gc(verbose = F)

  keep.mat = NULL
  if(!is.null(filter.mat)){
    # filter out edge if filter.mat contains > 0 or TRUE // keep if filter.mat <= 0 or FALSE
    keep.mat = matrix(T, nrow = nrow(expr_df), ncol =  nrow(expr_df), dimnames = list(rownames(expr_df), rownames(expr_df)))
    filter_row_genes = intersect(rownames(filter.mat), rownames(expr_df))
    filter_col_genes = intersect(colnames(filter.mat), rownames(expr_df))
    keep.mat[filter_row_genes, filter_col_genes] = keep.mat[filter_row_genes, filter_col_genes] * (filter.mat[filter_row_genes, filter_col_genes] <= 0)
    keep.mat[filter_col_genes, filter_row_genes] = keep.mat[filter_col_genes, filter_row_genes] * (t(filter.mat[filter_row_genes, filter_col_genes]) <= 0)
    diag(keep.mat) <- T
    keep.mat = SharedObject::share(keep.mat)
    tmp = gc(verbose = F)
  }

  get_expr_similarity <- function(expr_df, method='spearman', disc='equalfreq', nbins=sqrt(ncol(expr_df))){
    requireNamespace('infotheo', quietly = T)

    stopifnot(method %in% c('spearman', 'pearson', 'mi.empirical', 'mi.pearson', 'mi.spearman'))

    if(method %in% c('spearman', 'pearson')){
      if(sum(!is.finite(as.matrix(expr_df))) >0 ){
        pairwise.r = stats::cor(t(expr_df), method = method, use = "pairwise.complete.obs")
      } else {
        pairwise.r = stats::cor(t(expr_df), method = method)
      }
      replace_NA_in_matrix(pairwise.r, value = 0)
      mi = abs(pairwise.r)
      rm("pairwise.r")
      tmp = gc(verbose = F)
    } else if(method %in% c('mi.empirical')){
      requireNamespace('minet', quietly = T)
      expr_df = as.data.frame(t(expr_df)) # convert to samples x genes matrix, required for build.mim()
      if(disc != 'none')
        expr_df = infotheo::discretize(expr_df, disc = disc, nbins = nbins)
      tmp = gc(verbose = F)
      mi = minet::build.mim(dataset =  expr_df, estimator = method, disc = 'none')

      # normalize by sum of pairwise entropies
      h = apply(expr_df, 2, infotheo::entropy)
      pairwise_h_sum = utils::combn(h, 2, sum)
      pairwise_h_sum_mat = matrix(0, nrow=length(h), ncol=length(h), dimnames = list(names(h), names(h)))
      pairwise_h_sum_mat[lower.tri(pairwise_h_sum_mat)] = pairwise_h_sum
      diag(pairwise_h_sum_mat) = h
      rm(list = c("h", "pairwise_h_sum"))
      tmp = gc(verbose = F)
      pairwise_h_sum_mat = pairwise_h_sum_mat + t(pairwise_h_sum_mat)
      mi = 2 * mi / pairwise_h_sum_mat  # normalize MI
      rm("pairwise_h_sum_mat")
      tmp = gc(verbose = F)
    } else if(method %in% c('mi.pearson', 'mi.spearman')){
      requireNamespace('minet', quietly = T)
      expr_df = as.data.frame(t(expr_df)) # convert to samples x genes matrix, required for build.mim()
      mi.method = strsplit(method, split = "mi.")[[1]][2]
      mi = minet::build.mim(dataset =  expr_df, estimator = mi.method, disc = 'none')
      # normalize by diving by max
      diag(mi) = 0
      max_mi = max(mi)
      mi = mi / max_mi
      diag(mi) = 1
    }
    return(mi)
  }

  verbose_print <- function (msg, verbose = 1){
    if (verbose > 0) {
      print(sprintf("[%s] %s", format(Sys.time(), "%D %T"), msg))
    }
  }

  get_pairwise_assoc_per_iteraction <- function(it, frac.gene, frac.sample, use.keep.mat = F, verbose = F, seed = NULL){
    verbose_print(sprintf('spice iteration# %d', it), verbose = verbose)

    ### sample expression data
    if(!is.null(seed) && !is.na(seed) && is.finite(seed))
      set.seed(seed * 1000 + it)
    sampled_rows = sample(x = seq_len(nrow(expr_df)), size = round(nrow(expr_df) * frac.gene), replace = F)
    sampled_cols = sample(x = seq_len(ncol(expr_df)), size = round(ncol(expr_df) * frac.sample), replace = F)
    sampled_expr = expr_df[sampled_rows, sampled_cols, drop = F]

    ### compute pairwise similarity
    expr_sim = get_expr_similarity(expr_df = sampled_expr, method = method)
    rm("sampled_expr")
    tmp = gc(verbose = F)

    if(use.keep.mat){
      keep_genes = rownames(expr_sim)
      expr_sim[keep_genes, keep_genes] = keep.mat[keep_genes, keep_genes] * expr_sim[keep_genes, keep_genes]
    }

    ### compute maximum spanning
    kmst <- kruskal_mst_c(expr_sim, maximum = T, verbose = verbose)
    kmst_genes = rownames(expr_sim)
    kmst$from = kmst_genes[kmst$from]
    kmst$to = kmst_genes[kmst$to]
    sampled_mst <- igraph::graph_from_data_frame(kmst, directed=FALSE, vertices=kmst_genes)
    rm(list=c("kmst"))
    tmp = gc(verbose = F)

    ### update wgcna rank-prod
    wgcna.weights = from_sampled_assoc_matrix_to_all_assoc_vector(
      sampled_assoc = expr_sim,
      all_genes = rownames(expr_df),
      from_dist = F)
    wgcna.rank = data.table::frankv(wgcna.weights, order = -1L, na.last = "keep", ties.method = rank.ties)
    rm("wgcna.weights")
    tmp = gc(verbose = F)
    locked <- flock::lock(rankprod_lock)
    update_rankprod_matrix(rankprod_matrix, wgcna.rank)
    flock::unlock(locked)
    rm("wgcna.rank")
    tmp = gc(verbose = F)

    ### update spanning-tree rank-prod
    igraph::E(sampled_mst)$weight = 1 - abs(igraph::E(sampled_mst)$weight)  # convert correlation (or normalized mutual information) to distance
    sampled_distances = igraph::distances(sampled_mst)
    spanningtree.weights = from_sampled_assoc_matrix_to_all_assoc_vector(
      sampled_assoc = sampled_distances,
      all_genes = rownames(expr_df),
      from_dist = T)
    spanningtree.rank = data.table::frankv(spanningtree.weights, order = -1L, na.last = "keep", ties.method = rank.ties)
    rm("spanningtree.weights")
    tmp = gc(verbose = F)
    locked <- flock::lock(rankprod_lock)
    update_rankprod_matrix(rankprod_matrix, spanningtree.rank)
    flock::unlock(locked)
    rm("spanningtree.rank")
    tmp = gc(verbose = F)

    return(NA)
  }

  ### run iterations
  if(n.cores == 1){
    for(it in seq_len(iter)){
      get_pairwise_assoc_per_iteraction(it, frac.gene=frac.gene, frac.sample=frac.sample, use.keep.mat = !is.null(keep.mat), verbose = verbose, seed = seed)
    }
  } else {
    if(verbose>0){
      cl <- parallel::makeCluster(n.cores, outfile = "")
    } else {
      cl <- parallel::makeCluster(n.cores)
    }
    on.exit(parallel::stopCluster(cl))
    tmp <- parallel::clusterEvalQ(cl, {
      requireNamespace('flock', quietly = T)
      requireNamespace('data.table', quietly = T)
      requireNamespace('igraph', quietly = T)
    })
    shared_varlist = list( "get_pairwise_assoc_per_iteraction",
                           "expr_df",
                           "rankprod_matrix",
                           "get_expr_similarity")
    if(!is.null(keep.mat))
      shared_varlist = append(shared_varlist, "keep.mat")
    parallel::clusterExport(cl, varlist = shared_varlist, envir = environment())

    doParallel::registerDoParallel(cl)
    tmp <- foreach::foreach(it = seq_len(iter)) %dopar% {
      get_pairwise_assoc_per_iteraction(it, frac.gene, frac.sample, use.keep.mat = !is.null(keep.mat), verbose = verbose, seed = seed)
    }
    tmp = gc(verbose = F)
  }

  ### compute rank.prod from aggregated iterations
  #rank_products = exp(rankprod_matrix[,1] / rankprod_matrix[,2])
  rank_products = get_rankprod_vector_from_matrix(rankprod_matrix)
  rm("rankprod_matrix")
  tmp = gc(verbose = F)

  ### build network matrix based on rank-products
  net = matrix(0, nrow = nrow(expr_df), ncol = nrow(expr_df), dimnames = list(rownames(expr_df), rownames(expr_df)))
  max_possible_rank = choose(round(nrow(expr_df) * frac.gene), 2)
  if(weight.method == 'inverse.rank'){
    set_lower_triangle_vector(net, 1.0/rank_products)   # small rankprod => high weight
  } else if(weight.method == 'qnorm') {
    wnorm = abs(stats::qnorm(p=rank_products/(2*max_possible_rank))) # only left half of normal distribution
    max_possible_wnorm = abs(stats::qnorm(p=1/(2*max_possible_rank)))
    wnorm = wnorm / max_possible_wnorm
    set_lower_triangle_vector(net, wnorm)
    rm("wnorm")
    tmp = gc(verbose = F)
  }
  make_symmetric(net, method = "L")      # symmetric
  set_diag(net, 1)
  tmp = gc(verbose = F)

  ### adjust weights by powering them, similar to WGCNA scale-free thresholding
  if(adjust.weight == T){
    scale_free_fit <- function(k, nBreaks = 10){
      # this function was adapted from WGCNA::scaleFreeFitIndex()
      discretized.k = cut(k, nBreaks)       # which bin each degree(connectivity) falls in
      dk = tapply(k, discretized.k, mean)   # average degree in each bin
      p.dk = as.vector(tapply(k, discretized.k, length)/length(k))     # probability of each bin
      # for empty bins, use degree=midpoint and prob=0
      breaks1 = seq(from = min(k), to = max(k), length = nBreaks + 1)
      midpoints = (breaks1[1:nBreaks] + breaks1[2:(nBreaks+1)])/2
      dk = ifelse(is.na(dk), midpoints, dk)
      dk = ifelse(dk == 0, midpoints, dk)
      p.dk = ifelse(is.na(p.dk), 0, p.dk)
      # fit log(p.dk) ~ log(dk)
      log.dk = as.vector(log10(dk))
      log.p.dk = as.numeric(log10(p.dk + 1e-9))
      lm1 = stats::lm(log.p.dk ~ log.dk)
      return(list(rsquared = summary(lm1)$r.squared,
                  slope = summary(lm1)$coefficients[2, 1]))
    }
    scale_free_stats = sapply(adjust.weight.powers, function(cur_power){
      power_net = net ^ cur_power
      connectivities = rowSums(power_net, na.rm = T) - 1
      rm("power_net")
      tmp = gc(verbose = F)
      sff = scale_free_fit(k = connectivities, nBreaks = adjust.weight.bins)
      return(list(rsquared=sff$rsquared, slope=sff$slope))
    })
    rsquares = as.numeric(scale_free_stats[1,])
    rsquare_cut = min(max(rsquares), adjust.weight.rsquared)
    estimated_power_idx = which(rsquares>=rsquare_cut, arr.ind = T)[1]
    estimated_power = adjust.weight.powers[estimated_power_idx]
    estimated_rsquare = rsquares[estimated_power_idx]
    verbose_print(sprintf("Estimated power: %f (R2=%.2f)", estimated_power, estimated_rsquare), verbose = verbose)
    net = net ^ estimated_power
    tmp = gc(verbose = F)
  }

  ### make network sparse using minet::clr
  if(adjust.clr == T){
    net = minet::clr(net, skipDiagonal = T)
    tmp = gc(verbose = F)
  }

  ### remove data before return
  unlink(rankprod_lock)
  rm("expr_df")
  if(!is.null(keep.mat))
    rm(list = "keep.mat")
  tmp <- gc(verbose = F)

  return(net)
}
