#' SPICE: Spanning tree based inference of co-expression networks
#'
#' This function reconstructs a co-expression network from gene expression profiles.
#' @param expr data.frame or matrix. Processed gene expression data (gene x sample). See details.
#' @param method a character string indicating how to compute association between a pair of genes.
#' This must be one of "pearson" (default), "spearman", "mi.empirical", "mi.spearman", and "mi.pearson". See details.
#' @param iter numeric. Number of iterations.
#' @param frac.row numeric. Fraction of genes sampled without replacement in each iteration.
#' @param frac.col numeric. Fraction of samples sampled without replacement in each iteration.
#' @param n.cores numeric. The number of cores to use.
#' @param rank.ties a character string indicating how ties are treated in ranking. Accepted values include \code{"average", "first", "last", "random", "max", "min", "dense"}. For details, see \code{ties.method} parameter of the  \code{frank} function in  \code{data.table} package.
#' @param filter.mat  matrix or NULL. Which edges to exclude.
#' @param weight.method a character string indicating how the weights are assigned. It must be either \code{"qnorm"} (default) or \code{"inverse.rank"}.
#' @param adjust.weight logical. Should the weights be raised to a power to exhibit a scale-free topology.
#' @param powerVector numeric vector. Powers for which the scale free topology fit indices are to be calculated.
#' @param RsquaredCut numeric. Desired minimum scale free topology fitting index R^2.
#' @param removeFirst logical. Should the first bin be removed from the connectivity histogram while computing scale-free topology fitting indices?
#' @param nBreaks numeric. Number of bins in connectivity histograms.
#' @param adjust.clr logical. Should the network weights be adjusted using context likelihood or relatedness? See \code{minet::clr} for details.
#' @param seed integer or NULL. Random number generator seed.
#' @param verbose logical or numeric. Print messages if verbose > 0.
#'
#' @details
#' The rows of \code{expr} must have names.
#'
#'
#'
#' @return Returns a matrix representing a co-expression network.
#' @export
#' @examples
#' n_gene = 10
#' n_sample = 100
#' expr = matrix(rnorm(n_gene * n_sample),
#'               nrow = n_gene,
#'               ncol = n_sample,
#'               dimnames = list(sprintf("Gene%s", seq_len(n_gene)),
#'                               sprintf("Sample%s", seq_len(n_sample))))
#' spice_net = spice(expr, iter = 10, verbose = T)

spice <- function(expr,
                  method='pearson',
                  iter=100,
                  frac.row = 0.8,
                  frac.col= 0.8,
                  n.cores=1,
                  rank.ties = "average",
                  filter.mat = NULL,
                  weight.method='qnorm',
                  adjust.weight=T,
                  powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
                  RsquaredCut = 0.8,
                  removeFirst = FALSE,
                  nBreaks = 10,
                  adjust.clr = F,
                  seed = NULL,
                  verbose = F){
  suppressMessages(require('miscutil'))
  suppressMessages(require('parallel'))
  suppressMessages(require('WGCNA'))
  suppressMessages(require('minet'))
  suppressMessages(require('data.table'))
  suppressMessages(require('foreach'))
  suppressMessages(require('igraph'))
  suppressMessages(require('SharedObject'))
  suppressMessages(require('doParallel'))
  suppressMessages(require('flock'))

  stopifnot(method %in% c('spearman', 'pearson', 'mi.empirical', 'mi.spearman', 'mi.pearson'))
  stopifnot(weight.method %in% c('inverse.rank', 'qnorm'))

  rankprod_matrix = matrix(data = 0, nrow = choose(nrow(expr), 2), ncol = 2)
  rankprod_matrix = share(rankprod_matrix, copyOnWrite=F)
  tmp = gc(verbose = F)

  rankprod_lock = tempfile()
  expr_df = share(as.matrix(expr))
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
    keep.mat = share(keep.mat)
    tmp = gc(verbose = F)
  }

  get_expr_similarity <- function(expr_df, method='spearman', disc='equalfreq', nbins=sqrt(ncol(expr_df))){
    suppressMessages(require(infotheo))

    stopifnot(method %in% c('spearman', 'pearson', 'mi.empirical', 'mi.pearson', 'mi.spearman'))

    if(method %in% c('spearman', 'pearson')){
      if(sum(!is.finite(as.matrix(expr_df))) >0 ){
        pairwise.r = cor(t(expr_df), method = method, use = "pairwise.complete.obs")
      } else {
        pairwise.r = cor(t(expr_df), method = method)
      }
      replace_NA_in_matrix(pairwise.r, value = 0)
      mi = abs(pairwise.r)
      rm("pairwise.r")
      tmp = gc(verbose = F)
    } else if(method %in% c('mi.empirical')){
      suppressMessages(require(minet))
      expr_df = as.data.frame(t(expr_df)) # convert to samples x genes matrix, required for build.mim()
      if(disc != 'none')
        expr_df = infotheo::discretize(expr_df, disc = disc, nbins = nbins)
      tmp = gc(verbose = F)
      mi = minet::build.mim(dataset =  expr_df, estimator = method, disc = 'none')

      # normalize by sum of pairwise entropies
      h = apply(expr_df, 2, entropy)
      pairwise_h_sum = combn(h, 2, sum)
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
      suppressMessages(require(minet))
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

  get_pairwise_assoc_per_iteraction <- function(it, frac.row, frac.col, use.keep.mat = F, verbose = F, seed = NULL){
    verbose_print(sprintf('spice iteration# %d', it), verbose = verbose)
    if(!is.null(seed))
      set.seed(seed * 1000 + it)
    sampled_expr = miscutil::sample_df(x = expr_df,
                                       size.row = frac.row * 100,
                                       size.col = frac.col * 100,
                                       unit = 'percent',
                                       replace.row = F,
                                       replace.col = F)
    sampled_expr = sampled_expr[sample(rownames(sampled_expr), size = nrow(sampled_expr), replace = F),,drop=F]

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
    sampled_mst <- graph_from_data_frame(kmst, directed=FALSE, vertices=kmst_genes)
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

    E(sampled_mst)$weight = 1 - abs(E(sampled_mst)$weight)  # convert correlation (or normalized mutual information) to distance
    sampled_distances = distances(sampled_mst)
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

  if(n.cores == 0){
    for(it in seq_len(iter)){
      get_pairwise_assoc_per_iteraction(it, frac.row=frac.row, frac.col=frac.col, use.keep.mat = !is.null(keep.mat), verbose = verbose, seed = seed)
    }
  } else {
    if(verbose>0){
      cl <- parallel::makeCluster(n.cores, outfile = "")
    } else {
      cl <- parallel::makeCluster(n.cores)
    }
    on.exit(parallel::stopCluster(cl))
    tmp <- clusterEvalQ(cl, {
      suppressMessages(require("miscutil"))
      suppressMessages(require("flock"))
      suppressMessages(require("data.table"))
      suppressMessages(require("igraph"))
    })
    shared_varlist = list( "get_pairwise_assoc_per_iteraction",
                           "expr_df",
                           "rankprod_matrix",
                           "get_expr_similarity")
    if(!is.null(keep.mat))
      shared_varlist = append(shared_varlist, "keep.mat")
    clusterExport(cl, varlist = shared_varlist, envir = environment())

    doParallel::registerDoParallel(cl)
    tmp <- foreach(it = seq_len(iter)) %dopar% {
      get_pairwise_assoc_per_iteraction(it, frac.row, frac.col, use.keep.mat = !is.null(keep.mat), verbose = verbose, seed = seed)
    }
    tmp = gc(verbose = F)
  }

  #rank_products = exp(rankprod_matrix[,1] / rankprod_matrix[,2])
  rank_products = get_rankprod_vector_from_matrix(rankprod_matrix)
  rm("rankprod_matrix")
  tmp = gc(verbose = F)

  net = matrix(0, nrow = nrow(expr_df), ncol = nrow(expr_df), dimnames = list(rownames(expr_df), rownames(expr_df)))
  max_possible_rank = choose(round(nrow(expr_df) * frac.row), 2)
  if(weight.method == 'inverse.rank'){
    set_lower_triangle_vector(net, 1.0/rank_products)   # small rankprod => high weight
  } else if(weight.method == 'qnorm') {
    wnorm = abs(qnorm(p=rank_products/(2*max_possible_rank))) # only left half of normal distribution
    max_possible_wnorm = abs(qnorm(p=1/(2*max_possible_rank)))
    wnorm = wnorm / max_possible_wnorm
    set_lower_triangle_vector(net, wnorm)
    rm("wnorm")
    tmp = gc(verbose = F)
  }
  make_symmetric(net, method = "L")      # symmetric
  set_diag(net, 1)
  tmp = gc(verbose = F)

  if(adjust.weight == T){
    ### adjust weights by powering them, similar to WGCNA scale-free thresholding
    scale_free_stats = sapply(powerVector, function(cur_power){
      power_net = net ^ cur_power
      connectivities = rowSums(power_net, na.rm = T) - 1
      rm("power_net")
      tmp = gc(verbose = F)
      SFT1 = scaleFreeFitIndex(k = connectivities, nBreaks = nBreaks, removeFirst = removeFirst)
      return(list(Rsquared.SFT=SFT1$Rsquared.SFT, slope.SFT=SFT1$slope.SFT))
    })
    rsquares = as.numeric(scale_free_stats[1,])
    rsquare_cut = min(max(rsquares), RsquaredCut)
    estimated_power_idx = which(rsquares>=rsquare_cut, arr.ind = T)[1]
    estimated_power = powerVector[estimated_power_idx]
    estimated_rsquare = rsquares[estimated_power_idx]
    verbose_print(sprintf("Estimated power: %f (R2=%.2f)", estimated_power, estimated_rsquare), verbose = verbose)
    net = net ^ estimated_power
    tmp = gc(verbose = F)
  }

  if(adjust.clr == T){
    net = clr(net, skipDiagonal = T)
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
