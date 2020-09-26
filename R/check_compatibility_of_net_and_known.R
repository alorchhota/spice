check_row_col_compatibility_of_net_and_known <- function(net, known){
  if(!(is.matrix(net) || is.data.frame(net)))
    stop("net must be a matrix or a data.frame.")
  if(!(is.matrix(known) || is.data.frame(known)))
    stop('known must be a matrix or a data.frame.')
  if(any(dim(net) != dim(known)))
    stop('net and known must have same dimension.')

  if(xor(length(rownames(net))>0, length(colnames(net) >0)))
    stop("both rows and columns of net must be named or unnamed.")
  if(xor(length(rownames(known))>0, length(colnames(known) >0)))
    stop("both rows and columns of known must be named or unnamed.")
  if(xor(length(rownames(net))>0, length(rownames(known) >0)))
    stop("either both or none of net and known can have row names.")
  if(length(rownames(net)) >0 && length(rownames(known) >0)){
    if(nrow(net) != length(unique(rownames(net))) || nrow(known) != length(unique(rownames(known))))
      stop("duplicate rownames are not allowed in net and known.")
    if(!(setequal(rownames(net), colnames(net))) || !(setequal(rownames(known), colnames(known))))
      stop("rows and columns of net and known must have same set of genes.")
    if(!(setequal(rownames(net), rownames(known))))
      stop("net and known must have same set of genes")
  }
}

check_na.ignore_compatibility_of_net_and_known <- function(net, known, na.ignore){
  # assuming gene-orders of net and known are matched.
  if(!(na.ignore %in% c("known", "net", "any")))
    stop("values allowed for na.ignore: 'known', 'net', or 'any'.")
  if(na.ignore == 'net'){
    na_count = sum(is.na(known[!is.na(net)]))
    if(na_count > 0)
      stop("when na.ignore = 'net', known cannot have NA where net has non-NA.")
  } else if(na.ignore == 'known'){
    na_count = sum(is.na(net[!is.na(known)]))
    if(na_count > 0)
      stop("when na.ignore = 'known', net cannot have NA where known has non-NA.")
  }
}

check_matrix_symmetry <- function(x, varname){
  if(!is.matrix(x))
    stop(sprintf("%s must be a matrix.", varname))
  if (!isSymmetric(x))
    stop(sprintf("%s must be a symmetric matrix.", varname))
}

check_matrix_symmetry_of_net_and_known <- function(net, known){
  check_matrix_symmetry(net, varname = "net")
  check_matrix_symmetry(known, varname = "known")
}

check_value_range_of_net_and_known <- function(net, known, neg.treat){
  if(!(neg.treat %in% c("allow", "warn", "error")))
    stop("values allowed for neg.treat: 'allow', 'warn', or 'error'.")

  if(sum(known < 0, na.rm = T) || sum(known > 1, na.rm = T))
    stop('known must have values between 0 and 1')

  if(neg.treat == "warn" || neg.treat == "error"){
    raise_func = ifelse(neg.treat == "warn", warning, stop)
    if (sum(net < 0, na.rm = T)>0)
      raise_func("net has negative weight(s).")
  }
}
