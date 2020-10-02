#' Row-wise min-max
#'
#' This function finds row-wise min-max values in two columns of a given dataframe.
#' @param df data frame with at least two columns.
#' @param col1 character. must be a column of df.
#' @param col2 character. must be a column of df.
#' @param unique if TRUE, duplicate rows are removed from the returned data frame.
#' @details
#' Given a data frame (df) and two of its columns (col1 and col2),
#' this function re-arranges the content of these columns so that
#' the first and second column contains the rowwise minimum and
#' the rowwise maximum of the two columns, respectively.
#' Factor values are sorted by their labels.
#' If there are other columns in df, those columns remain unchanged.
#' if unique = TRUE, the output data frame contains only unique rows.
#'
#' Note: both columns must have same type of data.
#' @return
#' A data frame with the same set of columns as df.
#' @examples
#' my_df = data.frame(A=c(1,5,2,9,1), B = LETTERS[1:5], C=c(1,3,5,7,9))
#' my_minmax_df = rminmax_df(my_df, col1='A', col2 = 'C')

rminmax_df <- function(df, col1 = colnames(df)[1], col2 = colnames(df)[2], unique = F){
  # return a dataframe with rowwise min-max values i.e.,
  # col1 will contian the minimum values between col1 and col2, and
  # col2 will contian the maximum values between col1 and col2
  # if unique = T, return a data frame with duplicated rows removed.


  stopifnot(class(df) == 'data.frame')
  stopifnot(ncol(df) >= 2)
  stopifnot(col1 %in% colnames(df))
  stopifnot(col2 %in% colnames(df))
  stopifnot(class(df[,col1]) == class(df[,col2]))

  is_factor = class(df[,col1]) == 'factor'
  if(is_factor){
    stopifnot(all(levels(df[,col1]) == levels(df[,col2])))
    factor_levels = levels(df[,col1])
    df[,col1] = as.character(df[,col1])
    df[,col2] = as.character(df[,col2])
  }

  values1 = pmin(df[,col1], df[,col2])
  values2 = pmax(df[,col1], df[,col2])

  if(is_factor){
    values1 = factor(values1, levels = factor_levels)
    values2 = factor(values2, levels = factor_levels)
  }

  df[,col1] = values1
  df[,col2] = values2

  if(unique)
    df = unique(df)

  return(df)
}
