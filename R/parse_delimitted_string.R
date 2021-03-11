#' Function to parse delimitted string
#' 
#' This function splits the given string by a delimiter.
#' @param x character. String to be parsed.
#' @param delim character. The delimiter.
#' @param rm.empty logical. If TRUE, removes every empty string after parsing.
#' @return A vector of parsed character strings.
#' @export
#' @examples 
#' x = "one,two,three,"
#' parsed_strings = parse_delimitted_string(x, delim = ',', rm.empty = T)
parse_delimitted_string <- function(x, delim=',', rm.empty=T){
  require('stringr')
  stopifnot(class(x)=='character' && length(x)==1)
  parts = str_split(x, pattern = delim)[[1]]
  if(rm.empty==T){
    is_valid_parts = sapply(parts, function(s) nchar(s)>0)
    parts = parts[is_valid_parts]
  }
  return(parts)
}
