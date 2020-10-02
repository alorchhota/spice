#' Protein-protein interaction (PPIs) from STRING.
#'
#' This function generates a matrix of PPI scores available in STRING database
#' within a given set of genes.
#'
#' @param genes character vector. gene names (HGNC gene symbols).
#'
#' @param directed logical. Directed interactions?
#'
#' @param version character. Version of STRING database.
#'
#' @param species numeric. Species code for STRING (9606 for homo sapiens).
#'
#' @param score.threshold numeric. A value between 0 and 1000.
#' Any score below \code{score.threshold} will be replaced by 0.
#'
#' @param string.dir character. Local directory where STRING database
#' will be stored locally so that it can be used offline.
#' If \code{string.dir = ""}, a temporary directory is used.
#'
#' @param max.homology.bitscore numeric. Maximum homology score.
#' Use \code{0} to filter all interactions between homologous proteins, and
#' \code{Inf} not to filter any interaction.
#'
#' @param benchmark.pathway character. Pathway type to benchmark.
#' Examples: \code{"KEGG", "REACTOME", "BIOCARTA", "Disease",
#' "Pfam", "NCI", "ECOCYC"}.
#' Use \code{NULL} to avoid benchmarking.
#' If a pathway type is given, an interaction will have a non-zero score
#' only if both genes are present in at least one pathway.
#'
#' @param normalize.score logical. Should the scores be normalized?
#' If \code{TRUE}, the scores will be divided by 1000.
#' Consequently, each score will be between 0 and 1.

#' @details
#' Note: all the genes are expected as HGNC gene symbols.
#'
#' @return A 2-dimensional matrix with interactions between the given genes.
#' Note: any gene that cannot be mapped to a protein in STRING will be excluded.
#'
#' @export
#' @examples
#' genes = c("TP53", "RBM3", "SF3", "LIM12", "MDM4", "TMEM160", "TP53BP2", "MDM2", "PDR", "MEG3", "EGFR")
#' interactions = get_string_ppi_matrix(genes, version = "11")

get_string_ppi_matrix <- function(genes,
                                  directed = FALSE,
                                  version="11",
                                  species=9606,
                                  score.threshold=0,
                                  string.dir="",
                                  max.homology.bitscore=Inf,
                                  benchmark.pathway=NULL,
                                  normalize.score = T){

  require('STRINGdb')
  require('reshape2')

  ### create a stringdb object
  string_db <- STRINGdb$new( version=version,
                             species=species,
                             score_threshold=score.threshold,
                             input_directory=string.dir)

  ### map gene names to string ids
  bg_df = data.frame(gene = genes)
  bg_string_df = string_db$map( bg_df, 'gene', removeUnmappedRows = TRUE, quiet = T)
  bg_string_ids = unique(bg_string_df$STRING_id)
  bg_genes = unique(bg_string_df$gene)
  if(length(bg_genes) < 2){
    warning("less than 2 genes mapped to STRING.")
    ppi_score_mat = matrix(0,
                           nrow = length(bg_genes),
                           ncol = length(bg_genes),
                           dimnames = list(bg_genes, bg_genes))

    return(ppi_score_mat)
  }

  ### get background interactions
  all_interactions_df = string_db$get_interactions(bg_string_ids)
  all_interactions_df = unique(all_interactions_df[,c('from','to', 'combined_score'), drop=F])
  colnames(all_interactions_df) = c('proteinA','proteinB','score')

  ### filter interactions based on provided settings
  if(is.finite(max.homology.bitscore)){
    all_interactions_df = string_db$remove_homologous_interactions(interactions_dataframe = all_interactions_df,
                                                                   bitscore_threshold = max.homology.bitscore)
  }
  if(!is.null(benchmark.pathway)){
    all_interactions_df = string_db$benchmark_ppi(interactions_dataframe = all_interactions_df,
                                                  pathwayType = benchmark.pathway,
                                                  max_homology_bitscore = max.homology.bitscore)
    all_interactions_df = all_interactions_df[all_interactions_df$outcome == "TP",
                                              c('proteinA', 'proteinB', 'score'),
                                              drop=F]
  }
  if(!directed)
    all_interactions_df = rminmax_df(all_interactions_df,
                                     col1 = 'proteinA',
                                     col2 = 'proteinB',
                                     unique = T)

  ### add gene information in the interaction data.frame
  all_gene_interactions_df = merge(all_interactions_df, bg_string_df,
                                   by.x = 'proteinA',
                                   by.y = 'STRING_id')
  all_gene_interactions_df = merge(all_gene_interactions_df, bg_string_df,
                                   by.x = 'proteinB',
                                   by.y = 'STRING_id',
                                   suffixes = c('A','B'))
  all_gene_interactions_df = all_gene_interactions_df[,c('geneA','geneB',
                                                         'proteinA','proteinB','score'),
                                                      drop = F]

  ### convert gene interaction to a matrix
  known_ppi_score_mat = suppressWarnings(acast(data = all_gene_interactions_df,
                                               formula = geneA ~ geneB,
                                               value.var = "score",
                                               fun.aggregate = max,
                                               na.rm = T))
  known_ppi_score_mat[is.infinite(known_ppi_score_mat)] = NA
  rm(all_gene_interactions_df)

  ### create a matrix with all genes mapped in string
  ppi_score_mat = matrix(NA,
                         nrow = length(bg_genes),
                         ncol = length(bg_genes),
                         dimnames = list(bg_genes, bg_genes))
  ppi_score_mat[rownames(known_ppi_score_mat), colnames(known_ppi_score_mat)] = known_ppi_score_mat
  rm(known_ppi_score_mat)
  if(!directed)
    ppi_score_mat = pmax(ppi_score_mat, t(ppi_score_mat), na.rm = T)
  ppi_score_mat[is.na(ppi_score_mat)] = 0

  ### normalize score
  if(normalize.score)
    ppi_score_mat = ppi_score_mat/1000

  return(ppi_score_mat)
}
