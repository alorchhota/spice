#include <Rcpp.h>
#include <cmath>

// following instructions to use bigmemory in Rcpp
// from https://gallery.rcpp.org/articles/using-bigmemory-with-rcpp/
//
// The next line is all it takes to find the bigmemory
// headers -- thanks to the magic of Rcpp attributes,
// and as bigmemory now accesses Boost headers from the BH package,
// we need to make sure we do so as well in this Rcpp::depends comment.
//
// Boost headers can generate some warning with the default compilation
// options for R.  To suppress these, we can enable C++11 mode which gets
// us 'long long' types.
//
// If your compiler is to old, just disable / remove the following line
// [[Rcpp::plugins(cpp11)]]

// Ensure R uses the headers for Boost (a dependency) and bigmemory
//
// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>

#include <numeric>

using namespace Rcpp;

// [[Rcpp::export]]
void update_rankprod_matrix(SEXP rankprod_bigmem, const NumericVector ranks){
  // tell Rcpp what class to use for big.matrix objects
  Rcpp::XPtr<BigMatrix> rankprod(rankprod_bigmem);
  // Create the matrix accessor so we can get at the elements of the matrix.
  MatrixAccessor<double> rankprod_acc(*rankprod);

  size_t nr = rankprod->nrow();
  if(ranks.size() != nr)
    stop("sizes of rankprod and ranks are not compatible.");
  for(size_t i=0; i<nr; i++){
    double cur_rank = ranks.at(i);
    if(NumericVector::is_na(cur_rank))
      continue;
    // bigmemory accesses column first, then row
    // (i,j)-th entry = rankprod_acc[j][i] -- notice order of j and i.
    rankprod_acc[0][i] = rankprod_acc[0][i] + log(cur_rank);
    rankprod_acc[1][i] = rankprod_acc[1][i] + 1;
  }
}

// [[Rcpp::export]]
NumericVector from_sampled_assoc_matrix_to_all_assoc_vector(
    const NumericMatrix sampled_assoc,
    const CharacterVector all_genes,
    const bool from_dist=false){

  // ######### Function in R ############################################################
  // all_assoc = matrix(NA, nrow = length(all_genes), ncol = length(all_genes),
  //                    dimnames = list(all_genes, all_genes))
  // all_assoc[rownames(sampled_assoc), colnames(sampled_assoc)] = sampled_assoc
  // all_assoc_vector = all_assoc[lower.tri(all_assoc)]
  // if from_dist==true:  all_assoc_vector = 1 - all_assoc_vector / max(all_assoc_vector, na.rm=T)
  // ####################################################################################

  long sampled_nr = sampled_assoc.nrow(), sampled_nc = sampled_assoc.ncol();
  if(sampled_nr != sampled_nc)
    stop("sampled_assoc must be a square matrix.");
  if(all_genes.size() < sampled_nr)
    stop("all_genes.size() must be >= sampled_assoc.nrow().");

  // create map from gene to sampled_assoc row (and col) index
  CharacterVector sampled_row_genes = rownames(sampled_assoc);
  CharacterVector sampled_col_genes = colnames(sampled_assoc);
  std::map<String,long> gene_to_sampled_assoc_row_index;
  std::map<String,long> gene_to_sampled_assoc_col_index;
  for(long i=0; i<all_genes.size(); i++){
    gene_to_sampled_assoc_row_index[all_genes(i)] = -1; // init to -1
    gene_to_sampled_assoc_col_index[all_genes(i)] = -1; // init to -1
  }
  for(long i=0; i<sampled_nr; i++){
    gene_to_sampled_assoc_row_index[sampled_row_genes(i)] = i;
    gene_to_sampled_assoc_col_index[sampled_col_genes(i)] = i;
  }

  // create map from all_assoc_row_idx (or col) to sampled_assoc_row_idx (or col)
  std::map<long,long> all_to_sampled_assoc_row_index;
  std::map<long,long> all_to_sampled_assoc_col_index;
  for(long i=0; i<all_genes.size(); i++){
    String g = all_genes(i);
    all_to_sampled_assoc_row_index[i] = gene_to_sampled_assoc_row_index[g];
    all_to_sampled_assoc_col_index[i] = gene_to_sampled_assoc_col_index[g];
  }

  // create lower triangle association vector
  long n_genes = all_genes.size();
  NumericVector assoc_vec(n_genes * (n_genes-1)/2, NA_REAL);    // init to NA
  long vidx = -1;
  for(long j=0; j<n_genes; j++){
    long s_j = all_to_sampled_assoc_col_index[j];
    if(s_j < 0){
      vidx += (n_genes - j -1);
      continue;
    }
    for(long i=j+1; i<n_genes; i++){
      long s_i = all_to_sampled_assoc_row_index[i];
      if(s_i < 0){
        vidx ++;
        continue;
      }
      assoc_vec(++vidx) = sampled_assoc(s_i,s_j);
    }
  }

  if(from_dist){
    Function rmax("max");
    double max_val = as<double>(rmax(assoc_vec, _["na.rm"]=true));
    for(NumericVector::iterator it = assoc_vec.begin(); it != assoc_vec.end(); it++){
      *it = 1.0-*it/max_val;
    }
  }

  return assoc_vec;
}

// [[Rcpp::export]]
NumericVector get_rankprod_vector_from_matrix(SEXP rankprod_bigmem){
  // tell Rcpp what class to use for big.matrix objects
  Rcpp::XPtr<BigMatrix> rankprod(rankprod_bigmem);
  // Create the matrix accessor so we can get at the elements of the matrix
  MatrixAccessor<double> rankprod_acc(*rankprod);

  size_t nr = rankprod->nrow();
  NumericVector rankprod_vector(nr, NA_REAL);
  for(size_t i=0; i<nr; i++){
    double count = rankprod_acc[1][i];
    if(count > 0){
      rankprod_vector(i) = exp(rankprod_acc[0][i] / count);
    }
  }
  return rankprod_vector;
}
