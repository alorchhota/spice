#include <Rcpp.h>
using namespace Rcpp;

double maxC(double val1, double val2){
  if(NumericVector::is_na(val1))
    return(val2);
  if(NumericVector::is_na(val2))
    return(val1);
  return(std::max(val1, val2));
}

double minC(double val1, double val2){
  if(NumericVector::is_na(val1))
    return(val2);
  if(NumericVector::is_na(val2))
    return(val1);
  return(std::min(val1, val2));
}

void makeSymmetricMax(NumericMatrix x){
  size_t nrow = x.nrow();
  for (size_t i = 0; i < nrow-1; i++) {
    for (size_t j = i+1; j < nrow; j++) {
      x(i, j) = maxC(x(i, j), x(j, i));
      x(j, i) = x(i, j);
    }
  }
}

void makeSymmetricMin(NumericMatrix x){
  size_t nrow = x.nrow();
  for (size_t i = 0; i < nrow-1; i++) {
    for (size_t j = i+1; j < nrow; j++) {
      x(i, j) = minC(x(i, j), x(j, i));
      x(j, i) = x(i, j);
    }
  }
}

void makeSymmetricFromLowerTriangle(NumericMatrix x){
  size_t nrow = x.nrow();
  for (size_t i = 1; i < nrow; i++) {
    for (size_t j = 0; j < i; j++) {
      x(j, i) = x(i, j);
    }
  }
}

void makeSymmetricFromUpperTriangle(NumericMatrix x){
  size_t nrow = x.nrow();
  for (size_t i = 1; i < nrow; i++) {
    for (size_t j = 0; j < i; j++) {
      x(i, j) = x(j, i);
    }
  }
}

// [[Rcpp::export]]
void make_symmetric(NumericMatrix x, const String method="L"){
  // method == "L": retains lower triangular matrix
  // method == "U": retains upper triangular matrix
  // method == "max": max value
  // method == "min": min value
  if(method == "L"){
    makeSymmetricFromLowerTriangle(x);
  } else if(method == "U"){
    makeSymmetricFromUpperTriangle(x);
  } else if(method == "max"){
    makeSymmetricMax(x);
  } else if(method == "U"){
    makeSymmetricMin(x);
  }
}

// [[Rcpp::export]]
NumericVector get_lower_triangle_vector(const NumericMatrix x, const bool diagonal=false, const bool bycol=true){
  size_t nr = x.nrow(), nc = x.ncol();
  if(nr != nc)
    stop("x must be a square matrix.");

  size_t vector_length = diagonal ? nr * (nr + 1) / 2 : (nr - 1) * nr / 2;
  Rcpp::NumericVector lt(vector_length);
  size_t cur_len = 0;
  if(bycol == true){
    for(size_t j=0; j < nr; j++){
      size_t start_i = diagonal ? j : j+1;
      for(size_t i = start_i; i< nr; i++){
        lt.at(cur_len++) = x(i,j);
      }
    }
  } else {
    for(size_t i=0; i<nr; i++){
      size_t end_j = diagonal ? i+1 : i;
      for(size_t j = 0; j < end_j; j++){
        lt.at(cur_len++) = x(i,j);
      }
    }
  }
  return lt;
}


// [[Rcpp::export]]
void set_lower_triangle_vector(NumericMatrix x, const NumericVector lt, bool bycol = true){
  size_t nr = x.nrow(), nc = x.ncol();
  if(nr != nc)
    stop("x must be a square matrix.");

  bool diagonal = lt.size() == nr * (nr-1)/2 + nr;
  if(!diagonal && lt.size() != nr * (nr-1)/2){
    stop("sizes of x and v are not compatible.");
  }
  size_t cur_pos = 0;

  if(bycol == true){
    for(size_t j=0; j < nr; j++){
      size_t start_i = diagonal ? j : j+1;
      for(size_t i = start_i; i< nr; i++){
        x(i,j) = lt.at(cur_pos++);
      }
    }
  } else {
    for(size_t i=0; i<nr; i++){
      size_t end_j = diagonal ? i+1 : i;
      for(size_t j = 0; j < end_j; j++){
        x(i,j) = lt.at(cur_pos++);
      }
    }
  }
}


// [[Rcpp::export]]
void replace_NA_in_matrix(NumericMatrix x, const double value){
  for(NumericMatrix::iterator i = x.begin(); i != x.end(); ++i) {
    if(NumericMatrix::is_na(*i))
      *i = value;
  }
}

// [[Rcpp::export]]
void set_absolute_values_in_matrix(NumericMatrix x){
  for(NumericMatrix::iterator i = x.begin(); i != x.end(); ++i) {
    if(NumericMatrix::is_na(*i))
      *i = abs(*i);
  }
}

// [[Rcpp::export]]
void set_diag(NumericMatrix x, const double value){
  size_t nr = x.nrow(), nc = x.ncol();
  if(nr != nc)
    stop("x must be a square matrix.");
  for(size_t i=0; i<nr; i++){
    x(i,i) = value;
  }
}
