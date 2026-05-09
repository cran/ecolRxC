#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix IPF2_cpp(NumericMatrix matriz,
                       NumericVector vector_columna,
                       NumericVector vector_fila,
                       double precision = 0.01) {
  
  int nf = matriz.nrow();
  int nc = matriz.ncol();
  
  NumericMatrix X2 = clone(matriz);
  NumericMatrix X1(nf, nc);
  
  // ---- Primera iteración fuera del bucle (igual que en R) ------------------
  // Row scaling
  for (int i = 0; i < nf; i++) {
    double row_sum = 0.0;
    for (int j = 0; j < nc; j++) row_sum += X2(i, j);
    double scale = (row_sum == 0.0 || ISNAN(row_sum)) ? 0.0 : vector_fila[i] / row_sum;
    for (int j = 0; j < nc; j++) X1(i, j) = X2(i, j) * scale;
  }
  
  // Column scaling
  for (int j = 0; j < nc; j++) {
    double col_sum = 0.0;
    for (int i = 0; i < nf; i++) col_sum += X1(i, j);
    double scale = (col_sum == 0.0 || ISNAN(col_sum)) ? 0.0 : vector_columna[j] / col_sum;
    for (int i = 0; i < nf; i++) X2(i, j) = X1(i, j) * scale;
  }
  
  // Criterio de convergencia: max(dif.r, dif.c)  ← igual que en R
  double dif_r = 0.0, dif_c = 0.0;
  for (int i = 0; i < nf; i++) {
    double rs = 0.0;
    for (int j = 0; j < nc; j++) rs += X2(i, j);
    dif_r += std::abs(vector_fila[i] - rs);
  }
  for (int j = 0; j < nc; j++) {
    double cs = 0.0;
    for (int i = 0; i < nf; i++) cs += X2(i, j);
    dif_c += std::abs(vector_columna[j] - cs);
  }
  double dif = std::max(dif_r, dif_c);
  
  int iter = 0;
  
  // ---- Bucle iterativo ------------------------------------------------------
  while (dif > precision && iter < 1000) {
    
    // Row scaling
    for (int i = 0; i < nf; i++) {
      double row_sum = 0.0;
      for (int j = 0; j < nc; j++) row_sum += X2(i, j);
      double scale = (row_sum == 0.0 || ISNAN(row_sum)) ? 0.0 : vector_fila[i] / row_sum;
      for (int j = 0; j < nc; j++) X1(i, j) = X2(i, j) * scale;
    }
    
    // Column scaling
    for (int j = 0; j < nc; j++) {
      double col_sum = 0.0;
      for (int i = 0; i < nf; i++) col_sum += X1(i, j);
      double scale = (col_sum == 0.0 || ISNAN(col_sum)) ? 0.0 : vector_columna[j] / col_sum;
      for (int i = 0; i < nf; i++) X2(i, j) = X1(i, j) * scale;
    }
    
    // Convergencia: max(dif.r, dif.c)
    dif_r = 0.0; dif_c = 0.0;
    for (int i = 0; i < nf; i++) {
      double rs = 0.0;
      for (int j = 0; j < nc; j++) rs += X2(i, j);
      dif_r += std::abs(vector_fila[i] - rs);
    }
    for (int j = 0; j < nc; j++) {
      double cs = 0.0;
      for (int i = 0; i < nf; i++) cs += X2(i, j);
      dif_c += std::abs(vector_columna[j] - cs);
    }
    dif = std::max(dif_r, dif_c);
    
    iter++;
  }
  
  return X2;
}
