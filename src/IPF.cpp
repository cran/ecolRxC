#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix IPF_cpp(NumericMatrix matriz,
                      NumericVector vector_columna,
                      NumericVector vector_fila,
                      double precision = 1e-6) {
  
  int nf = matriz.nrow();
  int nc = matriz.ncol();
  
  NumericMatrix X2 = clone(matriz);
  NumericMatrix X1(nf, nc);
  
  // ---- Primera iteración fuera del bucle (igual que en R) ------------------
  // Row scaling: R1 = diag(vector.fila) %*% diag(1/rowSums(matriz))
  for (int i = 0; i < nf; i++) {
    double row_sum = 0.0;
    for (int j = 0; j < nc; j++) row_sum += X2(i, j);
    double scale = (row_sum == 0.0 || ISNAN(row_sum)) ? 0.0 : vector_fila[i] / row_sum;
    for (int j = 0; j < nc; j++) X1(i, j) = X2(i, j) * scale;
  }
  
  // Column scaling: S1 = diag(vector.columna) * diag(1/colSums(X1))
  for (int j = 0; j < nc; j++) {
    double col_sum = 0.0;
    for (int i = 0; i < nf; i++) col_sum += X1(i, j);
    double scale = (col_sum == 0.0 || ISNAN(col_sum)) ? 0.0 : vector_columna[j] / col_sum;
    for (int i = 0; i < nf; i++) X2(i, j) = X1(i, j) * scale;
  }
  
  // Criterio de convergencia: max(abs(X2 - matriz))  ← igual que en R
  double dif = 0.0;
  for (int i = 0; i < nf; i++)
    for (int j = 0; j < nc; j++)
      dif = std::max(dif, std::abs(X2(i, j) - matriz(i, j)));
  
  // ---- Bucle iterativo ------------------------------------------------------
  while (dif > precision) {
    
    // Guardar X2 anterior para calcular dif al final
    NumericMatrix X_prev = clone(X2);
    
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
    
    // Convergencia: max(abs(X2 - X_prev))
    dif = 0.0;
    for (int i = 0; i < nf; i++)
      for (int j = 0; j < nc; j++)
        dif = std::max(dif, std::abs(X2(i, j) - X_prev(i, j)));
  }
  
  // ---- Ajuste final (igual que en R) ----------------------------------------
  // X2 <- X2 * (vector.fila / rowSums(X2))
  for (int i = 0; i < nf; i++) {
    double row_sum = 0.0;
    for (int j = 0; j < nc; j++) row_sum += X2(i, j);
    double scale = (row_sum == 0.0 || ISNAN(row_sum)) ? 0.0 : vector_fila[i] / row_sum;
    for (int j = 0; j < nc; j++) {
      X2(i, j) *= scale;
      if (ISNAN(X2(i, j)) || std::isinf(X2(i, j))) X2(i, j) = 0.0;
    }
  }
  
  // X2 <- t(t(X2) * (vector.columna / colSums(X2)))
  for (int j = 0; j < nc; j++) {
    double col_sum = 0.0;
    for (int i = 0; i < nf; i++) col_sum += X2(i, j);
    double scale = (col_sum == 0.0 || ISNAN(col_sum)) ? 0.0 : vector_columna[j] / col_sum;
    for (int i = 0; i < nf; i++) {
      X2(i, j) *= scale;
      if (ISNAN(X2(i, j)) || std::isinf(X2(i, j))) X2(i, j) = 0.0;
    }
  }
  
  return X2;
}
