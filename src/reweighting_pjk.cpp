// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------------------
// Declaraciones forward de IPF_cpp e IPF2_cpp
// ---------------------------------------------------------------------------
NumericMatrix IPF_cpp(NumericMatrix matriz,
                      NumericVector vector_columna,
                      NumericVector vector_fila,
                      double precision = 1e-6);

NumericMatrix IPF2_cpp(NumericMatrix matriz,
                       NumericVector vector_columna,
                       NumericVector vector_fila,
                       double precision = 0.01);

// ---------------------------------------------------------------------------
// Helpers internos
// ---------------------------------------------------------------------------
static NumericMatrix arma_to_rcpp(const arma::mat& m) {
  int nr = m.n_rows, nc = m.n_cols;
  NumericMatrix out(nr, nc);
  for (int j = 0; j < nc; ++j)
    for (int i = 0; i < nr; ++i)
      out(i, j) = m(i, j);
  return out;
}

static arma::mat rcpp_to_arma(const NumericMatrix& m) {
  int nr = m.nrow(), nc = m.ncol();
  arma::mat out(nr, nc);
  for (int j = 0; j < nc; ++j)
    for (int i = 0; i < nr; ++i)
      out(i, j) = m(i, j);
  return out;
}

// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List reweighting_pjk_cpp(
    NumericVector pjk_flat,  // pjk.array aplanado con as.numeric()
    IntegerVector dims,      // dim(pjk.array) → c(J, K, I)
    int ref1,                // 1-based
    int ref2,                // 1-based
    std::string scale,       // "logit" o "probit"
    NumericMatrix x,         // [I x J]
    NumericMatrix y,         // [I x K]
    double tol = 1e-6
) {
  Function pnorm2d_fn = Environment::namespace_env("ecolRxC")["pnorm2d_mvtnorm"];
  Function qnorm_fn("qnorm");
  
  // ---- Dimensiones ---------------------------------------------------------
  int J = dims[0], K = dims[1], I = dims[2];
  int ref1_0 = ref1 - 1;
  int ref2_0 = ref2 - 1;
  
  arma::cube pjk_array(pjk_flat.begin(), J, K, I);
  arma::mat  xm(x.begin(), I, J, false);
  arma::mat  ym(y.begin(), I, K, false);
  
  std::vector<int> parties1, parties2;
  for (int j = 0; j < J; ++j) if (j != ref1_0) parties1.push_back(j);
  for (int k = 0; k < K; ++k) if (k != ref2_0) parties2.push_back(k);
  
  // ---- Marginales observados -----------------------------------------------
  arma::mat p_obs_rows(I, J), p_obs_cols(I, K);
  for (int i = 0; i < I; ++i) {
    double rs_x = arma::sum(xm.row(i));
    double rs_y = arma::sum(ym.row(i));
    for (int j = 0; j < J; ++j) p_obs_rows(i, j) = xm(i, j) / rs_x;
    for (int k = 0; k < K; ++k) p_obs_cols(i, k) = ym(i, k) / rs_y;
  }
  
  // ---- Bucle iterativo ------------------------------------------------------
  double suma0 = 0.0;
  double suma1 = arma::accu(pjk_array);
  double dif   = std::abs(suma0 - suma1);
  int    cont  = 0;
  
  while (dif > tol && cont < 10) {
    
    for (int jj : parties1) {
      for (int kk : parties2) {
        
        arma::mat suma_rows(2, I, arma::fill::zeros);
        arma::mat suma_cols(2, I, arma::fill::zeros);
        
        for (int i = 0; i < I; ++i) {
          double a = pjk_array(jj,     kk,     i);
          double b = pjk_array(jj,     ref2_0, i);
          double c = pjk_array(ref1_0, kk,     i);
          double d = pjk_array(ref1_0, ref2_0, i);
          suma_rows(0, i) = a + b;
          suma_rows(1, i) = c + d;
          suma_cols(0, i) = a + c;
          suma_cols(1, i) = b + d;
        }
        
        arma::vec ps1(I), ps2(I), weight(I);
        for (int i = 0; i < I; ++i) {
          double slice_sum = suma_rows(0,i) + suma_rows(1,i);
          if (scale == "logit") {
            ps1(i) = suma_rows(0,i) / suma_rows(1,i);
            ps2(i) = suma_cols(0,i) / suma_cols(1,i);
          } else {
            ps1(i) = suma_rows(0,i) / slice_sum;
            ps2(i) = suma_cols(0,i) / (suma_cols(0,i) + suma_cols(1,i));
          }
          weight(i) = arma::sum(xm.row(i)) * slice_sum;
        }
        double wsum = arma::sum(weight);
        if (wsum > 0) weight /= wsum;
        
        std::vector<bool> valid(I, false);
        for (int i = 0; i < I; ++i) {
          bool cond;
          if (scale == "logit") {
            cond = (ps1(i) != 0.0) && !std::isinf(ps1(i)) &&
              (ps2(i) != 0.0) && !std::isinf(ps2(i)) &&
              !std::isnan(ps1(i)) && !std::isnan(ps2(i));
          } else {
            cond = (ps1(i) != 0.0) && !std::isinf(ps1(i)) &&
              (ps2(i) != 0.0) && !std::isinf(ps2(i)) &&
              (ps1(i) != 1.0) && (ps2(i) != 1.0) &&
              !std::isnan(ps1(i)) && !std::isnan(ps2(i));
          }
          cond = cond &&
            (suma_rows(0,i) > 0.0) && (suma_rows(0,i) < 1.0) &&
            (suma_cols(0,i) > 0.0) && (suma_cols(0,i) < 1.0);
          valid[i] = cond;
        }
        
        std::vector<double> tr_rows, tr_cols, wv;
        for (int i = 0; i < I; ++i) {
          if (!valid[i]) continue;
          double tr1, tr2;
          if (scale == "logit") {
            tr1 = std::log(ps1(i));
            tr2 = std::log(ps2(i));
          } else {
            NumericVector qr = qnorm_fn(ps1(i));
            NumericVector qc = qnorm_fn(ps2(i));
            tr1 = qr[0]; tr2 = qc[0];
          }
          tr_rows.push_back(tr1);
          tr_cols.push_back(tr2);
          wv.push_back(weight(i));
        }
        
        double cor_e = 0.0;
        if (!tr_rows.empty()) {
          double wsum_v = 0.0;
          for (double w : wv) wsum_v += w;
          
          double mean1 = 0.0, mean2 = 0.0;
          for (size_t s = 0; s < tr_rows.size(); ++s) {
            mean1 += wv[s] * tr_rows[s];
            mean2 += wv[s] * tr_cols[s];
          }
          if (wsum_v > 0) { mean1 /= wsum_v; mean2 /= wsum_v; }
          
          double var1 = 0.0, var2 = 0.0, cov12 = 0.0;
          for (size_t s = 0; s < tr_rows.size(); ++s) {
            double d1 = tr_rows[s] - mean1;
            double d2 = tr_cols[s] - mean2;
            var1  += wv[s] * d1 * d1;
            var2  += wv[s] * d2 * d2;
            cov12 += wv[s] * d1 * d2;
          }
          if (wsum_v > 0) { var1 /= wsum_v; var2 /= wsum_v; }
          double denom = std::sqrt(var1 * var2);
          cor_e = (denom > 0) ? cov12 / denom : 0.0;
        }
        
        for (int i = 0; i < I; ++i) {
          if (!valid[i] || pjk_array(jj, kk, i) == 0.0) continue;
          NumericVector q1v = qnorm_fn(suma_rows(0, i));
          NumericVector q2v = qnorm_fn(suma_cols(0, i));
          NumericVector res = pnorm2d_fn(q1v[0], q2v[0], Named("rho") = cor_e);
          pjk_array(jj, kk, i) = res[0];
        }
        
      } // kk
    } // jj
    
    // Clamp a cero
    for (auto& v : pjk_array) if (v < 0.0) v = 0.0;
    
    // Marginales estimados
    arma::mat p_est_rows(I, J, arma::fill::zeros);
    arma::mat p_est_cols(I, K, arma::fill::zeros);
    for (int i = 0; i < I; ++i) {
      for (int j = 0; j < J; ++j) {
        for (int k = 0; k < K; ++k) {
          p_est_rows(i, j) += pjk_array(j, k, i);
          p_est_cols(i, k) += pjk_array(j, k, i);
        }
      }
    }
    
    // Reescalado columnas de referencia (ref1_0, kk)
    for (int kk : parties2) {
      for (int i = 0; i < I; ++i) {
        if (p_est_cols(i, kk) != 0.0) {
          pjk_array(ref1_0, kk, i) *= p_obs_cols(i, kk) / p_est_cols(i, kk);
        }
      }
    }
    
    // Reescalado filas de referencia (jj, ref2_0)
    for (int jj : parties1) {
      for (int i = 0; i < I; ++i) {
        if (p_est_rows(i, jj) != 0.0) {
          pjk_array(jj, ref2_0, i) *= p_obs_rows(i, jj) / p_est_rows(i, jj);
        }
      }
    }
    
    cont++;
    suma0 = suma1;
    suma1 = arma::accu(pjk_array);
    dif   = std::abs(suma0 - suma1);
    
  } // while
  
  // ---- IPF final por slice --------------------------------------------------
  for (int i = 0; i < I; ++i) {
    arma::mat slice_mat = pjk_array.slice(i) * arma::sum(ym.row(i));
    
    NumericMatrix slice_r = arma_to_rcpp(slice_mat);
    NumericVector y_i(K), x_i(J);
    for (int k = 0; k < K; ++k) y_i[k] = ym(i, k);
    for (int j = 0; j < J; ++j) x_i[j] = xm(i, j);
    
    NumericMatrix ipf_result = IPF_cpp(slice_r, y_i, x_i, tol);
    arma::mat ipf_mat = rcpp_to_arma(ipf_result);
    pjk_array.slice(i) = ipf_mat;
    
    double dif_r = 0.0, dif_c = 0.0;
    for (int j = 0; j < J; ++j) {
      dif_r += std::abs(xm(i,j) - arma::sum(ipf_mat.row(j)));
    }
    for (int k = 0; k < K; ++k) {
      dif_c += std::abs(ym(i,k) - arma::sum(ipf_mat.col(k)));
    }
    
    if (std::max(dif_r, dif_c) > 0.01) {
      NumericMatrix ipf2_result = IPF2_cpp(ipf_result, y_i, x_i);
      pjk_array.slice(i) = rcpp_to_arma(ipf2_result);
    }
    
  }
  
  // ---- Salida --------------------------------------------------------------
  NumericVector out_array(J * K * I);
  for (int i = 0; i < I; ++i) {
    for (int k = 0; k < K; ++k) {
      for (int j = 0; j < J; ++j) {
        out_array[j + J*k + J*K*i] = pjk_array(j, k, i);
      }
    }
  }
  
  return List::create(
    Named("vjk.array") = out_array,
    Named("iter")      = cont
  );
}
