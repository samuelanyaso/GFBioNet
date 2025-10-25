#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::NumericMatrix fastElemMult(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {
  int n = X.nrow();
  int p = X.ncol();
  int q = Y.ncol();
  
  Rcpp::NumericMatrix Z(n, p * q);
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      for (int k = 0; k < q; k++) {
        Z(i, j * q + k) = X(i, j) * Y(i, k);
      }
    }
  }
  
  return Z;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix fastScale(Rcpp::NumericMatrix X) {
  int n = X.nrow();
  int p = X.ncol();
  
  Rcpp::NumericMatrix Z(n, p);
  Rcpp::NumericVector col;
  
  for (int j = 0; j < p; j++) {
    // Extract column j
    col = X(Rcpp::_, j);
    
    // Center column
    double colmean = Rcpp::mean(col);
    Z(Rcpp::_, j) = col - colmean;
    
    // Scale column
    col = Z(Rcpp::_, j);
    double colstd = Rcpp::sd(col);
    Z(Rcpp::_, j) = col / colstd;
  }
  
  return Z;
}

// Correlation estimates
// x: data matrix
// Returns: Correlation matrix
// [[Rcpp::export]]
arma::mat cor2(const arma::mat& X) {
  arma::mat X_centered = X.each_row() - arma::mean(X, 0); // Center the matrix
  arma::mat X_std = X_centered.each_row() / arma::stddev(X, 0, 0); // Standardize
  arma::mat r = (X_std.t() * X_std) / (X.n_rows - 1); // Compute correlation matrix
  return r; 
}


// Compute correlation matrix with diagonal set to zero
arma::mat fast_correlation(const arma::mat& X) {
  arma::mat X_centered = X.each_row() - arma::mean(X, 0); // Center the matrix
  arma::mat X_std = X_centered.each_row() / arma::stddev(X, 0, 0); // Standardize
  arma::mat r = (X_std.t() * X_std) / (X.n_rows - 1); // Compute correlation matrix
  return r; 
}


// [[Rcpp::export]]
Rcpp::DataFrame getcor_fast(const arma::mat& X, const arma::mat& G,
                            const double cor_thresX, const double cor_thresG, const double cor_thresXG,
                            const std::vector<std::string>& xnames,
                            const std::vector<std::string>& gnames) {
  int M = X.n_cols;
  int TT = G.n_cols;
  int p = M * TT;
  
  // Precompute correlation matrices
  arma::mat corX = fast_correlation(X);
  arma::mat corG = fast_correlation(G);
  
  // Generate XG names
  std::vector<std::string> XGnames(p);
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < TT; ++j) {
      XGnames[i * TT + j] = xnames[i] + "_" + gnames[j];
    }
  }
  
  std::vector<std::string> XG1_all, XG2_all;
  std::vector<double> cor_all;
  std::vector<int> XG1ind_all, XG2ind_all;
  
  // 1. Correlated X
  for (int i = 0; i < M; ++i) {
    for (int j = i + 1; j < M; ++j) {
      double cx = corX(i, j);
      if (std::abs(cx) >= cor_thresX) {
        for (int t = 0; t < TT; ++t) {
          int idx1 = i * TT + t;
          int idx2 = j * TT + t;
          XG1_all.push_back(XGnames[idx1]);
          XG2_all.push_back(XGnames[idx2]);
          cor_all.push_back(cx);
          XG1ind_all.push_back(idx1 + 1);
          XG2ind_all.push_back(idx2 + 1);
        }
      }
    }
  }
  
  // 2. Correlated G
  for (int i = 0; i < TT; ++i) {
    for (int j = i + 1; j < TT; ++j) {
      double cg = corG(i, j);
      if (std::abs(cg) >= cor_thresG) {
        for (int m = 0; m < M; ++m) {
          int idx1 = m * TT + i;
          int idx2 = m * TT + j;
          XG1_all.push_back(XGnames[idx1]);
          XG2_all.push_back(XGnames[idx2]);
          cor_all.push_back(cg);
          XG1ind_all.push_back(idx1 + 1);
          XG2ind_all.push_back(idx2 + 1);
        }
      }
    }
  }
  
  // 3. Combined correlated X and G
  for (int xi = 0; xi < M; ++xi) {
    for (int xj = xi + 1; xj < M; ++xj) {
      double cx = corX(xi, xj);
      if (std::abs(cx) >= cor_thresX) {
        for (int gi = 0; gi < TT; ++gi) {
          for (int gj = gi + 1; gj < TT; ++gj) {
            double cg = corG(gi, gj);
            if (std::abs(cg) >= cor_thresG) {
              double corr = cx * cg;
              if (std::abs(corr) >= cor_thresXG) {
                int idx1 = xi * TT + gi;
                int idx2 = xj * TT + gj;
                XG1_all.push_back(XGnames[idx1]);
                XG2_all.push_back(XGnames[idx2]);
                cor_all.push_back(corr);
                XG1ind_all.push_back(idx1 + 1);
                XG2ind_all.push_back(idx2 + 1);
              }
            }
          }
        }
      }
    }
  }
  
  return Rcpp::DataFrame::create(
    Rcpp::Named("XG1") = XG1_all,
    Rcpp::Named("XG2") = XG2_all,
    Rcpp::Named("cor") = cor_all,
    Rcpp::Named("XG1ind") = XG1ind_all,
    Rcpp::Named("XG2ind") = XG2ind_all
  );
}



// [[Rcpp::export]]
Rcpp::List getcor_fast3(const std::vector<int> XG1ind_all,
                          const std::vector<int> XG2ind_all,
                        const int p) {
  
  // Determine correlated and uncorrelated predictors
  std::unordered_set<int> corpred_set;
  for (size_t i = 0; i < XG1ind_all.size(); ++i) {
    corpred_set.insert(XG1ind_all[i]);
    corpred_set.insert(XG2ind_all[i]);
  }
  
  std::vector<int> corpred(corpred_set.begin(), corpred_set.end());
  std::sort(corpred.begin(), corpred.end());
  
  std::vector<int> uncorpred;
  for (int i = 1; i <= p; ++i) {
    if (!corpred_set.count(i)) {
      uncorpred.push_back(i);
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("corpred") = corpred,
    Rcpp::Named("uncorpred") = uncorpred
  );
  
}
