##' Graph-guided, bias-corrected FDR control for \eqn{X \times G} interactions
##'
##' \code{GFBioNet} fits a Gaussian graphical model via node-wise regression
##' and estimates interaction effects between biological traits \eqn{X} and genomic 
##' factors \eqn{G}. The procedure uses correlation screening to build an interaction
##' graph, performs penalized regression (lasso or elastic net) at a path of
##' \eqn{\lambda} values, corrects for bias induced by correlated predictors,
##' and estimates the expected number of false positives to control the FDR.
##'
##' @title Fast bias-corrected node-wise regression with FDR control for \eqn{X \times G}
##' @param X A numeric matrix of size \eqn{N \times M} (samples \eqn{\times} traits).
##' @param G A numeric matrix of size \eqn{N \times T} (samples \eqn{\times} genomic factors).
##'   Interactions \eqn{X \times G} are formed internally (element-wise products by column).
##' @param standardize Logical flag for standardization of the columns of \code{X} and \code{G}, 
##' prior to fitting the model sequence. Default: \code{FALSE}.
##' @param basenet Optional \eqn{M \times M} adjacency (0/1 or -/+) matrix from an initial
##'   Gaussian graphical model. If \code{NULL} (default), the base network is estimated
##'   internally using node-wise penalized regression and symmetrized by \code{sym}.
##'   The diagonal is required to be zero.
##' @param nsignals1 Integer. If \code{basenet} is \code{NULL} To estimate the base network (stage 1),
##'   the rough upper bound on the number of non-null signals used to construct the \eqn{\lambda} range. 
##'   Default: \code{100}.
##' @param nsignals2 Integer; rough upper bound on the number of non-null signals used to
##'   construct the \eqn{\lambda} range. Default: \code{1000}.
##' @param nlambda Integer; number of penalty values in the path. Default: \code{30}.
##' @param lambda.min,lambda.max Optional numeric bounds for the penalty path. If either
##'   is \code{NULL}, bounds are computed by \code{lambda.lower}/\code{lambda.upper}
##'   heuristics given \code{N}, \code{M}, \code{T}, \code{alpha}, and \code{nsignals}.
##' @param alpha Elastic-net mixing parameter passed to \pkg{glmnet}: \code{alpha = 1}
##'   gives lasso; \code{0 < alpha < 1} gives elastic net (ridge component included).
##'   Bias is estimated accordingly (lasso and EN formulas are handled internally).
##' @param sym Character; symmetrization rule for the undirected base network:
##'   \code{"or"} (default) uses the union of directed edges, \code{"and"} uses the
##'   intersection.
##' @param qfdr1 Numeric in (0,1); target FDR level for the base-network stage (stage 1).
##'   Default: \code{0.1}.
##' @param qfdr2 Numeric in (0,1); target FDR level for the interaction stage (stage 2). Default: \code{0.1}.
##' @param r1 Numeric in (0,1); “strong-correlation” threshold used to thin highly
##'   correlated detections and to define the exclusion set when computing bias-adjusted 
##'   false-positive counts. Default: \code{0.5}.
##' @param r2 Numeric in (0,1); reserved for stricter correlation handling. Default: \code{0.7}.
##' @param topL_mu Integer; number of strongest neighbors (by absolute correlation) used
##'   to accumulate the lasso bias term \eqn{\mu} per detected predictor. For elastic net,
##'   bias uses the analytical \eqn{\delta} factor and does not require top-\eqn{L}
##'   truncation. Default: \code{3}.
##' @param ncores Integer; number of CPU cores for parallel node-wise fits. Defaults to
##'   \code{max(1, parallel::detectCores() - 1)}.
##'
##' @details
##' \strong{Algorithm (high level)}:
##' \enumerate{
##'   \item Fit/obtain a base network over \eqn{X} (Gaussian graphical model via node-wise
##'         penalized regression), then symmetrize by \code{sym}.
##'   \item Build the interaction design \eqn{X \times G} and a correlation screening graph
##'         combining: (i) \eqn{X}-only correlations replicated across \eqn{G}, (ii) \eqn{G}-only
##'         correlations replicated across \eqn{X}, and (iii) cross terms.
##'   \item For each node \eqn{i} and for each \eqn{\lambda}:
##'     \enumerate{
##'       \item Fit \pkg{glmnet} (lasso/EN) including \emph{main effects} (neighboring \eqn{X}
##'             and all \eqn{G}) and \emph{interactions} pre-screened to neighbors of \eqn{i}.
##'       \item Thin highly correlated detections using \code{r1} (keep stronger by absolute
##'             coefficient).
##'       \item Compute expected false positives in two pools (uncorrelated vs correlated),
##'             applying bias correction \eqn{\mu}:
##'             \itemize{
##'               \item Lasso: \eqn{\mu_z = \lambda \, \mathrm{sign}(\beta_t)\, \mathrm{cor}(z,t)}.
##'               \item Elastic net: \eqn{\delta_t = \frac{(1 + \hat{b}_{OLS,t})\, \alpha \lambda}{1 + \alpha \lambda}},\;
##'                     \eqn{\mu_z = \delta_t \cdot \mathrm{cor}(z,t)} with \eqn{\hat{b}_{OLS,t}}
##'                     computed from the pre-thinning residual for exact parity with the baseline.
##'             }
##'       \item Accumulate the expected false positives and detections to estimate the FDR
##'             along the path and select an optimal \eqn{\lambda}.
##'     }
##' }
##'
##' The implementation keeps the baseline math intact while using cached adjacency lists,
##' vectorized variance updates, and sparse triplet assembly for speed.
##'
##' @return A named \code{list} with:
##' \item{fit}{Signals detected.}
##' \item{signals}{Length-\code{nlambda} list of sparse \eqn{M \times M} coefficient
##'   matrices (directed node-wise fits after correlation-thinning at each \eqn{\lambda}).}
##' \item{basenet}{The symmetrized base adjacency used for node-wise interaction screening.}
##' \item{lambda}{Numeric vector of length \code{nlambda} with penalty values.}
##' \item{opt.index}{Index of the selected \eqn{\lambda} based on the FDR curve.}
##' \item{dfres}{Data frame with columns: \code{lambda}, \code{exp.false} (expected false
##'   positives), \code{detections}, and \code{fdr}.}
##' \item{beta.fit}{Sparse matrix: coefficients at the selected \eqn{\lambda}.}
##' \item{path.fit}{Absolute-value adjacency used to form the undirected graph by
##'   \code{sym}.}
##'
##' @section Notes:
##' \itemize{
##'   \item \code{alpha = 1} (lasso) and \code{0 < alpha < 1} (elastic net) are supported.
##'   \item Columns of \code{X} and \code{G} should be on comparable scales if effect-size
##'         comparability is desired; interactions are standardized internally.
##'   \item The diagonal of \code{basenet} must be zero; if \code{basenet = NULL} the function
##'         estimates it internally (recommended for end-to-end runs).
##' }
##'
##' @author Samuel Anyaso-Samuel
##' @export
##' @useDynLib GFBioNet
##' @import Matrix
##' @import glmnet
##' @import methods
##' @import RcppArmadillo
##' @importFrom matrixStats colVars
##' @importFrom stats pnorm qnorm sd var
##' @importFrom parallel detectCores
##' @importFrom methods as
##' @importFrom Rcpp sourceCpp
GFBioNet <- function(
    X, G,
    standardize=FALSE,
    basenet   = NULL,
    nsignals1 = 100,
    nsignals2  = 1000,
    nlambda   = 30,
    lambda.min=NULL, lambda.max=NULL,
    alpha     = 1,
    sym       = "or",
    qfdr1     = 0.1,
    qfdr2     = 0.1,
    r1        = 0.5,
    r2        = 0.7,
    topL_mu   = 3,  
    ncores    = max(1L, parallel::detectCores(logical = TRUE) - 1L)
) {

  # -------------------------------
  # Shapes & canonical column names
  # -------------------------------
  N <- nrow(X)
  M <- ncol(X)
  TT <- ncol(G)
  stopifnot(M > 1, TT >= 1, N >= 2)
  
  if(standardize){
    X <- fastScale(X)
    G <- fastScale(G)
  }
  
  ## define variable names
  Xnames <- paste0("X", seq_len(M))
  Gnames <- paste0("G", seq_len(TT))
  if(is.null(colnames(X))){
    cxnames <- data.frame(tnames=Xnames, dummy=Xnames)
  } else {
    cxnames <- data.frame(tnames=colnames(X), dummy=Xnames)
  }
  colnames(X) <- Xnames
  
  if(is.null(colnames(G))){
    cgnames <- data.frame(tnames=Gnames, dummy=Gnames)
  } else {
    cgnames <- data.frame(tnames=colnames(G), dummy=Gnames)
  }
  colnames(G) <- Gnames
  
  
  # X-major scaffold which DEFINES the column order in XG
  xreps   <- rep(Xnames, each = TT)     # "X1","X1",...,"X2","X2",...
  greps   <- rep(Gnames, times = M)     # "G1","G2",...,"GTT","G1",...
  XGnames <- paste(xreps, greps, sep = "_")
  
  # Global index space (1..p) that we use EVERYWHERE
  p        <- M * TT
  idx_of   <- function(j, k) (j - 1L) * TT + k
  j_of     <- function(idx) ((idx - 1L) %/% TT) + 1L
  k_of     <- function(idx) ((idx - 1L) %%  TT) + 1L
  name_of  <- function(idx) paste0(Xnames[j_of(idx)], "_", Gnames[k_of(idx)])
  
  # Quick sanity: the first TT columns must be X1_G1..X1_GTT and map via idx_of(1,k)
  stopifnot(identical(XGnames[idx_of(1L, seq_len(TT))], paste0("X1_", Gnames)))
  
  # --------------------------------------
  # Step 1) Base network (Gaussian GGM)
  # --------------------------------------
  if (is.null(basenet)) {
    ggm <- fitGGM(
      X = X, alpha = alpha, sym = sym,
      nsignals = nsignals1, qfdr = qfdr1, nlambda = nlambda,
      lambda.min = NULL, lambda.max = NULL,
      ncores = max(1L, ncores %/% 2L)
    )
    basenet <- ggm$path.fit
  }
  if (sum(Matrix::diag(basenet)) != 0) stop("Diagonal of adjacency matrix should be 0.")
  base_adj <- basenet
  Matrix::diag(base_adj) <- 0L
  
  # --------------------------------------
  # Step 2) Build full XG once (X-major)
  # --------------------------------------
  XG <- fastElemMult(X, G)   # N x (M*TT), already column-major in R
  XG <- fastScale(XG)        # standardize for glmnet / correlations
  colnames(XG) <- XGnames               # enforce the X-major naming
  
  # ------------------------------------------------------------------
  # Step 3) Interaction-correlation screening graph (your exact logic)
  #   1) X-only edges       (replicate each (j1,j2) across k=1..TT)
  #   2) G-only edges       (replicate each (k1,k2) across j=1..M)
  #   3) cross edges        (j1,k1) -- (j2,k2) with weight corX*corG
  # ------------------------------------------------------------------
  cor.thresX  <- max(0.35, r1 - 0.05)
  cor.thresG  <- cor.thresX
  cor.thresXG <- cor.thresX
  
  corX <- cor2(X); diag(corX) <- 0
  corG <- cor2(G); diag(corG) <- 0
  
  ## 3.1) X-only
  cX_idx <- which(abs(corX) >= cor.thresX & upper.tri(corX), arr.ind = TRUE)
  if (!is.matrix(cX_idx)) cX_idx <- matrix(cX_idx, ncol = 2)
  j1v <- as.integer(cX_idx[, 1]); j2v <- as.integer(cX_idx[, 2])
  cXw <- corX[cbind(j1v, j2v)]
  
  if (length(cXw)) {
    kvec        <- rep(seq_len(TT),   times = length(cXw)) # k cycles fastest
    j1rep       <- rep(j1v,           each  = TT)
    j2rep       <- rep(j2v,           each  = TT)
    Xonly_i     <- idx_of(j1rep, kvec)
    Xonly_j     <- idx_of(j2rep, kvec)
    Xonly_w     <- rep(cXw,           each  = TT)
  } else {
    Xonly_i <- Xonly_j <- integer(0); Xonly_w <- numeric(0)
  }
  
  ## 3.2) G-only
  cG_idx <- which(abs(corG) >= cor.thresG & upper.tri(corG), arr.ind = TRUE)
  if (!is.matrix(cG_idx)) cG_idx <- matrix(cG_idx, ncol = 2)
  k1v <- as.integer(cG_idx[, 1]); k2v <- as.integer(cG_idx[, 2])
  cGw <- corG[cbind(k1v, k2v)]
  
  if (length(cGw)) {
    jvec        <- rep(seq_len(M),    times = length(cGw)) # j cycles fastest
    k1rep       <- rep(k1v,           each  = M)
    k2rep       <- rep(k2v,           each  = M)
    Gonly_i     <- idx_of(jvec, k1rep)
    Gonly_j     <- idx_of(jvec, k2rep)
    Gonly_w     <- rep(cGw,           each  = M)
  } else {
    Gonly_i <- Gonly_j <- integer(0); Gonly_w <- numeric(0)
  }
  
  ## 3.3) Cross part (exact product rule), streamed to limit memory
  Cross_i <- Cross_j <- integer(0); Cross_w <- numeric(0)
  if (length(cXw) && length(cGw)) {
    block <- max(256L, 8192L %/% max(1L, length(cGw)))
    for (st in seq.int(1L, length(cXw), by = block)) {
      en <- min(length(cXw), st + block - 1L); bx <- st:en
      mask <- (abs(cXw[bx]) %o% abs(cGw)) > cor.thresXG
      if (!any(mask)) next
      sel <- which(mask, arr.ind = TRUE)
      r   <- bx[sel[, 1]]; s <- sel[, 2]
      Cross_i <- c(Cross_i, idx_of(j1v[r], k1v[s]))
      Cross_j <- c(Cross_j, idx_of(j2v[r], k2v[s]))
      Cross_w <- c(Cross_w,  cXw[r] * cGw[s])
    }
  }
  
  # Collate, keep upper-triangle, apply cor.thresXG
  XG1ind <- c(Xonly_i, Gonly_i, Cross_i)
  XG2ind <- c(Xonly_j, Gonly_j, Cross_j)
  Ww     <- c(Xonly_w, Gonly_w, Cross_w)
  keep   <- which(XG1ind < XG2ind & abs(Ww) > cor.thresXG)
  XG1ind <- XG1ind[keep]; XG2ind <- XG2ind[keep]; Ww <- Ww[keep]
  
  # For auditing (optional): the table in your original string format
  corXG <- data.frame(
    XG1    = XGnames[XG1ind],
    XG2    = XGnames[XG2ind],
    cor    = Ww,
    XG1ind = XG1ind,
    XG2ind = XG2ind,
    row.names = NULL, stringsAsFactors = FALSE
  )
  
  # Correlated vs uncorrelated predictors (GLOBAL SETS)
  preds.vec  <- seq_len(p)
  corpred0   <- sort(unique(c(XG1ind, XG2ind)))
  uncorpred0 <- setdiff(preds.vec, corpred0)
  
  # Build adjacency lists (GLOBAL, symmetric). We never slice corXG again later.
  make_adj <- function(XG1ind, XG2ind, Ww, p, r1) {
    src <- c(XG1ind, XG2ind)
    dst <- c(XG2ind, XG1ind)
    wts <- c(Ww,     Ww)
    o   <- order(src)
    src <- src[o]; dst <- dst[o]; wts <- wts[o]
    bins <- split(seq_along(src), src)
    
    any_idx <- vector("list", p)
    any_w   <- vector("list", p)
    r1_idx  <- vector("list", p)
    r1_w    <- vector("list", p)
    
    for (s in names(bins)) {
      s <- as.integer(s)
      ii <- bins[[as.character(s)]]
      d  <- dst[ii]; w <- wts[ii]
      any_idx[[s]] <- d; any_w[[s]] <- w
      m <- abs(w) >= r1
      r1_idx[[s]] <- if (any(m)) d[m] else integer(0)
      r1_w[[s]]   <- if (any(m)) w[m] else numeric(0)
    }
    list(any_idx = any_idx, any_w = any_w, r1_idx = r1_idx, r1_w = r1_w)
  }
  adj <- make_adj(XG1ind, XG2ind, Ww, p, r1)
  
  ###################################################
  ## Compute lambda range
  if(is.null(lambda.min) | is.null(lambda.max)){
    lambda1 <- lambda.lower(N=N, M=M, TT=TT, fdr=0.2, nsignals=nsignals2, alpha=alpha)
    lambda2 <- lambda.upper(N=N, M=M, TT=TT, alpha=alpha)
  } else {
    lambda1 <- lambda.min
    lambda2 <- lambda.max
  }
  lambda <- seq(from=lambda1, to=lambda2, length.out=nlambda)
  lambda <- sort(lambda, decreasing=TRUE)
  
  
  
  # --------------------------------------
  # Step 4) Per-node fits (parallelizable)
  #   IMPORTANT: we keep **global ids** everywhere.
  #   When we need a design submatrix, we select columns by global ids directly.
  # --------------------------------------
  # runner   <- if (ncores > 1L && .Platform$OS.type != "windows") parallel::mclapply else lapply
  # node_ids <- seq_len(M)
  
  # -------- Parallel over nodes (columns)
  node_ids <- seq_len(M)
  
  worker_fun <- function(i) {
    xi <- X[, i]
    
    # Neighbor taxa of node i (integers 1..M)
    sig_taxa <- which(base_adj[, i] != 0L)
    M0 <- length(sig_taxa)
    if (M0 == 0L) {
      return(list(signals2 = vector("list", nlambda),
                  exp.fp3   = rep(0, nlambda),
                  detections3 = rep(0L, nlambda)))
    }
    
    # --- GLOBAL ids for interactions used in this node:
    #     For each j in sig_taxa, include all k=1..TT using the X-major indexer.
    preds.vec1 <- as.vector(unlist(lapply(sig_taxa, function(j) idx_of(j, seq_len(TT)))))
    
    # Double-check: column names match exactly the global selection
    stopifnot(identical(colnames(XG)[preds.vec1],
                        as.vector(unlist(lapply(sig_taxa, function(j) paste0(Xnames[j], "_", Gnames))))))
    
    # Local counts for FP partitioning (STILL GLOBAL ids)
    local_uncor_global <- intersect(uncorpred0, preds.vec1)
    local_cor_global   <- intersect(corpred0,  preds.vec1)
    n.uncorpred        <- length(local_uncor_global)
    n.corpred          <- length(local_cor_global)
    
    # Design matrix for glmnet: [neighbors of X_i], all G, and the selected interactions (in global order)
    X_design <- cbind(X[, sig_taxa, drop = FALSE], G, XG[, preds.vec1, drop = FALSE])
    
    fit <- glmnet::glmnet(
      x = X_design, y = xi, alpha = alpha,
      thresh = 1e-8, maxit = 10000,
      intercept = FALSE, standardize = FALSE,
      lambda = lambda
    )
    coef_path <- as.matrix(stats::coef(fit))[-1, , drop = FALSE]
    betas  <- coef_path[seq_len(M0), , drop = FALSE]
    deltas <- coef_path[M0 + seq_len(TT), , drop = FALSE]
    gammas <- coef_path[M0 + TT + seq_len(length(preds.vec1)), , drop = FALSE]  # NOTE: columns are in EXACT preds.vec1 order
    
    # Outputs per lambda
    expected.false3 <- numeric(nlambda)
    detections3     <- integer(nlambda)
    signals_path    <- vector("list", nlambda)
    
    for (k in seq_len(nlambda)) {
      beta.k   <- betas[,  k]
      delta.k  <- deltas[, k]
      gamma.k  <- gammas[, k]                          # aligned with preds.vec1
      # Embed into a GLOBAL coefficient vector length p (cheap & unambiguous)
      gamma.global <- numeric(p)
      gamma.global[preds.vec1] <- gamma.k
      
      # ---- (A) Thin correlated detections: keep the stronger within |corr| ≥ r1 ----
      # --- Keep-stronger among r1-correlated detections -----------------------
      nz_global <- which(gamma.global != 0)            # GLOBAL ids
      if (length(nz_global) <= 1L) {
        kept_global <- nz_global
      } else {
        
        ## estimate the correlation between the detected markers
        cemp1 <- cor2(XG[, nz_global, drop=FALSE])
        cemp2 <- apply(cemp1, 2, function(x) abs(x) >= r1);
        cemp3 <- which(cemp2 & upper.tri(cemp1), arr.ind = T)
        
        ## for each pair of correlated variables, retain the stronger, and drop the weaker.
        ## more than two detected markers may be highly correlated, but with low probability,
        ## because we are interested in large values of lambda.
        if(nrow(cemp3) > 0){
          
          remp <- as.numeric()
          for(ii in 1:nrow(cemp3)){
            temp00 <- nz_global[cemp3[ii,]]
            temp0 <- gamma.global[temp00]
            remp <- c(remp,  temp00[which.min(abs(temp0))])
          }
          kept_global <- setdiff(nz_global, remp)
          rm(temp0, temp00, remp)
        } else {
          kept_global <- nz_global
        }
      }
      detections3[k] <- length(kept_global)
      
      
      # ---- (B) Expected FPs: UNCORRELATED pool (use only GLOBAL ids) ----------
      exp_p1 <- 0
      if (n.uncorpred > 0L) {
        uncor_nz_global <- intersect(local_uncor_global, nz_global)
        uncor_z_global  <- setdiff(local_uncor_global, uncor_nz_global)
        
        y_main <- X[, sig_taxa, drop = FALSE] %*% beta.k + G %*% delta.k
        if (length(uncor_nz_global)) {
          y_res <- xi - y_main - XG[, uncor_nz_global, drop = FALSE] %*% gamma.global[uncor_nz_global]
        } else {
          y_res <- xi - y_main
        }
        sd_res <- stats::sd(y_res); if (sd_res <= 0) sd_res <- .Machine$double.eps
        
        # zeros in the uncorrelated pool
        if (length(uncor_z_global)) {
          exp_p1 <- 2 * pnorm(-lambda[k] * alpha*sqrt(N) / sd_res) * length(uncor_z_global)
        }
        # nonzeros (rank-1 variance updates; vectorized)
        if (length(uncor_nz_global)) {
          vy  <- as.numeric(stats::var(y_res))
          vx  <- matrixStats::colVars(XG[, uncor_nz_global, drop = FALSE])
          cy  <- as.numeric(crossprod(y_res - mean(y_res),
                                      scale(XG[, uncor_nz_global, drop = FALSE], center = TRUE, scale = FALSE))) / (N - 1)
          bj  <- gamma.global[uncor_nz_global]
          sd_upd <- sqrt(pmax(0, vy + bj^2 * vx + 2 * bj * cy))
          sd_upd[sd_upd <= 0] <- .Machine$double.eps
          exp_p1 <- exp_p1 + sum(2 * pnorm(-lambda[k] * alpha*sqrt(N) / sd_upd))
        }
      }
      
      # ---- (C) Expected FPs: CORRELATED pool (GLOBAL ids only) ----------------
      exp_p2 <- 0
      if (n.corpred > 0L) {
        # xrm = any GLOBAL feature (in this node's set) that is |corr|≥r1 to a kept NZ
        xrm_global <- integer(0)
        if (length(kept_global)) {
          for (tg in kept_global) {
            nb <- adj$r1_idx[[tg]]
            if (!is.null(nb) && length(nb)) xrm_global <- c(xrm_global, nb)
          }
          xrm_global <- unique(setdiff(xrm_global, kept_global))
          # keep only features that are part of this node's interaction set
          xrm_global <- intersect(xrm_global, preds.vec1)
        }
        
        ## the markers removed in first part of case 1
        xrm1 <- nz_global[which(!(nz_global %in% kept_global))]
        xrm_global <- sort(unique(c(xrm1,xrm_global)))
        n.xrm <- length(xrm_global)
        
        
        mu_global <- numeric(p)  # GLOBAL index vector of biases (same size as gamma.global)
        
        if (length(kept_global)) {
          for (tg in kept_global) {
            cand <- adj$any_idx[[tg]]
            if (is.null(cand) || !length(cand)) next
            
            # restrict to this node's interaction set and exclude tg & xrm already handled outside
            cand <- intersect(cand, preds.vec1)
            cand <- setdiff(cand, tg)
            if (!length(cand)) next
            
            # exact correlations; XG columns are standardized => cor = crossprod/(N-1)
            cor_exact <- as.numeric(crossprod(XG[, tg, drop = FALSE],
                                              XG[, cand, drop = FALSE])) / (N - 1)
            
            # LASSO: original behavior (top-L by |corr|, use λ * sign(γ_t))
            if (alpha == 1) {
              sgn  <- sign(gamma.global[tg])
              ord  <- order(-abs(cor_exact))
              take <- head(ord, topL_mu)
              if (length(take)) {
                mu_global[cand[take]] <- mu_global[cand[take]] + lambda[k] * sgn * cor_exact[take]
              }
            } else {
              # ELASTIC NET (alpha in (0,1)): compute δ_t = (1 + b_ols(t))*(αλ)/(1+αλ)
              # b_ols(t) computed efficiently via a single y_res_int:
              #   b_ols(t) = γ_t + mean(x_t * y_res_int) / mean(x_t^2)
              xt      <- XG[, tg]
              denom   <- mean(xt^2); # if (denom <= 0) denom <- .Machine$double.eps
              gamma_t <- gamma.global[tg]
              y_res_int <- xi - (XG[, -tg, drop = FALSE] %*% gamma.global[-tg])
              bols_t <- mean(xt * as.numeric(y_res_int)) / denom
              
              # if (is.null(y_res_int)) {
              #   # no selected interactions: y_res_int = xi
              #   bols_t <- gamma_t + mean(xt * xi) / denom
              # } else {
              #   bols_t <- gamma_t + mean(xt * as.numeric(y_res_int)) / denom
              # }
              A      <- alpha * lambda[k]
              delta  <- (1 + bols_t) * A / (1 + A)
              
              # IMPORTANT: for EN we do NOT use sign(gamma_t); δ already carries the effect.
              # Also, per Version 1, do not truncate to topL_mu under EN (keep all neighbors).
              mu_global[cand] <- mu_global[cand] + delta * cor_exact
            }
          }
        }
        # --------- END NEW μ (bias) block ----------
        
        
        
        # index.2 = GLOBAL features in this node's correlated set with coef OR bias nonzero
        index2_global <- union(intersect(local_cor_global, nz_global),
                               which(mu_global != 0))
        index2_global <- intersect(index2_global, preds.vec1)       # keep only this node's features
        
        # Residual for "both zero" terms: subtract main + ALL index2 contributions
        y_tmp <- xi - (X[, sig_taxa, drop = FALSE] %*% beta.k + G %*% delta.k)
        if (length(index2_global)) {
          y_tmp <- y_tmp - XG[, index2_global, drop = FALSE] %*% gamma.global[index2_global]
        }
        sd_tmp <- stats::sd(y_tmp); if (sd_tmp <= 0) sd_tmp <- .Machine$double.eps
        
        # remove xrm before the shifted tails; count of "both zero"
        index2_kept <- setdiff(index2_global, xrm_global)
        n.index2    <- length(index2_kept)
        exp_p2 <- exp_p2 + 2 * pnorm(-lambda[k] * alpha*sqrt(N) / sd_tmp) *
          max(0L, n.corpred - n.index2 - n.xrm)
        
        # shifted tails for coef or bias nonzero (vectorized)
        if (n.index2) {
          vy  <- as.numeric(stats::var(y_tmp))
          vx  <- matrixStats::colVars(XG[, index2_kept, drop = FALSE])
          cy  <- as.numeric(crossprod(y_tmp - mean(y_tmp),
                                      scale(XG[, index2_kept, drop = FALSE], center = TRUE, scale = FALSE))) / (N - 1)
          bj  <- gamma.global[index2_kept]
          sd_prime <- sqrt(pmax(0, vy + bj^2 * vx + 2 * bj * cy))
          sd_prime[sd_prime <= 0] <- .Machine$double.eps
          mu_vec <- mu_global[index2_kept]
          exp_p2 <- exp_p2 + sum(
            pnorm((alpha*lambda[k] - mu_vec) * sqrt(N) / sd_prime, lower.tail = FALSE) +
              pnorm(-(alpha*lambda[k] + mu_vec) * sqrt(N) / sd_prime, lower.tail = TRUE)
          )
        }
      }
      expected.false3[k] <- exp_p1 + exp_p2
      
      # ---- Pack outputs (main effects + kept interactions, all named) --------
      # main effects are those rows without "_" in the glmnet coef names
      nz_all  <- coef_path[, k, drop = TRUE]
      nz_main <- nz_all[nz_all != 0]
      is_main <- vapply(strsplit(names(nz_main), "_"), length, integer(1L)) == 1L
      
      kept_int_ids <- kept_global
      kept_int     <- gamma.global[kept_int_ids]
      names(kept_int) <- vapply(kept_int_ids, name_of, character(1))
      
      signals_path[[k]] <- c(nz_main[is_main], kept_int)
    } # end lambda loop
    
    list(signals2 = signals_path,
         exp.fp3  = expected.false3,
         detections3 = detections3)
  }
  
  
  # Run in parallel
  node_results <- if (ncores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(node_ids, worker_fun, mc.cores = ncores)
  } else {
    lapply(node_ids, worker_fun)
  }
  
  
  
  ## replace NULL elements (cases with edges) with zeros
  roo <- function(XX, nlambda){
    XX <- lapply(XX, function(x) {
      if (is.null(x)) {
        rep(0, nlambda)
      } else {
        x
      }
    })
    return(XX)
  } 
  
  
  signals2_temp <- sapply(node_results, "[[",  1, simplify = F)
  
  expected.false2 <- sapply(node_results, "[[", 2, simplify = T)
  expected.false2 <- t(expected.false2)
  
  ps2 <- sapply(node_results, "[[", 3, simplify = T)
  ps2 <- t(ps2)
  
  signals2 <- vector("list",length = M)
  signals2 <- replicate(nlambda, signals2, simplify = FALSE)
  
  for(i in 1:M){
    
    sig.taxa <- Xnames[which(basenet[, i] != 0)]
    M0 <- length(sig.taxa)
    if(M0 == 0) next
    
    ## signals detected based on original (non permuted) data
    for(k in 1:nlambda){
      signals2[[k]][[i]] <- signals2_temp[[i]][[k]]
    } 
  }
  
  
  tmp_signals2 <- lapply(signals2, function(x1){
    r0 <- sapply(x1, function(x){
      if(length(x) > 0){
        t1 <- names(x)
        t0 <- sapply(strsplit(t1, "_"), function(y) length(y) > 1, simplify = T)
        t2 <- t1[t0]
        t2
      }
    })
  })
  
  ndetections2 <- sapply(tmp_signals2, function(x1){
    jj <- mapply(function(x, y){
      if(length(y) > 0)  paste(x, y, sep="_")
    }, Xnames, x1)
    jj <- unlist(jj)
    if(!is.null(jj)) {
      jj <- foo11(aa=jj) ## adjusts for double counting
    }
    return(length(jj))
  }, simplify = TRUE)
  
  #################################################################################
  ## estimates the FDR
  df <- data.frame(d=1:nlambda, lambda=lambda,
                   exp.false2=pmin(colSums(expected.false2), colSums(ps2)),
                   detections2=colSums(ps2),
                   detections4=ndetections2)
  
  ## FDR calculated after BC
  df$fdr2 <- df$exp.false2/pmax(1, df$detections2)
  
  
  ## Ensure the calculated FDR is monotonically decreasing
  df <- df[order(df$lambda, decreasing = F), ]
  df$fdr2 <- f00b(df$fdr2 )
  df <- df[order(df$lambda, decreasing = T), ]
  
  # utilize the optimal lambda based on FDR calculated using permutation
  opt.index <- max(min(which(df$fdr2 >= qfdr2))-1, 1)
  opt.lambda <- lambda[opt.index]
  
  
  ## Post-processing
  if(is.numeric(opt.index)){
    fit <- triplets_process(signals=signals2[[opt.index]], cxnames = cxnames,cgnames = cgnames, sym=sym)
  } else {
    fit <- NULL
  }
    
  
  res <- list(fit=fit, signals=signals2, basenet=basenet, lambda=lambda,
              opt.index = opt.index, dfres = df)
  return(res)
}
