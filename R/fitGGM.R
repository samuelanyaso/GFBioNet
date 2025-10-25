##' Fast bias-corrected node-wise regression for Gaussian graphical models
##'
##' Estimate a Gaussian Graphical Model (GGM) via node-wise penalized regression
##' (lasso / elastic net) with a fast, bias-corrected procedure for counting
##' expected false positives and selecting \eqn{\lambda} by an estimated FDR curve.
##' Compared to a naïve implementation, this function caches correlations,
##' vectorizes variance updates, and assembles coefficient paths sparsely for speed.
##'
##' @title Estimate an undirected GGM by fast node-wise penalized regression with FDR control
##' @param X Numeric matrix of size \eqn{N \times M} (samples \eqn{\times} variables).
##'   Columns are the nodes of the GGM.
##' @param standardize Logical flag for standardization of the columns of \code{X}, 
##'   prior to fitting the model sequence. Default: \code{FALSE}.
##' @param nsignals Rough upper bound on the number of non-nulls used to set the
##'   \eqn{\lambda} range when \code{lambda.min} or \code{lambda.max} is \code{NULL}.
##'   Default \code{1000}.
##' @param nlambda Length of the penalty path. Default \code{30}.
##' @param lambda.min,lambda.max Optional numeric endpoints of the \eqn{\lambda} path.
##'   If either is \code{NULL}, the range is computed by internal heuristics
##'   \code{lambda.GGM.lower}/\code{lambda.GGM.upper} given \code{N}, \code{M},
##'   \code{alpha}, and \code{nsignals}.
##' @param alpha Elastic-net mixing parameter passed to \pkg{glmnet}:
##'   \code{alpha = 1} is lasso; \code{0 < alpha < 1} is elastic net.
##'   Default \code{1}.
##' @param sym Symmetrization rule for the undirected graph constructed from
##'   directed node-wise fits; \code{'or'} (union, default) or \code{'and'}
##'   (intersection).
##' @param qfdr Target FDR level used to pick the optimal \eqn{\lambda} along the
##'   path. Default \code{0.1}.
##' @param r1 Threshold (in absolute correlation) used to thin highly-correlated
##'   detections within a node-wise fit and to define the exclusion set when
##'   computing bias-adjusted false-positive counts. Default \code{0.4}.
##' @param corX Optional precomputed correlation matrix of \code{X} (diagonal will
##'   be zeroed). If \code{NULL}, it is computed internally. Supplying it can save time
##'   when calling the function repeatedly on the same data. Default \code{NULL}.
##' @param cor.thres Optional correlation screening threshold used to form
##'   per-node candidate sets; if \code{NULL}, set to \code{max(0.2, r1 - 0.2)}.
##' @param ncores Number of CPU cores for parallel node-wise regressions. Default
##'   \code{max(1, parallel::detectCores(logical = TRUE) - 1)}.
##'
##' @details
##' \strong{Algorithm (high level):}
##' \enumerate{
##' \item For each node \eqn{i}, regress \eqn{X_i} on the remaining columns \eqn{X_{-i}}
##'       using \pkg{glmnet} over a path of \eqn{\lambda} values (\code{alpha} controls
##'       lasso vs elastic net).
##' \item Within each \eqn{\lambda}, thin pairs of selected predictors whose absolute
##'       correlation exceeds \code{r1} by keeping the predictor with the larger
##'       absolute coefficient.
##' \item Split predictors into an “uncorrelated” pool and a “correlated” pool
##'       (based on \code{cor.thres}). For each pool, compute the expected number of
##'       false positives using Gaussian tail bounds. For the correlated pool, apply
##'       a bias term \eqn{\mu} to the tails:
##'       \itemize{
##'         \item Lasso: \eqn{\mu_z=\lambda\,\mathrm{sign}(\beta_t)\,\mathrm{cor}(z,t)}.
##'         \item Elastic net (\eqn{0<\alpha<1}): \eqn{\delta_t=\frac{(1+\hat b_{OLS,t})\,\alpha\lambda}{1+\alpha\lambda}},\;
##'               \eqn{\mu_z=\delta_t\,\mathrm{cor}(z,t)}, where \eqn{\hat b_{OLS,t}} is a
##'               simple OLS update using the current residual (fast rank-1 formulas).
##'       }
##' \item Sum the expected false positives across nodes to get \code{exp.false} at
##'       each \eqn{\lambda}, and compute \code{fdr = exp.false / max(1, detections)}.
##' \item Choose the largest \eqn{\lambda} such that the estimated FDR does not exceed
##'       \code{qfdr}, and build the undirected adjacency by \code{sym}.
##' }
##'
##' The implementation uses cached adjacency lists from \code{corX}, vectorized
##' variance/SD updates for rank-1 changes in residuals, and sparse triplet
##' assembly of coefficient matrices to reduce memory traffic.
##'
##' @return A list with components:
##' \item{lambda}{Numeric vector of \eqn{\lambda} values (length \code{nlambda}).}
##' \item{opt.index}{Index of the selected \eqn{\lambda} based on the FDR curve.}
##' \item{signals}{Length-\code{nlambda} list of sparse \eqn{M \times M} matrices
##'   containing directed node-wise coefficients at each \eqn{\lambda} after
##'   correlation thinning.}
##' \item{dfres}{Data frame with columns \code{lambda}, \code{exp.false},
##'   \code{detections}, and \code{fdr}.}
##' \item{beta.fit}{Sparse \eqn{M \times M} matrix of coefficients at the selected
##'   \eqn{\lambda}.}
##' \item{path.fit}{Absolute-value adjacency used to form the undirected graph by
##'   \code{sym}.}
##'
##' @author Samuel Anyaso-Samuel
##' @export
##' @useDynLib GFBioNet
##' @import Matrix
##' @import glmnet
##' @import methods
##' @importFrom matrixStats colVars
##' @importFrom stats pnorm qnorm sd var
##' @importFrom parallel detectCores
##' @importFrom Rcpp sourceCpp
fitGGM <- function(X, standardize = FALSE, nsignals = 1000, nlambda = 30, lambda.min = NULL, lambda.max = NULL,
    alpha = 1, sym = "or", qfdr = 0.1, r1 = 0.4, corX = NULL, cor.thres = NULL, ncores = max(1, parallel::detectCores(logical = TRUE) -
        1)) {

    N <- nrow(X)
    M <- ncol(X)

    if (standardize) {
        X <- fastScale(X)
    }

    ## define column names
    Xnames <- paste0("X", seq_len(M))
    if (is.null(colnames(X))) {
        cxnames <- data.frame(tnames = Xnames, dummy = Xnames)
    } else {
        cxnames <- data.frame(tnames = colnames(X), dummy = Xnames)
    }
    colnames(X) <- Xnames

    # default correlation threshold
    if (is.null(cor.thres))
        cor.thres <- max(0.2, r1 - 0.2)

    ################################################### define lambda for GGM (same as V1)
    if (is.null(lambda.min) | is.null(lambda.max)) {
        lambda1 <- lambda.GGM.lower(N = N, M = M, fdr = 0.2, nsignals = nsignals, alpha = alpha)
        lambda2 <- lambda.GGM.upper(N = N, M = M, alpha = alpha)
    } else {
        lambda1 <- lambda.min
        lambda2 <- lambda.max
    }
    lambda <- sort(seq(from = lambda2, to = lambda1, length.out = nlambda), decreasing = TRUE)

    # -------- Correlation: compute once, cache structure (same thresholds as V1)
    if (is.null(corX)) {
        corX <- cor2(X)
        diag(corX) <- 0
    } else {
        diag(corX) <- 0
    }
    cor_mask <- abs(corX) >= cor.thres  # for μ candidates
    r1_mask <- abs(corX) >= r1  # for “strong-corr” thinning and xrm

    # adjacency lists for fast lookups
    cor2_list <- lapply(seq_len(M), function(j) which(cor_mask[, j]))
    r1_list <- lapply(seq_len(M), function(j) which(r1_mask[, j]))

    # global partitions (exclude i per-node later)
    cor3X <- vapply(cor2_list, length, integer(1))
    uncorpred_global <- which(cor3X == 0L)
    corpred_global <- which(cor3X != 0L)

    # -------- Pre-allocate outputs
    exp.fal <- matrix(0, nrow = M, ncol = nlambda)
    ndetections <- matrix(0L, nrow = M, ncol = nlambda)

    # collect sparse coefficient triplets (kept after thinning) per lambda
    coef_triplets <- vector("list", nlambda)
    for (k in seq_len(nlambda)) coef_triplets[[k]] <- list(i = integer(0), j = integer(0), x = numeric(0))
    append_triplets <- function(lst, ii, jj, xx) {
        lst$i <- c(lst$i, ii)
        lst$j <- c(lst$j, jj)
        lst$x <- c(lst$x, xx)
        lst
    }

    # parallel scaffold
    node_ids <- seq_len(M)
    X <- as.matrix(X)

    worker_fun <- function(i) {
        xi <- X[, i]
        idx_minus_i <- setdiff(node_ids, i)
        Xtemp <- X[, idx_minus_i, drop = FALSE]

        # per-node partitions (mirror V1)
        uncorpred <- setdiff(uncorpred_global, i)
        corpred <- setdiff(corpred_global, i)
        n.unc <- length(uncorpred)
        n.cor <- length(corpred)

        # one glmnet fit for all lambdas
        fit <- glmnet::glmnet(x = Xtemp, y = xi, alpha = alpha, thresh = 1e-08, maxit = 10000, intercept = FALSE,
            standardize = FALSE, lambda = lambda)
        B <- as.matrix(stats::coef(fit))[-1, , drop = FALSE]  # (M-1) x nlambda

        res_exp_fal <- numeric(nlambda)
        res_ndet <- integer(nlambda)
        local_trip <- vector("list", nlambda)
        for (k in seq_len(nlambda)) local_trip[[k]] <- list(i = integer(0), j = integer(0), x = numeric(0))

        for (k in seq_len(nlambda)) {
            beta_k <- B[, k]
            if (all(beta_k == 0)) {
                # No selection at this lambda
                sd_xi <- stats::sd(xi)
                if (sd_xi <= 0)
                  sd_xi <- .Machine$double.eps
                expfal_p1 <- if (n.unc > 0)
                  2 * pnorm(-lambda[k] * alpha * sqrt(N)/sd_xi) * n.unc else 0
                expfal_p2 <- if (n.cor > 0)
                  2 * pnorm(-lambda[k] * alpha * sqrt(N)/sd_xi) * n.cor else 0
                res_exp_fal[k] <- expfal_p1 + expfal_p2
                res_ndet[k] <- 0L
                next
            }

            # pre-thinning nonzeros (est.nonzero.markers0 in V1)
            est_nz0 <- idx_minus_i[which(beta_k != 0)]

            # r1-thinning: keep stronger within |corr|>=r1 pairs
            if (length(est_nz0) > 1L) {
                sub_mask <- r1_mask[est_nz0, est_nz0, drop = FALSE]
                if (any(sub_mask & upper.tri(sub_mask, diag = FALSE))) {
                  pairs <- which(sub_mask & upper.tri(sub_mask, diag = FALSE), arr.ind = TRUE)
                  babs <- abs(beta_k[match(est_nz0, idx_minus_i)])
                  losers <- integer(0)
                  for (pp in seq_len(nrow(pairs))) {
                    a <- pairs[pp, 1]
                    b <- pairs[pp, 2]
                    losers <- c(losers, if (babs[a] < babs[b]) est_nz0[a] else est_nz0[b])
                  }
                  est_nz1 <- setdiff(est_nz0, unique(losers))
                } else {
                  est_nz1 <- est_nz0
                }
            } else {
                est_nz1 <- est_nz0
            }
            res_ndet[k] <- length(est_nz1)

            # store kept coefficients as triplets (sparse)
            if (length(est_nz1)) {
                ii <- est_nz1
                jj <- rep.int(i, length(ii))
                xx <- beta_k[match(ii, idx_minus_i)]
                local_trip[[k]]$i <- c(local_trip[[k]]$i, ii)
                local_trip[[k]]$j <- c(local_trip[[k]]$j, jj)
                local_trip[[k]]$x <- c(local_trip[[k]]$x, xx)
            }

            ## ---------- Part 1: Uncorrelated pool ---------- FIX #1: use the **pre-thinning** residual
            ## (V1 uses beta.temp before thinning)
            y_full <- if (length(est_nz0))
                xi - X[, est_nz0, drop = FALSE] %*% beta_k[match(est_nz0, idx_minus_i)] else xi
            y_full <- as.numeric(y_full)

            expfal_p1 <- 0
            if (n.unc > 0L) {
                # FIX #2: split using est_nz0 (not est_nz1) to mirror V1's `index.1`
                unc_nz <- intersect(uncorpred, est_nz0)
                unc_z <- setdiff(uncorpred, unc_nz)

                if (length(unc_z)) {
                  sd_full <- stats::sd(y_full)
                  if (sd_full <= 0)
                    sd_full <- .Machine$double.eps
                  expfal_p1 <- 2 * pnorm(-lambda[k] * alpha * sqrt(N)/sd_full) * length(unc_z)
                }
                if (length(unc_nz)) {
                  vy <- as.numeric(stats::var(y_full))
                  vx <- matrixStats::colVars(X[, unc_nz, drop = FALSE])
                  cy <- as.numeric(crossprod(y_full - mean(y_full), scale(X[, unc_nz, drop = FALSE], center = TRUE,
                    scale = FALSE)))/(N - 1)
                  bj <- beta_k[match(unc_nz, idx_minus_i)]
                  sd_upd <- sqrt(pmax(0, vy + bj^2 * vx + 2 * bj * cy))
                  sd_upd[sd_upd <= 0] <- .Machine$double.eps
                  expfal_p1 <- expfal_p1 + sum(2 * pnorm(-lambda[k] * alpha * sqrt(N)/sd_upd))
                }
            }

            ## ---------- Part 2: Correlated pool ----------
            expfal_p2 <- 0
            if (n.cor > 0L)
                {
                  # xrm as in V1: neighbors (|corr|>=r1) of **kept** detections, plus the losers
                  xrm <- integer(0)
                  if (length(est_nz1)) {
                    xrm <- unique(unlist(r1_list[est_nz1], use.names = FALSE))
                    xrm <- setdiff(xrm, c(i, est_nz1))
                  }
                  if (length(est_nz0)) {
                    lost <- setdiff(est_nz0, est_nz1)
                    if (length(lost))
                      xrm <- sort(unique(c(xrm, lost)))
                  }
                  n_xrm <- length(xrm)

                  # μ accumulation only for variables correlated (>= cor.thres) with kept detections,
                  # and **excluding {i} ∪ xrm** exactly like V1 does on est.nonzero.cors.
                  mu_all <- numeric(M)
                  if (length(est_nz1)) {
                    if (alpha == 1) {
                      bt <- beta_k[match(est_nz1, idx_minus_i)]
                      sgn <- sign(bt)
                      for (tt in seq_along(est_nz1)) {
                        t_idx <- est_nz1[tt]
                        z <- cor2_list[[t_idx]]
                        if (!length(z))
                          next
                        z <- setdiff(z, c(i, xrm))  # FIX #3: drop {i} ∪ xrm from bias targets
                        if (!length(z))
                          next
                        mu_all[z] <- mu_all[z] + lambda[k] * sgn[tt] * corX[z, t_idx]
                      }
                    } else if (alpha > 0 && alpha < 1) {
                      # EN bias with one-pass OLS update from the **pre-thinning** residual (y_full)
                      A <- alpha * lambda[k]
                      # Provide beta_t when computing b_ols(t) = beta_t + mean(y_full * x_t) /
                      # mean(x_t^2)
                      beta_full_t <- numeric(M)
                      beta_full_t[idx_minus_i] <- beta_k
                      for (t_idx in est_nz1) {
                        xt <- X[, t_idx]
                        denom <- mean(xt^2)
                        if (denom <= 0)
                          denom <- .Machine$double.eps
                        bols_t <- beta_full_t[t_idx] + mean(y_full * xt)/denom  # FIX #4: uses y_full (pre-thinning)
                        delta <- (1 + bols_t) * A/(1 + A)
                        z <- cor2_list[[t_idx]]
                        if (!length(z))
                          next
                        z <- setdiff(z, c(i, xrm))  # FIX #3 again for EN branch
                        if (!length(z))
                          next
                        mu_all[z] <- mu_all[z] + delta * corX[z, t_idx]
                      }
                    }
                  }

                  # FIX #5: index.2 must be built with **est_nz0** (pre-thinning), just like V1
                  nz_idx2 <- union(intersect(corpred, est_nz0), which(mu_all != 0))

                  # V1: y.temp = xi - X[, index.2] %*% beta[index.2] (do NOT drop xrm here)
                  y_tmp <- if (length(nz_idx2))
                    xi - X[, nz_idx2, drop = FALSE] %*% {
                      # beta at global positions: map back from (M-1)
                      bk <- numeric(M)
                      bk[idx_minus_i] <- beta_k
                      bk[nz_idx2]
                    } else xi
                  y_tmp <- as.numeric(y_tmp)

                  # Then remove xrm from the pool for counting shifted tails
                  nz_idx2_keep <- setdiff(nz_idx2, xrm)
                  n2 <- length(nz_idx2_keep)

                  # “both zero” count
                  if (n.cor - n2 - n_xrm > 0) {
                    sd_tmp <- stats::sd(y_tmp)
                    if (sd_tmp <= 0)
                      sd_tmp <- .Machine$double.eps
                    expfal_p2 <- expfal_p2 + 2 * pnorm(-lambda[k] * alpha * sqrt(N)/sd_tmp) * (n.cor - n2 -
                      n_xrm)
                  }

                  # shifted tails where coef!=0 or μ!=0
                  if (n2 > 0) {
                    vy <- as.numeric(stats::var(y_tmp))
                    vx <- matrixStats::colVars(X[, nz_idx2_keep, drop = FALSE])
                    cy <- as.numeric(crossprod(y_tmp - mean(y_tmp), scale(X[, nz_idx2_keep, drop = FALSE],
                      center = TRUE, scale = FALSE)))/(N - 1)
                    # bj is the *fitted coefficient* (zero for bias-only entries); map from (M-1) where
                    # needed
                    bj <- numeric(length(nz_idx2_keep))
                    map <- match(intersect(nz_idx2_keep, idx_minus_i), nz_idx2_keep)
                    if (length(map))
                      bj[map] <- beta_k[match(nz_idx2_keep[map], idx_minus_i)]

                    sd_prime <- sqrt(pmax(0, vy + bj^2 * vx + 2 * bj * cy))
                    sd_prime[sd_prime <= 0] <- .Machine$double.eps
                    mu_vec <- mu_all[nz_idx2_keep]

                    expfal_p2 <- expfal_p2 + sum(pnorm((alpha * lambda[k] - mu_vec) * sqrt(N)/sd_prime, lower.tail = FALSE) +
                      pnorm(-(alpha * lambda[k] + mu_vec) * sqrt(N)/sd_prime, lower.tail = TRUE))
                  }
                }  # end Part 2

            res_exp_fal[k] <- expfal_p1 + expfal_p2
        }  # end lambda loop

        list(expfal = res_exp_fal, ndet = res_ndet, triplets = local_trip)
    }  # end worker_fun

    # Run in parallel
    results <- if (ncores > 1L && .Platform$OS.type != "windows") {
        parallel::mclapply(node_ids, worker_fun, mc.cores = ncores)
    } else {
        lapply(node_ids, worker_fun)
    }

    # collect node-wise summaries + merge triplets per lambda
    for (i in seq_len(M)) {
        exp.fal[i, ] <- results[[i]]$expfal
        ndetections[i, ] <- results[[i]]$ndet
        for (k in seq_len(nlambda)) {
            tr <- results[[i]]$triplets[[k]]
            if (length(tr$x)) {
                coef_triplets[[k]] <- append_triplets(coef_triplets[[k]], tr$i, tr$j, tr$x)
            }
        }
    }

    # FDR & selection (unchanged)
    df <- data.frame(id = seq_len(nlambda), lambda = lambda, exp.false = pmin(colSums(exp.fal), colSums(ndetections)),
        detections = colSums(ndetections))
    df$fdr <- df$exp.false/pmax(1L, df$detections)
    opt.index <- max(min(which(df$fdr >= qfdr)) - 1L, 1L)
    opt.lambda <- lambda[opt.index]

    # build sparse coefficient matrices for each lambda from triplets
    signals <- vector("list", nlambda)
    for (k in seq_len(nlambda)) {
        tr <- coef_triplets[[k]]
        if (length(tr$x)) {
            signals[[k]] <- sparseMatrix(i = tr$i, j = tr$j, x = tr$x, dims = c(M, M), dimnames = list(Xnames,
                Xnames))
        } else {
            signals[[k]] <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(M, M), dimnames = list(Xnames,
                Xnames))
        }
    }

    ## update object names
    signals <- sapply(signals, function(x) {
        rownames(x) <- colnames(x) <- cxnames$tnames
        return(x)
    }, simplify = FALSE)

    beta.fit <- signals[[opt.index]]
    path.fit <- abs(beta.fit)
    if (sym == "and") {
        path.fit <- sign(path.fit * t(path.fit))
    } else if (sym == "or") {
        path.fit <- sign(path.fit + t(path.fit))
    }


    list(lambda = lambda, opt.index = opt.index, signals = signals, dfres = df, beta.fit = beta.fit, path.fit = path.fit)
}
