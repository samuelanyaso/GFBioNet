################################################################# combining results of parallel loop
comb <- function(x, ...) {
    lapply(seq_along(x), function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

################################################################# lambda sequence for GFBioNet
lambda.upper <- function(N, M, TT, alpha = 1) {
    Mstar <- M * (M - 1)  ## use this, if you do not adjust for double counting in actual calc. of theo. FDR
    # Mstar <- M*(M-1)/2 ## use this, if you adjusts double counting in actual calc. of theo. FDR
    p0 <- 0.5/(Mstar * TT)  # Bonferroni threshold (p-value)
    zscore <- qnorm(1 - p0/2)  # corresponding z-score
    res <- zscore/(alpha * sqrt(N))
    return(res)
}

lambda.lower <- function(N, M, TT, fdr, nsignals = 1000, alpha = 1) {
    temp <- (fdr * nsignals)/(2 * M * (M - 1) * (TT))  ## use this, if you do not adjust for double counting in actual calc. of theo. FDR
    # temp <- (fdr*nsignals)/(M*(M-1)*(TT)) ## use this, if you adjusts double counting in actual calc.
    # of theo. FDR
    temp <- qnorm(1 - temp)
    temp <- temp/(alpha * sqrt(N))
    return(temp)
}


################################################### define lambda for GGM
lambda.GGM.upper <- function(N, M, alpha = 1) {
    # Mstar <- M*(M-1)/2
    Mstar <- M * (M - 1)  # without correction for double counts
    p0 <- 0.5/(Mstar)  # Bonferroni threshold (p-value)
    zscore <- qnorm(1 - p0/2)  # corresponding z-score
    res <- zscore/(alpha * sqrt(N))
    return(res)
}

lambda.GGM.lower <- function(N, M, fdr, nsignals = 100, alpha = 1) {
    # temp <- (fdr*nsignals)/(M*(M-1))
    temp <- (fdr * nsignals)/(2 * M * (M - 1))  # without correction for double counts
    temp <- qnorm(1 - temp)
    temp <- temp/(alpha * sqrt(N))
    return(temp)
}



############################################################################### takes a vector
############################################################################### c('X1_X2_G1',
############################################################################### 'X2_X1_G1', 'X1_X2_G3')
############################################################################### returns unique edge-gene
############################################################################### interaction:
############################################################################### c('X1_X2_G1',
############################################################################### 'X1_X2_G3')
foo11 <- function(aa) {
    aa1 <- data.frame(do.call("rbind", strsplit(aa, "_")))
    aa1[, 1:2] <- apply(aa1[, 1:2], 2, function(x) as.numeric(substr(x, 2, nchar(x))))
    aa1$dum1 <- ifelse(aa1[, 1] < aa1[, 2], aa1[, 1], aa1[, 2])
    aa1$dum2 <- ifelse(aa1[, 2] > aa1[, 1], aa1[, 2], aa1[, 1])
    aa1[, 1] <- paste("X", aa1$dum1, sep = "")
    aa1[, 2] <- paste("X", aa1$dum2, sep = "")
    aa1$dum1 <- aa1$dum2 <- NULL
    aa1 <- apply(aa1, 1, paste0, collapse = "_")
    return(unique(aa1))
}

############################################################################### FDR calculated via
############################################################################### permutation may not
############################################################################### behave well when the
############################################################################### number of true signals
############################################################################### detected is small. Use
############################################################################### this function to ensure
############################################################################### the calculated FDR based
############################################################################### on increasing values of
############################################################################### lambda, is
############################################################################### non-increasing.

# Iterate through the vector to make it non-increasing
f00b <- function(vec) {

    is_non_increasing <- function(x) {
        all(diff(x) <= 0)
    }

    # Check if the vector is non-increasing
    if (!is_non_increasing(vec)) {
        # ensure the vector to be non-increasing
        for (i in 2:length(vec)) {
            if (vec[i] > vec[i - 1]) {
                vec[i] <- vec[i - 1]
            }
        }
    }
    return(vec)
}


############################################################################### Evaluate and process the
############################################################################### results for the triplets
############################################################################### Note: signals from BC
############################################################################### version returns both
############################################################################### main & interaction
############################################################################### effects However signals
############################################################################### from permutation version
############################################################################### returns only interaction
############################################################################### effects The function
############################################################################### handles both cases


triplets_process <- function(signals, cxnames, cgnames, sym = "or") {
    Xnames <- cxnames$dummy
    Gnames <- cgnames$dummy

    betas <- mapply(function(x, y) {
        if (!(is.null(x) | length(x) == 0)) {
            remp = data.frame(xi = y, xj = names(x), value = x)
            is_int = vapply(strsplit(remp$xj, "_"), length, integer(1L)) == 2L
            remp = remp[is_int, ]  ## keep only interactions
        } else {
            remp = NULL
        }
        return(remp)
    }, signals, Xnames, SIMPLIFY = F)
    betas <- do.call(rbind, betas)
    rownames(betas) <- NULL
    temp <- do.call(rbind, sapply(betas$xj, strsplit, split = "_"))
    betas <- data.frame(xi = betas$xi, xj = temp[, 1], G = temp[, 2], value = betas$value)



    if (nrow(betas) > 0) {
        betas$XX1 <- as.numeric(sapply(betas$xi, function(x) substr(x, 2, nchar(x))))
        betas$XX2 <- as.numeric(sapply(betas$xj, function(x) substr(x, 2, nchar(x))))
        betas$X1 <- ifelse(betas$XX1 < betas$XX2, betas$XX1, betas$XX2)
        betas$X2 <- ifelse(betas$XX2 > betas$XX1, betas$XX2, betas$XX1)
        betas$xi <- paste("X", betas$X1, sep = "")
        betas$xj <- paste("X", betas$X2, sep = "")
        betas$triplet <- apply(betas[, c("xi", "xj", "G")], 1, function(x) paste(x, collapse = "_", sep = ""))
        betas <- betas[order(betas$triplet), ]
        betas[, c("XX1", "XX2", "X1", "X2")] <- NULL
        rownames(betas) <- NULL

        ## For duplicate triplets
        betas <- by(betas, betas$triplet, function(x) {
            if (nrow(x) > 1) {
                xx <- x$value
                abs_val <- abs(xx)
                if (sym == "and") {
                  xx <- xx[1] * (abs_val[1] < abs_val[2]) + xx[2] * (abs_val[2] < abs_val[1])  ## AND criterion
                } else {
                  xx <- xx[1] * (abs_val[1] >= abs_val[2]) + xx[2] * (abs_val[2] >= abs_val[1])  ## OR criterion
                }
                x <- x[1, , drop = FALSE]
                x$value <- xx
            }
            return(x)
        }, simplify = F)
        betas <- do.call(rbind, betas)
        betas$triplet <- NULL
        rownames(betas) <- NULL



        ## update true variable names
        betas$xi <- cxnames$tnames[match(betas$xi, cxnames$dummy)]
        betas$xj <- cxnames$tnames[match(betas$xj, cxnames$dummy)]
        betas$G <- cgnames$tnames[match(betas$G, cgnames$dummy)]
        colnames(betas) <- c("trait1", "trait2", "genomic_factor", "value")
    } else {
        betas <- NULL
    }

    return(betas)
}
