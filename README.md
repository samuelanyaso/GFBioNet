# GFBioNet (Genomic Factors regulating BIOlogical NETworks)

Identify **genomic factors that regulate a biological network**.

Given two data modalities—**X** (biological traits; e.g., genes, taxa, metabolites) and **G** (genomic factors; e.g., SNPs)—**GFBioNet**:

1. learns a **Gaussian Graphical Model (GGM)** over traits (edges = conditional dependencies), and
2. finds **trait–trait–genomic** triplets $(X_i, X_j, G_k)$ indicating which genomic factors modulate which network edges.

A key contribution is **robust false discovery control** that remains valid **even with strongly correlated predictors**, via a fast bias-correction scheme and node-wise regression along a regularization path (lasso / elastic net).

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("samuelanyaso/GFBioNet")
```

---

## What the package does

* **Stage 1 – Network learning (traits only):**
  `fitGGM()` estimates a sparse GGM on X using node-wise regression and controls FDR across the path to select a tuning parameter.
* **Stage 2 – Regulatory mapping (traits × genomic factors):**
  `GFBioNet()` tests trait–genomic main effects and **trait×genomic** interaction terms *conditioned on the learned network*, corrects for correlation-induced bias (lasso or elastic net), and returns discoveries as **triplets** $(X_i, X_j, G_k)$ with controlled FDR.

Both stages return a **path of solutions**, an **estimated FDR curve**, and the **optimal index** selected by the FDR cutoff.

---

## Data requirements

* `X`: numeric matrix $N \times M$ (rows = samples; columns = traits).
* `G`: numeric matrix $N \times T$ (rows = samples; columns = genomic factors).
* No missing values. Centering/scaling is recommended.

---

## Quick start (simulated example)

```r
library(GFBioNet)

X     # matrix for traits
G     # matrix for genomic factors

## -----------------------------
## Stage 1: learn the trait GGM
## -----------------------------
fit1 <- fitGGM(
  X            = X,
  nsignals     = 1000,
  nlambda      = 50,
  alpha        = 1,     # 1 = lasso, (0,1) = elastic net
  qfdr         = 0.10,    # FDR cap for stage 1
  r1           = 0.4,     # high-corr threshold for thinning
  ncores       = max(1L, parallel::detectCores() - 1L)
)

# Selected tuning parameter
fit1$opt.index
# Undirected adjacency (0/1) of the learned network
A <- fit1$path.fit

## --------------------------------------------
## Stage 2: map genomic factors onto the edges
## --------------------------------------------
fit2 <- GFBioNet(
  X            = X,
  G            = G,
  basenet      = A,        # from Stage 1
  nsignals2     = 10e4,
  nlambda      = 50,
  alpha        = 1,      # must match modeling choice
  qfdr2        = 0.10,     # FDR cap for stage 2
  r1           = 0.5,
  topL_mu      = 3,
  ncores       = max(1L, parallel::detectCores() - 1L)
)

# Selected tuning parameter
fit2$opt.index

# Detected signals
triplets <- fit2$fit
head(triplets)
```

**Interpretation:** Each row of `triplets` is a detected **(trait1, trait2, genomic factor)** regulator with its estimated coefficient at the chosen $\lambda$.

---
