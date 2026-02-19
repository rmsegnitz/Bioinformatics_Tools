#######################################################################################
# Implementation of parallel processing within limma's core normalization functions.
# Author: Max Segnitz, msegnitz@uw.edu
########################################################################################


# Modified from limma  
# Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). 
# “limma powers differential expression analyses for RNA-sequencing and microarray studies.” 
# Nucleic Acids Research, 43(7), e47. doi:10.1093/nar/gkv007.

# Started November 2025

# © Richard M Segnitz 2025 (modified from limma)
# License: This software is licensed under GNU General Public License, and may
# be modified and/or redistributed under GNU GPL version 3 or later. License details
# can be found in the accompanying this script, or at  (http://www.gnu.org/licenses/).
# 


# This was developed via chatGPT and was described as follows

# The aim is to add parallelism are inside the gene-by-gene array-weights routine 
# (.arrayWeightsGeneByGene) and the per-gene parts of voom() 
# (the parts that compute per-gene fits / residuals). There are two important algorithmic 
# constraints:
  
# The original .arrayWeightsGeneByGene uses a progressive (gene-by-gene) update: it updates the
# information matrix and the parameter vector gam as each gene is processed. That makes a simple 
# mclapply replacing the loop not numerically equivalent, because the weights aw change within 
# the loop and those updated aw values influence later gene fits.

# However, the mathematics can be rearranged to compute each gene’s contribution 
# (the small matrices / vectors that are added to the global info2 and dl) independently, then 
# sum those contributions and update gam. Doing that in a repeated outer-iteration 
# (recompute contributions using the current aw, aggregate, update gam, recompute, ...) 
# will produce an algorithm that is parallelizable across genes while preserving the 
# iterative nature of the original algorithm. This is what I implement below.


# I implemented a parallel-aware replacement for the gene-by-gene routine that:
#   Computes per-gene contributions in parallel 
# (using BiocParallel::bplapply if available; falls back to parallel::mclapply or plain lapply),
# Aggregates them,
# Updates gam and aw,
# Repeats until convergence (or maxiter is reached).
# 
# I also provide a voomWithQualityWeights_parallel() wrapper that uses the existing 
# voom() but calls the parallel arrayWeights-style routine so you get parallel 
# array-weight estimation inside the same overall flow.



## Parallel-friendly .arrayWeightsGeneByGene replacement
.arrayWeightsGeneByGene_parallel <- function(E, design=NULL, weights=NULL, var.design=NULL,
                                             prior.n=10, BPPARAM=NULL, maxiter=50L, tol=1e-5, trace=FALSE) {
  ngenes <- nrow(E)
  narrays <- ncol(E)
  if(is.null(design)) design <- matrix(1, narrays, 1)
  nparams <- ncol(design)
  
  if(is.null(var.design)) {
    Z2 <- contr.sum(narrays)
  } else {
    Z2 <- var.design
  }
  Z <- cbind(1, Z2)
  ngam <- ncol(Z2)
  
  # initial gam, aw
  gam <- rep_len(0, ngam)
  aw <- rep_len(1, narrays)
  info_prior <- prior.n * crossprod(Z2) # prior part of info2
  
  # helper: per-gene contribution given current aw
  per_gene_fun <- function(i, E_loc, design_loc, weights_loc, Z_loc, Z2_loc, aw_loc) {
    # returns list(info_adj (matrix ngam x ngam), dl (vector ngam), ok logical)
    y <- E_loc[i, ]
    if(is.null(weights_loc)) w_full <- aw_loc else w_full <- aw_loc * weights_loc[i, ]
    # check NAs
    if(anyNA(y)) {
      obs <- is.finite(y)
      nobs <- sum(obs)
      if(nobs <= 2L) return(list(info_adj = matrix(0, ncol(Z2_loc), ncol(Z2_loc)),
                                 dl = rep(0, ncol(Z2_loc)), ok = FALSE))
      X <- design_loc[obs, , drop=FALSE]
      ysub <- y[obs]
      wsub <- w_full[obs]
      fit <- lm.wfit(X, ysub, wsub)
      if(fit$df.residual < 2L) return(list(info_adj = matrix(0, ncol(Z2_loc), ncol(Z2_loc)),
                                           dl = rep(0, ncol(Z2_loc)), ok = FALSE))
      d <- rep(0, length(w_full)); d[obs] <- wsub * fit$residuals^2
      h1 <- rep(0, length(w_full)); h1[obs] <- 1 - hat(fit$qr)
    } else {
      fit <- lm.wfit(design_loc, y, w_full)
      d <- w_full * fit$residuals^2
      h1 <- 1 - hat(fit$qr)
      if(fit$df.residual < 2L) return(list(info_adj = matrix(0, ncol(Z2_loc), ncol(Z2_loc)),
                                           dl = rep(0, ncol(Z2_loc)), ok = FALSE))
    }
    s2 <- mean(fit$effects[-(1:fit$rank)]^2)
    if(s2 < 1e-15) return(list(info_adj = matrix(0, ncol(Z2_loc), ncol(Z2_loc)),
                               dl = rep(0, ncol(Z2_loc)), ok = FALSE))
    # info contribution: crossprod(Z, h1 * Z)
    info_full <- crossprod(Z_loc, h1 * Z_loc)
    # compute info adjustment that corresponds to original code:
    # info2 update used: info2 <- info2 + info[-1,-1] - (info[-1,1]/info[1,1]) %*% info[1,-1]
    # where info is crossprod(Z, h1*Z)
    info11 <- info_full[1, 1]
    if(info11 == 0) {
      # numerical guard
      info_adj <- matrix(0, ncol(Z2_loc), ncol(Z2_loc))
    } else {
      info_1m <- info_full[1, -1, drop = FALSE]       # 1 x ngam
      info_m1 <- info_full[-1, 1, drop = FALSE]       # ngam x 1
      info_mm <- info_full[-1, -1, drop = FALSE]      # ngam x ngam
      info_adj <- info_mm - (info_m1 / info11) %*% info_1m
    }
    # z = d/s2 - h1
    z <- d / s2 - h1
    dl <- crossprod(Z2_loc, z)  # vector length ngam
    list(info_adj = info_adj, dl = as.numeric(dl), ok = TRUE)
  }
  
  # choose parallel backend
  use_biocpar <- FALSE
  if(!is.null(BPPARAM)) {
    if(requireNamespace("BiocParallel", quietly = TRUE)) use_biocpar <- TRUE
  }
  
  if(use_biocpar) {
    bapply <- function(X, FUN) BiocParallel::bplapply(X, FUN, BPPARAM=BPPARAM)
  } else {
    # fallback: use mclapply if available (unix), else lapply
    if(.Platform$OS.type == "unix" && requireNamespace("parallel", quietly = TRUE)) {
      bapply <- function(X, FUN) parallel::mclapply(X, FUN, mc.preschedule = TRUE, mc.cores = getOption("mc.cores", 2L))
    } else {
      bapply <- function(X, FUN) lapply(X, FUN)
    }
  }
  
  # Outer iterative loop: recompute contributions using current aw, update gam, repeat
  for(iter in seq_len(maxiter)) {
    if(trace) cat("arrayWeights parallel: iteration", iter, "\n")
    # a single pass: compute contributions for all genes in parallel
    idx <- seq_len(ngenes)
    per_gene_wrapper <- function(i) per_gene_fun(i, E, design, weights, Z, Z2, aw)
    lst <- bapply(idx, per_gene_wrapper)
    
    # aggregate
    # initialize sums
    dl_sum <- rep(0, ngam)
    info_sum <- matrix(0, ngam, ngam)
    for(res in lst) {
      if(!is.null(res) && isTRUE(res$ok)) {
        dl_sum <- dl_sum + res$dl
        info_sum <- info_sum + res$info_adj
      }
    }
    
    info2_total <- info_prior + info_sum
    # solve for new gam. Guard against singular
    # use tryCatch to fall back to pseudoinverse if needed
    gam_new <- tryCatch({
      drop(solve(info2_total, dl_sum))
    }, error = function(e) {
      warning("info matrix singular in parallel arrayWeights; using solve with tol (svd) fallback")
      svd_inv <- function(M) {
        s <- svd(M)
        d <- s$d
        tol.svd <- max(dim(M)) * max(d) * .Machine$double.eps
        d_inv <- ifelse(d > tol.svd, 1/d, 0)
        s$v %*% (diag(d_inv, length(d_inv)) %*% t(s$u))
      }
      as.numeric(svd_inv(info2_total) %*% dl_sum)
    })
    
    # check convergence on gam
    if(max(abs(gam_new - gam)) < tol) {
      gam <- gam_new
      aw <- drop(exp(Z2 %*% (-gam)))
      if(trace) cat("Converged at iteration", iter, "\n")
      break
    }
    # update and continue
    gam <- gam_new
    aw <- drop(exp(Z2 %*% (-gam)))
    if(iter == maxiter && trace) cat("Reached maxiter without full convergence\n")
  }
  
  # return array weights
  names(aw) <- colnames(E)
  aw
}

## Wrapper that mirrors arrayWeights(..., method="genebygene") but calls parallel implementation
arrayWeights_parallel <- function(object, design=NULL, weights=NULL, var.design=NULL, var.group=NULL,
                                  prior.n=10, BPPARAM=NULL, maxiter=50L, tol=1e-5, trace=FALSE) {
  y <- getEAWP(object)
  E <- y$exprs
  ngenes <- nrow(E)
  narrays <- ncol(E)
  
  if(ngenes < 2L) return(rep(1, narrays))
  
  if(is.null(design)) {
    design <- matrix(1, narrays, 1)
    p <- 1L
  } else {
    design <- as.matrix(design)
    QR <- qr(design)
    p <- QR$rank
    if(p < ncol(design)) design <- design[, QR$pivot[1:p], drop = FALSE]
  }
  
  if(narrays - p < 2L) return(rep(1, narrays))
  
  if(is.null(weights)) weights <- y$weights
  
  if(!is.null(var.group)) {
    var.group <- droplevels(as.factor(var.group))
    contrasts(var.group) <- contr.sum(levels(var.group))
    var.design <- model.matrix(~var.group)
    var.design <- var.design[,-1, drop = FALSE]
  }
  
  # call parallel gene-by-gene
  .arrayWeightsGeneByGene_parallel(E, design=design, weights=weights, var.design=var.design,
                                   prior.n=prior.n, BPPARAM=BPPARAM, maxiter=maxiter, tol=tol, trace=trace)
}

## Top-level voomWithQualityWeights that uses the parallel arrayWeights
voomWithQualityWeights_parallel <- function(counts, design=NULL, lib.size=NULL, normalize.method="none",
                                            plot=FALSE, span=0.5, var.design=NULL, var.group=NULL,
                                            method="genebygene", maxiter=50, tol=1e-5, trace=FALSE,
                                            BPPARAM=NULL, ...) {
  if(plot) {
    oldpar <- par(mfrow=c(1,2))
    on.exit(par(oldpar))
  }
  
  # initial voom
  v <- voom(counts, design=design, lib.size=lib.size, normalize.method=normalize.method,
            plot=FALSE, span=span, ...)
  
  # choose arrayWeights implementation
  if(method == "genebygene") {
    aw <- arrayWeights_parallel(v, design=design, weights=NULL, var.design=var.design, var.group=var.group,
                                prior.n=10, BPPARAM=BPPARAM, maxiter=maxiter, tol=tol, trace=trace)
  } else {
    # for reml or other methods fall back to existing arrayWeights
    aw <- arrayWeights(v, design=design, weights=NULL, var.design=var.design, var.group=var.group,
                       prior.n=10, method=method, maxiter=maxiter, tol=tol, trace=trace)
  }
  
  # re-run voom with array weights
  v <- voom(counts, design=design, weights=aw, lib.size=lib.size, normalize.method=normalize.method,
            plot=plot, span=span, ...)
  
  # update array weights with new v (again parallel if chosen)
  if(method == "genebygene") {
    aw <- arrayWeights_parallel(v, design=design, weights=NULL, var.design=var.design, var.group=var.group,
                                prior.n=10, BPPARAM=BPPARAM, maxiter=maxiter, tol=tol, trace=trace)
  } else {
    aw <- arrayWeights(v, design=design, weights=NULL, var.design=var.design, var.group=var.group,
                       prior.n=10, method=method, maxiter=maxiter, tol=tol, trace=trace)
  }
  
  v$weights <- t(aw * t(v$weights))
  v$targets$sample.weights <- aw
  
  if(plot) {
    col <- NULL
    barplot(aw, names=1:length(aw), main="Sample-specific weights", ylab="Weight", xlab="Sample", col=col)
    abline(h=1, col=2, lty=2)
  }
  
  v
}


#################################
#  Notes & caveats to document. #
#################################

# 1) Algorithmic behaviour & numerical differences: The parallel function uses outer iterations that 
# recompute per-gene contributions in parallel and update gam in each iteration. This is not verbatim 
# identical to the original progressive single-pass update, but it retains the same statiscal objective 
# (the same information and score contributions are being aggregated) and converges toward a similar 
# solution. In practice it should produce nearly identical weights, but you must validate on your datasets.
# 
# 2) Convergence: We implemented an iterative outer-loop with maxiter and tol. You can increase maxiter 
# or tighten tol to match original behavior more closely. The trace=TRUE option helps monitor progress.
# 
# 3) Parallel backend: BiocParallel::MulticoreParam() on Unix/macOS; use BiocParallel::SnowParam() on 
# Windows (BiocParallel is cross-platform and recommended).
# If BiocParallel is not available, the code falls back to parallel::mclapply where supported, or lapply. 
# You should install & use BiocParallel for best results.
# 
# 4) Memory: The parallel version briefly holds per-gene results in memory for the aggregation step; 
# for extremely large datasets that may increase memory usage. If memory is tight, use a Snow backend with 
# chunked jobs or reduce the number of workers.



#################################
#    Next steps for testing     #
#################################

# Testing: Test on a small dataset first.
# 
# 1) Compare results with original voomWithQualityWeights():
#   
# 2) Compare v$targets$sample.weights from both versions (e.g., all.equal with tolerances).
# 
# 3) Confirm downstream lmFit/eBayes results are stable.
