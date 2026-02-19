# Compute covariate variances

# Takes a sample by feature matrix S of size f x n (eg expression matrix), and a design matrix D for n samples.
# Reports the variance explained by each covariate individually, as well as cumulatively as ordered in D


designVariance_explained <- function(S, D, groups, center = TRUE) {
    stopifnot(ncol(S) == nrow(D))
    stopifnot(length(groups) == ncol(D))
    
    if (center) {
      S <- S - rowMeans(S)
    }
    
    total_var <- sum(S^2)
    group_names <- unique(groups)
    
    indiv <- numeric(length(group_names))
    cumul <- numeric(length(group_names))
    
    for (i in seq_along(group_names)) {
      cols <- which(groups == group_names[i])
      Dj <- D[, cols, drop = FALSE]
      Sj_hat <- qr.fitted(qr(Dj), t(S))
      indiv[i] <- sum(Sj_hat^2) / total_var
    }
    term_order<-order(indiv, decreasing = T)
    group_names<-
    group_names[term_order]
    
    for (i in seq_along(group_names)) {
      cols <- which(groups %in% group_names[1:i])
      Dj <- D[, cols, drop = FALSE]
      Sj_hat <- qr.fitted(qr(Dj), t(S))
      cumul[i] <- sum(Sj_hat^2) / total_var
    }
    
    data.frame(
      term = group_names,
      var_explained_individual = indiv[term_order],
      var_explained_cumulative = cumul
    )
  }

