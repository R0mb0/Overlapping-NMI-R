# Calculating the Overlapping Normalized Mutual Information (NMI) between two covers
# Robust and compatible with lists of character vectors (each a community)

get_membership_matrix <- function(communities, all_nodes) {
  mat <- matrix(0, nrow=length(all_nodes), ncol=length(communities))
  rownames(mat) <- all_nodes
  for (j in seq_along(communities)) {
    idx <- match(communities[[j]], all_nodes)
    idx <- idx[!is.na(idx)]
    if (length(idx) > 0) {
      mat[idx, j] <- 1
    }
  }
  mat
}

NMI <- function(cover1, cover2) {
  all_nodes <- sort(unique(c(unlist(cover1), unlist(cover2))))
  X <- get_membership_matrix(cover1, all_nodes)
  Y <- get_membership_matrix(cover2, all_nodes)
  n <- length(all_nodes)
  safe_log2 <- function(x) ifelse(x > 0, log2(x), 0)
  cond_entropy <- function(A, B) {
    kA <- ncol(A)
    kB <- ncol(B)
    H <- 0
    for (i in 1:kA) {
      minH <- Inf
      for (j in 1:kB) {
        Nij <- sum(A[,i] & B[,j])
        if (Nij == 0) next
        Ni <- sum(A[,i])
        Nj <- sum(B[,j])
        pij <- Nij / n
        pi <- Ni / n
        pj <- Nj / n
        Hij <- 0
        if (pij > 0 && (pi * pj) > 0)
          Hij <- Hij - (pij) * safe_log2(pij / (pi * pj))
        if ((Ni-Nij) > 0 && (pi * (1-pj)) > 0)
          Hij <- Hij - ((Ni-Nij)/n) * safe_log2(((Ni-Nij)/n) / (pi*(1-pj)))
        if ((Nj-Nij) > 0 && ((1-pi)*pj) > 0)
          Hij <- Hij - ((Nj-Nij)/n) * safe_log2(((Nj-Nij)/n) / ((1-pi)*pj))
        if ((n-Ni-Nj+Nij) > 0 && ((1-pi)*(1-pj)) > 0)
          Hij <- Hij - ((n-Ni-Nj+Nij)/n) * safe_log2(((n-Ni-Nj+Nij)/n) / ((1-pi)*(1-pj)))
        if (!is.nan(Hij) && Hij < minH)
          minH <- Hij
      }
      if (is.finite(minH)) H <- H + minH
    }
    H / kA
  }
  H_XY <- cond_entropy(X, Y)
  H_YX <- cond_entropy(Y, X)
  NMI_value <- 1 - 0.5 * (H_XY + H_YX)
  NMI_value
}