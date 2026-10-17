# Normalised Laplacian of adjacency A.
# n = 1: unnormalised L = D - A
# n = 2: symmetric    L = I - D^{-1/2} A D^{-1/2}   (default)
# n = 3: random-walk  L = I - D^{-1} A
#
# Isolated nodes (deg = 0) are guarded so divisions stay finite.
LAPnorm <- function(A, n = 2) {
  if (!(n %in% c(1, 2, 3))) n <- 2

  deg <- rowSums(A)

  if (n == 1) {
    L <- -A
    diag(L) <- deg

  } else if (n == 2) {
    # L[i,j] = -A[i,j] / sqrt(deg[i] * deg[j]); diagonal set to 1
    d <- 1 / sqrt(pmax(deg, 1e-12))
    L <- (-A * outer(d, d))
    diag(L) <- 1

  } else if (n == 3) {
    # L[i,j] = -A[i,j] / deg[i]; diagonal set to 1
    L <- (-A / pmax(deg, 1e-12))
    diag(L) <- 1
  }

  return(L)
}
