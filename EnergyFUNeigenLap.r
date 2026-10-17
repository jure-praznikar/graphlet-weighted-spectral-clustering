# Energy function F = S + T - P (see paper / README).
# Closure: matrices A and baseline spectrum lambdas0 are taken from the
# caller environment (LFR_Lap.r). hjn() minimises, so we return -F.
energy_function <- function(par) {
  # G0 weight is fixed at 1; par = (w1, ..., w8)
  weights <- c(1, par)

  # Soft sparsity penalty: alpha rises sharply once total weight exceeds 4
  # (logistic with steepness k = 5). P = alpha * sum(weights).
  total_w <- sum(weights)
  k <- 5
  alpha <- 1 / (1 + exp(-k * (total_w - 4)))
  P <- alpha * sum(weights)

  # Weighted sum of graphlet adjacency matrices, then scale to max = 1
  AOPT <- 0
  for (i in seq_along(weights)) {
    AOPT <- AOPT + weights[i] * A[, , i]
  }
  m <- max(AOPT)
  if (m > 0) AOPT <- AOPT / m

  # Symmetric normalised Laplacian spectrum (smallest EVmax eigenvalues)
  L <- LAPnorm(AOPT, 2)
  EVmax <- 20
  ev <- eigs_sym(L, k = EVmax, which = "SM")
  lambdas <- sort(ev$values)[1:EVmax]

  # Largest eigengap; force G >= 2 by discarding gap(lambda_2 - lambda_1)
  raw_gaps <- diff(lambdas)
  raw_gaps[1] <- NA
  G <- which.max(raw_gaps)
  gap <- raw_gaps[G]

  beta <- 1
  # S: large gap after a block of small eigenvalues (spectral separation)
  S <- log(gap * G^(beta) / (mean(lambdas[1:G]) + 1.0e-6) + 1)
  # Equivalent form:
  # S <- log(gap * G^(beta + 1) / (sum(lambdas[1:G]) + 1.0e-6) + 1)

  # T: eigenvalues lowered relative to the G0-only baseline lambdas0
  T <- log(mean(lambdas0[1:G]) / (mean(lambdas[1:G]) + 1.0e-6))

  # Log on S and T keeps them on a comparable scale to P when eigenvalues
  # are very small (otherwise gap / mean(lambda) would dominate).
  x <- S + T - P
  return(-x)  # minimise -F  <=>  maximise F
}
