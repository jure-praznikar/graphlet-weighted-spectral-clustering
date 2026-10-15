
# Closure: A is captured from outer environment !!!!!
energy_function <- function(par) {
  weights <- c(1, par)
  ## sparsity_penalty
  total_w <- sum(weights)
  # k controls how steep the transition is at 4. 
  # k = 5 makes it look almost like a sharp step function.
  k <- 5
  alpha <- (1 / (1 + exp(-k * (total_w - 4)))) #4
  # sparsity_penalty <- alpha * sum((weights))
  P <- alpha * sum((weights)) # sparsity_penalty

  
  ## weighted sum of matrices
  AOPT <- 0
  for (i in seq_along(weights)) {
    AOPT <- AOPT + weights[i] * A[,,i]
  }
  # normalize
  m <- max(AOPT)
  if (m > 0) AOPT <- AOPT / m
  ## Normlasied Laplacian
  L<-LAPnorm(AOPT,2) # Symmetric normalized: L = I - D^{-1/2} A D^{-1/2}
  EVmax<-20 # specified number of eigenvalues and eigenvectors
  ev <- eigs_sym(L, k = EVmax, which = "SM") # "SM"  The k eigenvalues with smallest magnitude.
  lambdas <- sort(ev$values)[1:EVmax]

  raw_gaps <- diff(lambdas)
  raw_gaps[1] <- NA  # Force the algorithm to ignore lambda_2 - lambda_1 (forces G >= 2)
  G <- which.max(raw_gaps) # index of largest gap
  gap <- raw_gaps[G] # value of larges gape
  #
  beta<-1
  # Spectral separation
  S <- log( gap*G^(beta) / (mean(lambdas[1:G])+1.0E-6) + 1)
  #this is the same >> S <- log( gap*G^(beta+1) / (sum(lambdas[1:G])+1.0E-6) + 1)
  
  # Vertical Shift
  T <- log( mean(lambdas0[1:G])/(mean(lambdas[1:G])+1.0E-6) )
  
  # final function / for S and T log is used !!! can happen
  # that eogen values are ver small and gap/egen_value become large number
  # compared to P (sparsity), so log transform is used
  x <- S + T - P
  return(-x) # hjn is "minimise by default" but we search for MAX >> minus sign
}

