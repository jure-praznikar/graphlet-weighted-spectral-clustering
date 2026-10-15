##LAPnorm <- function(A,n) {
##      #check if is 1,2,3, default is 2 symmetric normalised
##      if (!(n %in% c(1,2,3))) n<-2    
##      ## Degree matrix
##      deg <- rowSums(A)
##      D <- diag(deg)
##      if (n==1) {
##        ## 1. Unnormalized Laplacian
##        L <- D - A
##      } else if (n==2) {    
##        # symmetric normalised
##        D_inv_sqrt <- diag(1 / sqrt(deg))
##        L <- diag(nrow(A)) - D_inv_sqrt %*% A %*% D_inv_sqrt
##     } else if (n==3) {     
##        ## Random-walk normalized Laplacian
##        D_inv <- diag(1 / deg)
##        L <- diag(nrow(A)) - D_inv %*% A
##    }   
##    return(L)
##}


## faster
LAPnorm <- function(A, n = 2) {
  if (!(n %in% c(1, 2, 3))) n <- 2
  
  deg <- rowSums(A)
  
  if (n == 1) {
    # Unnormalized: L = D - A
    L <- -A
    diag(L) <- deg
    
  } else if (n == 2) {
    # Symmetric normalized: L = I - D^{-1/2} A D^{-1/2}
    # equivalent to: L[i,j] = -A[i,j] / sqrt(deg[i] * deg[j])
    d <- 1 / sqrt(deg)
    L <- (-A * outer(d, d))   # vectorized, no matrix mult
    diag(L) <- 1
    
  } else if (n == 3) {
    # Random-walk: L = I - D^{-1} A
    # equivalent to: L[i,j] = -A[i,j] / deg[i]
    L <- (-A / deg)             # broadcasts row-wise
    diag(L) <- 1
  }
  
  return(L)
}

