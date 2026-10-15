# initial values of weights
initial_pop<-rep(0.1,8)

## Hooke and Jeeves Pattern Search Optimization
result <- hjn(par = initial_pop,
          fn = energy_function,
          control = list(maxfeval = 500, trace=0, eps=1.0E-3),
          lower = c(0, 0, 0, 0, 0, 0, 0, 0),
          upper = c(1, 1, 1, 1, 1, 1, 1, 1)*3
          )
#
cat('Best fitness value: ',result$value,'\n')
cat("Weights: ",sprintf("%.3f",result$par),'\n')

## Apply solution
AOPT <- A[,,1]
for (i in 1:8) {
      # linear combination, or sum, final result is one matrix > AOPT
      AOPT <- AOPT + result$par[i] * A[,,i+1]
}

# normalize (max=1)
m <- max(AOPT)
if (m > 0) AOPT <- AOPT / m


