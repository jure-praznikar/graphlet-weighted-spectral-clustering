# Hooke–Jeeves pattern search over graphlet weights w1..w8
# (G0 weight stays fixed at 1 inside energy_function / AOPT below).
initial_pop <- rep(0.1, 8)

result <- hjn(
  par = initial_pop,
  fn = energy_function,
  control = list(maxfeval = 500, trace = 0, eps = 1.0e-3),
  lower = rep(0, 8),
  upper = rep(3, 8)
)

cat("Best fitness value: ", result$value, "\n")
cat("Weights: ", sprintf("%.3f", result$par), "\n")

# Build the optimised adjacency: AOPT = 1*A0 + w1*A1 + ... + w8*A8
AOPT <- A[, , 1]
for (i in 1:8) {
  AOPT <- AOPT + result$par[i] * A[, , i + 1]
}

# Same max-normalisation as inside energy_function
m <- max(AOPT)
if (m > 0) AOPT <- AOPT / m
