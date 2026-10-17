rm(list = ls())

# =============================================================================
# Libraries and helpers
# =============================================================================
library(optimx)
library(reticulate)
library(RSpectra)
library(igraph)

source("myLAPnorm.r")           # LAPnorm(): normalised Laplacian
source("EnergyFUNeigenLap.r")  # energy_function F = S + T - P (to maximise)

# =============================================================================
# 1. Read LFR network and ground-truth communities
# =============================================================================
f1 <- "graph_1668509426.txt"
f2 <- "membership_1668509426.txt"

# Edge list is 0-indexed; igraph expects 1-based vertex ids
temp <- as.matrix(read.table(f1))
temp <- temp + 1
g <- graph_from_edgelist(temp, directed = FALSE)
g <- simplify(g)
gORG <- g
adj <- as.matrix(as_adjacency_matrix(g, type = "both"))

GroundTruth <- read.table(f2)[, 2]
cat("GroundTruth communities:", length(unique(GroundTruth)), "\n")
N <- dim(adj)[1]
cat("Network size:", N, "\n")

# =============================================================================
# 2. Graphlet adjacency matrices via GRADCO (Python)
# =============================================================================
use_python("/usr/bin/python3")  # adjust to your environment if needed

gradco <- import("gradco")
np     <- import("numpy")

adj_np <- np$array(adj, dtype = "int32")
counter <- gradco$Counter(adj_np)
counter$count()

# Index 0 = plain adjacency A0; indices 1–8 = graphlet adjacency matrices
graphlet_mats <- lapply(0:8, function(i) {
  temp  <- counter$get_graphlet_adjacency(i)
  dense <- temp$todense()
  py_to_r(dense)
})
names(graphlet_mats) <- paste0("A", 0:8)

# Stack into a 3D array: A[,,1] = A0, ..., A[,,9] = A8
A <- array(0, c(N, N, 9))
for (i in 1:9) {
  A[, , i] <- graphlet_mats[[i]]
}

# =============================================================================
# 3. Baseline spectral clustering on A0 only
# =============================================================================
AOPT <- A[, , 1]
L <- LAPnorm(AOPT, 2)  # symmetric normalised Laplacian

# Smallest EVmax eigenvalues (RSpectra returns them unordered)
EVmax <- 20
ev  <- eigs_sym(L, k = EVmax, which = "SM")
idx <- order(ev$values)              # ascending: lambda_1 <= ... <= lambda_EVmax
lambdas  <- ev$values[idx][1:EVmax]
lambdas0 <- lambdas                  # keep baseline spectrum for energy_function / plot
evecs    <- ev$vectors[, idx]

# Eigengap heuristic: ignore gap(lambda_2 - lambda_1) so k >= 2
gaps <- diff(lambdas)
gaps[1] <- NA
k <- which.max(gaps)
cat("Number of clusters (k):", k, "\n")

# Ng–Jordan–Weiss: first k eigenvectors, row-normalise, then k-means
evecs <- evecs[, 1:k, drop = FALSE]
row_norms  <- sqrt(rowSums(evecs^2))
evecs_norm <- evecs / (row_norms + 1e-10)
set.seed(123)
CL <- kmeans(evecs_norm, k, iter.max = 500, nstart = 10)
cat("A0 NMI :", igraph::compare(CL$cluster, GroundTruth, method = "nmi"), "\n")

# =============================================================================
# 4. Optimise graphlet weights (Hooke–Jeeves)
# =============================================================================
cat("############ OPTIMISE ############# \n")

# Scale G1–G8 to max = 1 so weights are comparable (G0 stays 0/1)
for (j in 2:9) {
  A[, , j] <- A[, , j] / max(A[, , j])
}

set.seed(123)
# optHJN.r maximises energy_function over w1..w8 (G0 weight fixed at 1)
# and writes AOPT = G0 + w1*G1 + ... + w8*G8 (then normalised to max = 1)
source("optHJN.r")

# =============================================================================
# 5. Spectral clustering on the optimised combination AOPT
# =============================================================================
L <- LAPnorm(AOPT, 2)
ev  <- eigs_sym(L, k = EVmax, which = "SM")
idx <- order(ev$values)
lambdas    <- ev$values[idx][1:EVmax]
lambdasOPT <- lambdas
evecs      <- ev$vectors[, idx]

gaps <- diff(lambdas)
gaps[1] <- NA
k <- which.max(gaps)
cat("Optimal number of clusters (k):", k, "\n")

evecs <- evecs[, 1:k, drop = FALSE]
row_norms  <- sqrt(rowSums(evecs^2))
evecs_norm <- evecs / (row_norms + 1e-10)
set.seed(123)
CL <- kmeans(evecs_norm, k, iter.max = 500, nstart = 10)

cat("FINAL NMI :", igraph::compare(CL$cluster, GroundTruth, method = "nmi"), "\n")
cat("Lambda2(0, OPT, Ratio): ",
    lambdas0[2], lambdasOPT[2], lambdasOPT[2] / lambdas0[2], "\n")

# Full weight vector: G0 fixed at 1, then optimised w1..w8
weights <- c(1, result$par)
names(weights) <- paste0("G", 0:8)

# =============================================================================
# 6. Plot weights and spectra
# =============================================================================
plot_file <- "lfr.pdf"
pdf(file = plot_file, width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(6, 4, 3, 1))

# Left: barplot of the nine weights (G0 is fixed at 1)
barplot(
  weights,
  main = "Graphlet weights",
  ylab = "weight",
  ylim = c(0, max(weights) * 1.15 + 0.05),
  las = 2
)

# Right: eigenvalues from lambda_2 onward (drop lambda_1 ~ 0)
y0 <- lambdas0[-1]
y1 <- lambdasOPT[-1]
plot(
  y0,
  type = "p",
  pch = 16,
  col = "black",
  ylim = c(min(y0, y1) - 0.05, max(y0, y1) + 0.05),
  xlab = expression(paste("index, starting at ", lambda[2])),
  ylab = "eigenvalue",
  main = "Laplacian spectrum"
)
points(y1, pch = 16, col = "red")
legend(
  "bottomright",
  legend = c("G0 only", "optimised"),
  col = c("black", "red"),
  pch = 16,
  bty = "n"
)

dev.off()
cat("Wrote", plot_file, "\n")
