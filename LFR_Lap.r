rm(list=ls())

### Libraries ##################################3
library(optimx)
library(reticulate)
library(RSpectra)
library(igraph)
# Normlasied Laplacian
source('myLAPnorm.r') 
#        F <- S + T - P (see Paper)
source("EnergyFUNeigenLap.r") # energy_function which needs to be optimised

### read LFR network  #########
f1<-"graph_1668509426.txt" 
f2<-"membership_1668509426.txt"
temp<-read.table(f1)
temp<-as.matrix(temp)
temp<-temp+1 # because it strats at zero, edge lst needs values 1,2,3
g<-graph_from_edgelist(temp,directed = FALSE)
g<-simplify(g)
gORG<-g
adj<-as.matrix(as_adjacency_matrix(g,type="both"))
## MEMBERSHIP - ground truth communities
temp<-read.table(f2)
GroundTruth<-temp[,2] # Ground truth
cat('GroundTruth',length(unique(GroundTruth)),'\n')
N<-dim(adj)[1]
cat('Network size ',N,'\n')

### CALL Pyhonf form R !!!!!!!!!!!!!!
######### GRADCO PYTHON ##########################################################
# ── 0. Python environment ─────────────────────────────────────────────────────
use_python("/usr/bin/python3")   # adjust to your env 
# ── 1. Import modules ─────────────────────────────────────────────────────────
gradco  <- import("gradco")   
np     <- import("numpy")
#
adj_np <- np$array(adj, dtype = "int32")
# ── 3. Run GRADCO exactly as in your Python snippet ───────────────────────────
counter <- gradco$Counter(adj_np)  # grdco.Counter(adj)
counter$count()                   # counter.count()
# ── 4. Pull all 9 graphlet adjacency matrices into R ─────────────────────────
# index 0 = A0 (plain adjacency), 1–8 = graphlet matrices
graphlet_mats <- lapply(0:8, function(i) {
  temp  <- counter$get_graphlet_adjacency(i)   # counter.get_graphlet_adjacency(i)
  dense <- temp$todense()                       # temp.todense()  → numpy matrix
  py_to_r(dense)                                # numpy matrix → R matrix
})
names(graphlet_mats) <- paste0("A", 0:8)

#############################################################################
# 3D matrix
A<-array(0, c(N, N, 9));
A[,,1]<-graphlet_mats$A0 #ADJ 0 matrix
A[,,2]<-graphlet_mats$A1 #ADJ 1 matrix
A[,,3]<-graphlet_mats$A2 #ADJ 2 matrix
A[,,4]<-graphlet_mats$A3 #ADJ 3 matrix
A[,,5]<-graphlet_mats$A4 #ADJ 4 matrix
A[,,6]<-graphlet_mats$A5 #ADJ 5 matrix
A[,,7]<-graphlet_mats$A6 #ADJ 6 matrix
A[,,8]<-graphlet_mats$A7 #ADJ 7 matrix
A[,,9]<-graphlet_mats$A8 #ADJ 8 matrix

### ADJ 0 matrix
AOPT<-A[,,1]
## Normlasied Laplacian
L<-LAPnorm(AOPT,2) # call function LAPnorm / (2) Symmetric normalized
# 1. Perform eigen decomposition
EVmax<-20 # we do not need more than first 20 eigen values
ev <- eigs_sym(L, k = EVmax, which = "SM") # "SM"  The k eigenvalues with smallest magnitude.
# R returns eigenvalues in decreasing order. 
# We must sort them (and the eigenvectors) in ASCENDING order.
idx <- order(ev$values)
lambdas <- ev$values[idx][1:EVmax]
lambdas0<-lambdas
evecs <- ev$vectors[, idx]
# The gap is the difference between consecutive sorted eigenvalues
gaps <- diff(lambdas)
gaps[1]<-NA
k <- which.max(gaps)
cat("Number of clusters (k):", k, "\n")
# 3. Select the first k eigenvectors (smallest eigenvalues)
evecs <- evecs[, 1:k]
# normalise vectors
row_norms <- sqrt(rowSums(evecs^2))
evecs_norm <- evecs / (row_norms + 1e-10)
set.seed(123)
CL<-kmeans(evecs_norm, k,iter.max = 500, nstart = 10)
# NMI base line for A0 matrix only
cat('A0 NMI :',igraph::compare(CL$cluster,GroundTruth, method = "nmi"),'\n')

########## OPTIMISE section
cat('############ OPTIMISE ############# \n')
## normalise matrices
for (j in 2:9) {
     A[,,j]<-A[,,j]/max(A[,,j])
}
set.seed(123)
## Hooke and Jeeves Pattern Search Optimization
source("optHJN.r") # see EnergyFUNeigenLap.r for energy_function
# the output from optHJN.r is AOPT=w0*A0 + w1*A1 + ... w8*A8
# weights has ben optmised to get max of energy_function

## Normlasied Laplacian
L<-LAPnorm(AOPT,2) # AOPT is optimised linera combination Ai, sum, just one matrix
# 1. Perform eigendecomposition
EVmax<-20
ev <- eigs_sym(L, k = EVmax, which = "SM") # eigen(L)
# R returns eigenvalues in decreasing order. 
# We must sort them (and the eigenvectors) in ASCENDING order.
idx <- order(ev$values)
lambdas <- ev$values[idx][1:EVmax]
lambdasOPT<-lambdas
evecs <- ev$vectors[, idx]
# The gap is the difference between consecutive sorted eigenvalues
gaps <- diff(lambdas)

gaps[1]<-NA
k <- which.max(gaps)
cat("Optimal number of clusters (k):", k, "\n")
# 3. Select the first k eigenvectors (smallest eigenvalues)
evecs <- evecs[, 1:k]
# normalise vectors
row_norms <- sqrt(rowSums(evecs^2))
evecs_norm <- evecs / (row_norms + 1e-10)
set.seed(123)
CL<-kmeans(evecs_norm, k, iter.max = 500, nstart = 10)
cat('FINAL NMI :',igraph::compare(CL$cluster,GroundTruth, method = "nmi"),'\n')
cat('Lambda2(0, OPT, Ratio): ',lambdas0[2],lambdasOPT[2],lambdasOPT[2]/lambdas0[2],'\n')
finW<-c(1,result$par)
## plot
pdf(file='lfr.pdf', width=10, height=5)
    names(finW)<-c("G0","G1","G2","G3","G4","G5","G6","G7","G8")
    par(mfrow = c(1, 2))
    barplot(finW,main="LFR",ylim=c(0,2))
    Ymax<-max(lambdas0[2:EVmax],lambdasOPT[2:EVmax])+0.05
    Ymin<-min(lambdas0[2:EVmax],lambdasOPT[2:EVmax])-0.05
    plot(lambdas0[2:EVmax], typ='p', pch=16, col='black', ylim=c(Ymin,Ymax) )
    lines(lambdasOPT[2:EVmax],typ='p',pch=16,col='red')
dev.off()
