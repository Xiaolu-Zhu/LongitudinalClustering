rm(list=ls())
library(plyr);
library(fda);#bspline basis
library(Matrix);
library(matrixcalc);
library(igraph);
library(Rcpp);
library(RcppArmadillo);

set.seed(100)
sourceCpp("admmmcp_code.cpp");
source("clustering_functions.R")
#read data
data <- read.table("sampledata.txt")
#sample size
n = 60
#TT: max number of repeated measurements
TT = 10;

gamma1 = 0.005
gamma2 = 0.7

##########################################################################################
##############create design matrix(B spline basis) and second order penality for smoothing
nknots = 3;
order = 3;
p =  order + nknots;
#time points 
timerange = seq(0, 1, length.out = TT);

basis = dlply(data, .(id), function(xx) bsplineS(xx$time, knots_eq3(timerange, k = order, m = nknots), norder = order))
X = bdiag(basis)
X = as.matrix(X);

C <- matrix(0, nrow=nknots+order-2, ncol=nknots+order)
for (j in 1:(nknots+order-2)){
  d_j <- c(rep(0,j-1),1,-2,1,rep(0,(nknots+order)-3-(j-1)))
  e_j <- c(rep(0,j-1), 1 ,rep(0,(nknots+order)-3-(j-1)))
  C <- C + e_j%*%t(d_j)
}

D = t(C)%*%C;
diagD <- kronecker(diag(1, n), D);


####################### ADMM algorithm###########################
index = t(combn(n,2));
B_ini0 = update_B_ini(X, diagD, as.vector(data$y), n, gamma1, index, lambda0 = gamma1)

sol_final = prclust_admm(X, y=as.vector(data$y), diagD, B_ini0, index,
                            gamma1 = gamma1, gamma2 = gamma2, theta=1, tau = 2, n, p,  max_iter=1000,
                            eps_abs=1e-4, eps_rel=1e-4)

#result:
Ad_final <- create_adjacency(sol_final$V, n);
G_final <- graph.adjacency(Ad_final, mode = 'upper')
#clustering membership
cls_final <- components(G_final);
#number of clusters
k_final <- cls_final$no;


