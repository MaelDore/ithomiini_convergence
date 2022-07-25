
###################################
#      Author: Maël Doré          #
#  Contact: mael.dore@gmail.com   #
###################################

# Script adapted from Revell, 2009: DOI = 10.1111/j.1558-5646.2009.00804.x


Revell_phyl_pca <- function(C,X,mode) {
  # find out how many columns and taxa we have
  m <- ncol(X); n <- nrow(X);
  # compute inverse of C with C being the phylogenetic variance-covariance matrix extracted from a tree (vcv)
  invC <- solve(C);
  # compute vector of ancestral states
  one <- matrix(1,n,1);
  a <- t(t(one)%*%invC%*%X)*sum(sum(invC))^-1;
  # compute evolutionary VCV matrix
  V <- t(X-one%*%t(a))%*%invC%*%(X-one%*%t(a))*(n-1)^-1;
  evolVCV <- V # Store evolutionary VCV matrix
  # if correlation matrix
  if (mode=="corr") {
    # standardize X
    X = X/(one%*%t(sqrt(diag(V))));
    # change V to correlation matrix
    V = V/(sqrt(diag(V))%*%t(sqrt(diag(V))));
    evol_corr <- V # Store evolutionary correlation matrix
    # recalculate a
    a <- t(t(one)%*%invC%*%X)*sum(sum(invC))^-1;
  }
  # eigenanalyze
  es = eigen(V);
  result <- NULL; 
  result$Eval <- diag(es$values);
  result$Evec <- es$vectors;
  # objects used in the computation process
  result$ancestral <- a  # Ancestral states of character/traits
  result$evolVCV <- evolVCV  # Evolutionary variance/covariance matrix
  if (mode=="corr") {result$evol_corr <- evol_corr} # Evolutionary correlation matrix
  result$X <- X  # Trait data used in the computation standardized (or not)
  # compute scores in the species space
  result$S <- (X-one%*%t(a))%*%result$Evec;
  # compute cross covariance matrix
  # and loadings
  Ccv <- t(X-one%*%t(a))%*%invC%*%result$S/(n-1);
  result$L <- matrix(,m,m);
  for(i in 1:m)
    for(j in 1:m)
      result$L[i,j] <- Ccv[i,j]/sqrt(V[i,i]*result$Eval[j,j]);
  Revell_phyl_pca <- result;
  # done
}