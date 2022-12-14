# Modified version of ecodist::MRM function that can provide the statistics obtained from permutated data

# ?ecodist::MRM

# Credits for original version: Sarah Goslee
# Credits for modified version: Maël Doré

MRM_verbose <- function (formula = formula(data), data, 
                         nperm = 1000, method = "linear", 
                         mrank = FALSE, perm_output = FALSE) 
{
  m <- match.call(expand.dots = FALSE)
  m2 <- match(c("formula", "data"), names(m), nomatch = 0)
  m <- m[c(1, m2)]
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  m <- as.matrix(m)
  n <- (1 + sqrt(1 + 8 * nrow(m)))/2
  if (abs(n - round(n)) > 1e-07)
    stop("Matrix not square.\n")
  n <- round(n)
  if (ncol(m) < 2) 
    stop("Not enough data. \n")
  
  if (method == "linear")
  {
    if (mrank) 
    {
      m <- apply(m, 2, rank)
    }
    for (thiscol in seq_len(ncol(m)))
    {
      tempmat <- full(m[, thiscol])
      m[, thiscol] <- tempmat[col(tempmat) > row(tempmat)]
    }
    X <- m[, 2:ncol(m), drop = FALSE]
    X <- cbind(rep(1, nrow(X)), X)
    Y <- m[, 1, drop = FALSE]
    nd <- nrow(X)
    XX <- crossprod(X)
    XX <- solve(XX)
    XY <- crossprod(X, Y)
    YY <- crossprod(Y)
    b <- XX %*% XY
    rownames(b) <- c("Int", colnames(X)[2:ncol(X)])
    bXY <- crossprod(b, XY)
    SSE <- YY - bXY
    SSTO <- YY - sum(Y)^2/nd
    SSR = SSTO - SSE
    R2 <- 1 - SSE/SSTO
    R2 <- as.vector(R2)
    p <- ncol(X)
    F <- (SSR/(p - 1))/(SSE/(nd - p))
    R2.pval <- NA
    b.pval <- rep(NA, ncol(X))
    F.pval <- NA
    
    # Build result output if no permutations
    results <- list(coef = cbind(b, pval = b.pval),
                    r.squared = c(R2 = R2, pval = R2.pval),
                    F.test = c(F = F, F.pval = F.pval))
    
    if (nperm > 0)
    {
      var_names <- row.names(b)
      R2.all <- numeric(nperm)
      b.all <- numeric(nperm * p)
      F.all <- numeric(nperm)
      cresults <- .C("mrmperm", as.double(as.vector(X)), 
                     as.double(as.vector(Y)), as.integer(p), as.integer(nd), 
                     as.integer(n), as.integer(nperm), R2.all = as.double(R2.all), 
                     b.all = as.double(b.all), F.all = as.double(F.all), 
                     as.double(numeric(n * n)), as.integer(numeric(n)), 
                     as.double(as.vector(XX)), as.double(numeric(p)), 
                     as.double(0), as.double(numeric(p)), PACKAGE = "ecodist")
      
      R2.all <- cresults$R2.all
      R2 <- R2.all[1]
      R2.Q95 <- unname(quantile(x = R2.all, probs = 0.95))
      R2.mean <- mean(R2.all)
      R2.pval <- length(R2.all[R2.all >= R2.all[1]])/nperm
      
      F.all <- cresults$F.all
      F <- F.all[1]
      F.Q95 <- unname(quantile(x = F.all, probs = 0.95))
      F.mean <- mean(F.all)
      F.pval <- length(F.all[F.all >= F.all[1]])/nperm
      
      b.all <- matrix(cresults$b.all, nrow = nperm, ncol = p, 
                      byrow = TRUE)
      b <- b.all[1, ]
      names(b) <- var_names
      b.Q95 <- unlist(apply(X = b.all, MARGIN = 2, FUN = quantile, probs = 0.95))
      b.mean <- unlist(apply(X = b.all, MARGIN = 2, FUN = mean))
      b.pval <- apply(b.all, 2, function(x) length(x[abs(x) >= 
                                                       abs(x[1])])/nperm)
      
      # Build result output with permutations
      results <- list(coef = cbind(beta.coef = b, mean.NULL = b.mean, Q95 = b.Q95, pval = b.pval),
                      r.squared = c(R2 = R2, mean.NULL = R2.mean, Q95 = R2.Q95, R2.pval = R2.pval),
                      F.test = c(F = F, mean.NULL = F.mean, Q95 = F.Q95, F.pval = F.pval))
      if (perm_output)
      {
        results <- append(results, values = list(permutations = list(b.perm = b.all, R2.perm = R2.all, F.perm = F.all)))
      }
    }
    
  } else {
    
    if (method == "logistic")
    {
      X <- m[, 2:ncol(m), drop = FALSE]
      Y <- m[, 1, drop = FALSE]
      colnames(Y) <- "Y"
      newdata <- data.frame(Y = Y, X)
      fit1 <- glm(Y ~ ., data = newdata, family = binomial(link = "logit"))
      b <- coefficients(fit1)
      dev <- summary(fit1)$deviance
      dev.df <- summary(fit1)$df.residual
      b.pval <- NA
      dev.pval <- NA
      
      # Build result output if no permutations
      results <- list(coef = cbind(b, pval = b.pval), 
                      dev = c(resid.dev = dev, resid.df = dev.df, dev.pval = dev.pval))
      
      if (nperm > 0) {
        b.all <- matrix(NA, nrow = nperm, ncol = length(b))
        b.all[1, ] <- b
        dev.all <- rep(NA, nperm)
        dev.all[1] <- dev
        for (i in 2:nperm) {
          newSample <- sample(n)
          newY <- full(Y)
          newY <- newY[newSample, newSample]
          newY <- lower(newY)
          newdata <- data.frame(Y = newY, X = X)
          newfit <- glm(Y ~ ., data = newdata, family = binomial(link = "logit"))
          b.all[i, ] <- coefficients(newfit)
          dev.all[i] <- summary(newfit)$deviance
        }
        
        b.Q95 <- unlist(apply(X = b.all, MARGIN = 2, FUN = quantile, probs = 0.95))
        b.mean <- unlist(apply(X = b.all, MARGIN = 2, FUN = mean))
        b.pval <- apply(b.all, 2, function(x) length(x[abs(x) >= abs(x[1])])/nperm)
        
        dev.Q95 <- unname(quantile(x = dev.all, probs = 0.95))
        dev.mean <- mean(dev.all)
        dev.pval <- length(dev.all[dev.all >= dev.all[1]])/nperm
        
        # Build result output with permutations
        results <- list(coef = cbind(beta.coef = b, mean.NULL = b.mean, Q95 = b.Q95, pval = b.pval), 
                        dev = c(resid.dev = dev, resid.df = dev.df, mean.NULL = dev.mean, Q95 = dev.Q95, dev.pval = dev.pval))
        
        if (perm_output)
        {
          results <- append(results, values = list(permutations = list(b.perm = b.all, dev.perm = dev.all)))
        }
      }
      
    } else {
      stop("method must be 'linear' or 'logistic'\n")
    }
  }
  results
}
