## Core functions of pugwash to perform association

##################################################################
## SINGLE VARIANT TEST
## INPUT VARIABLES:
##   n        : total # of individuals
##   NS       : number of called samples
##   AC       : allele count
##   MAF      : minor allele frequency
##   vids     : indices from 1:n after AF/AC threshold
##   genos    : genotype matrix (after AF/AC threshold)
## EXPECTED OUTPUT : list(p, addcols, addnames) for each genos row
##   p        : p-value
##   addcols  : additional columns (with proper column names)
##################################################################

## single.b.firth :
## Use Firth bias-reduced logistic regression to perform association
## By: Clement Ma
## 
## Adapted from 'logistf' R package (v1.10)
## By:  Ploner M, Dunkler D, Southworth H, Heinze G
## TRAITS  : BINARY
## RETURNS : PVALUE, BETA, SEBETA, CHISQ
## MISSING VALUES : IGNORED

##################################################################
## HELPER FUNCTION:  fast.logistf.fit
## This function calculates the maximum penalized likelihood estimates
##################################################################
 
fast.logistf.fit <- function (x, y, weight = NULL, offset = NULL, firth = TRUE, col.fit = NULL, 
    init = NULL, control) {
    n <- nrow(x)
    k <- ncol(x)
    if (is.null(init)) 
        init = rep(0, k)
    if (is.null(col.fit)) 
        col.fit = 1:k
    if (is.null(offset)) 
        offset = rep(0, n)
    if (is.null(weight)) 
        weight = rep(1, n)
    if (col.fit[1] == 0) 
        maxit <- 0
    if (missing(control)) 
        control <- fast.logistf.control()
    maxit <- control$maxit
    maxstep <- control$maxstep
    maxhs <- control$maxhs
    lconv <- control$lconv
    gconv <- control$gconv
    xconv <- control$xconv
    beta <- init
    l.change <- 5
    iter <- 0
    pi <- as.vector(1/(1 + exp(-x %*% beta - offset)))
    loglik <- sum(weight[y == 1] * log(pi[y == 1])) + sum(weight[y == 
        0] * log(1 - pi[y == 0]))
    if (firth) {
        #XW2 <- crossprod(x, diag(weight * pi * (1 - pi))^0.5)
        XW2 <- t(x * (weight * pi * (1-pi))^0.5)
        Fisher <- crossprod(t(XW2))
        loglik <- loglik + 0.5 * determinant(Fisher)$modulus[1]
    }
    evals <- 1
    repeat {
        loglik.old <- loglik
        beta.old <- beta
        #XW2 <- crossprod(x, diag(weight * pi * (1 - pi))^0.5)
        XW2 <- t(x * (weight * pi * (1-pi))^0.5)
        #Fisher <- crossprod(t(XW2))
        #covs <- fast.invFisher(Fisher)
        #H <- crossprod(XW2, covs) %*% XW2
        myQR <- qr(t(XW2))
        Q <- qr.Q(myQR)
        h <- (Q*Q) %*% rep(1, ncol(Q))        
        if (firth) 
            #U.star <- crossprod(x, weight * (y - pi) + diag(H) * (0.5 - pi))
            U.star <- crossprod(x, weight * (y - pi) + h * (0.5 - pi))
        else U.star <- crossprod(x, weight * (y - pi))
        XX.covs <- matrix(0, k, k)
        if (col.fit[1] != 0) {
            #XX.XW2 <- crossprod(x[, col.fit, drop = FALSE], diag(weight * pi * (1 - pi))^0.5)
            XX.XW2 <- t(x[, col.fit, drop=FALSE] * (weight * pi * (1-pi))^0.5)
            XX.Fisher <- crossprod(t(XX.XW2))
            XX.covs[col.fit, col.fit] <- fast.invFisher(XX.Fisher)   ###### HERE IS THE PROBLEM!!!
        }
        if(all(is.na(XX.covs)) == T) {
            break
        }  
        delta <- as.vector(XX.covs %*% U.star)
        delta[is.na(delta)] <- 0
        mx <- max(abs(delta))/maxstep
        if (mx > 1) 
            delta <- delta/mx
        evals <- evals + 1
        if (maxit > 0) {
            iter <- iter + 1
            beta <- beta + delta
            for (halfs in 1:maxhs) {
                pi <- as.vector(1/(1 + exp(-x %*% beta - offset)))
                loglik <- sum(weight[y == 1] * log(pi[y == 1])) + 
                  sum(weight[y == 0] * log(1 - pi[y == 0]))
                if (firth) {
                  #XW2 <- crossprod(x, diag(weight * pi * (1 - pi))^0.5)
                  XW2 <- t(x * (weight * pi * (1-pi))^0.5)  
                  Fisher <- crossprod(t(XW2))
                  loglik <- loglik + 0.5 * fast.logDet(Fisher)
                }
                if(is.na(loglik)==T) {
                    break
                }
                evals <- evals + 1
                l.change <- loglik - loglik.old
                if (loglik > loglik.old) 
                  break
                beta <- beta - delta * 2^(-halfs)
            }
            if(is.na(loglik)==T) {
                break
            }
        }
        if(is.na(loglik)==T) {
            break
        }
        if (iter == maxit | ((max(abs(delta)) <= xconv) & (all(abs(U.star[col.fit]) < 
            gconv)) & (all(l.change < lconv)))) 
            break
    }
    # Error catching (if chol(x) not positive definite)
    if(all(is.na(XX.covs))==T | is.na(loglik)==T) {
        var <- XX.covs
        list(beta = NA, var = var, pi = NA, hat.diag = NA, 
        loglik = NA, iter = NA, evals = NA, conv = c(NA, 
            NA, NA))
    } else {
        var <- XX.covs
        list(beta = beta, var = var, pi = pi, hat.diag = h, 
            loglik = loglik, iter = iter, evals = evals, conv = c(l.change, 
                max(abs(U.star)), max(abs(delta))))
    }
}


##################################################################
## HELPER FUNCTION:  fast.logistf.control
## This function provides the convergence control parameters
##################################################################
fast.logistf.control <- function (maxit = 50, maxhs = 15, maxstep = 15, lconv = 1e-05, 
    gconv = 1e-05, xconv = 1e-05) 
{
    list(maxit = maxit, maxhs = maxhs, maxstep = maxstep, lconv = lconv, 
        gconv = gconv, xconv = xconv)
}

##################################################################
## HELPER FUNCTION:  fast.logDet
##################################################################
fast.logDet <- function (x) {
    my.chol <- tryCatch(chol(x),error=function(e) {NA})
    if (all(is.na(my.chol))==T) {
        return(NA)
    } else {
        return (2 * sum(log(diag(my.chol))))
    }
}   

##################################################################
## HELPER FUNCTION:  fast.invFisher
##################################################################
fast.invFisher <- function(x) {
  my.chol <- tryCatch(chol(x),error=function(e) {NA})
    if (all(is.na(my.chol))==T) {
        return(NA)
    } else {
        return (chol2inv(my.chol))
    }
  #ifelse(is.na(my.chol), NA, chol2inv(my.chol))
}

