#======================================================================================
# Local likelihood estimation for covariance functions with spatially-varying
# parameters: the convoSPAT() package for R
# Mark D. Risser / The Ohio State University / 2014-2015
#======================================================================================

#======================================================================================
# Local/Global Parameter Estimation
#======================================================================================

#======================================================================================
# Global parameter estimation
#======================================================================================
# For a given correlation matrix, the following functions calculate the global
# parameters based on what is NOT specified to be spatially-varying.
#======================================================================================

#===================================
# First, models without kappa
# Estimates: nugget, variance

#ROxygen comments ----
#' Constructor functions for global parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of
#' global variance parameters tausq, sigmasq with a fixed correlation
#' matrix (smoothness is fixed).
#'
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#' @param Corr The correlation matrix.
#' @param nugg2.var Fixed values for the covariance of the second nugget term.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_global_loglik1( data, Xmat, Corr, nugg2.var )
#' }
#'
#' @export

make_global_loglik1 <- function( data, Xmat, Corr, nugg2.var ){

  fixed <- c(FALSE, FALSE)
  params <- fixed
  function(p){
    params[!fixed] <- p

    # Parameters
    tausq <- params[1]
    sigmasq <- params[2]

    N <- dim(data)[1]
    p <- dim(data)[2]

    Edecomp <- eigen(sigmasq*Corr + diag(rep(tausq,N)) + nugg2.var)
    Dmat.temp <- diag(Edecomp$values)
    Vmat.temp <- Edecomp$vectors

    Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

    XCX <- t(Xmat) %*% Cov.inv %*% Xmat
    XCX.inv <- chol2inv( chol(XCX) )

    Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

    # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
    loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(data)%*%Ptemp%*%data)) )

    if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
    return(loglikelihood)

  }
}

# Estimates: variance

#ROxygen comments ----
#' Constructor functions for global parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of
#' global variance parameter sigmasq with a fixed correlation
#' matrix (smoothness is fixed). The nugget variance is taken
#' to be spatially-varing.
#'
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#' @param Corr The correlation matrix.
#' @param obs.nuggets A vector containing the spatially-varying nuggets
#' corresponding to each data location.
#' @param nugg2.var Fixed values for the covariance of the second nugget term.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_global_loglik2( data, Xmat, Corr, obs.nuggets, nugg2.var )
#' }
#'
#' @export

make_global_loglik2 <- function( data, Xmat, Corr, obs.nuggets, nugg2.var ){

  fixed <- c(FALSE)
  params <- fixed
  function(p){
    params[!fixed] <- p

    # Parameters
    sigmasq <- params[1]

    N <- dim(data)[1]
    p <- dim(data)[2]

    Edecomp <- eigen(sigmasq*Corr + diag(obs.nuggets) + nugg2.var)
    Dmat.temp <- diag(Edecomp$values)
    Vmat.temp <- Edecomp$vectors

    Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

    XCX <- t(Xmat) %*% Cov.inv %*% Xmat
    XCX.inv <- chol2inv( chol(XCX) )

    Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

    # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
    loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(data)%*%Ptemp%*%data)) )

    if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
    return(loglikelihood)

  }
}


# Estimates: nugget

#ROxygen comments ----
#' Constructor functions for global parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of
#' global variance parameter tausq with a fixed correlation
#' matrix (smoothness is fixed). The process variance is taken
#' to be spatially-varing.
#'
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#' @param Corr The correlation matrix matrix.
#' @param obs.variance A vector containing the spatially-varying variance
#' corresponding to each data location.
#' @param nugg2.var Fixed values for the covariance of the second nugget term.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_global_loglik3( data, Xmat, Corr, obs.variance, nugg2.var )
#' }
#'
#' @export

make_global_loglik3 <- function( data, Xmat, Corr, obs.variance, nugg2.var ){

  fixed <- c(FALSE)
  params <- fixed
  function(p){
    params[!fixed] <- p

    tausq <- params[1]

    N <- dim(data)[1]
    p <- dim(data)[2]

    Edecomp <- eigen(diag(sqrt(obs.variance)) %*% Corr %*% diag(sqrt(obs.variance))
                     + diag(rep(tausq,N)) + nugg2.var )
    Dmat.temp <- diag(Edecomp$values)
    Vmat.temp <- Edecomp$vectors

    Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

    XCX <- t(Xmat) %*% Cov.inv %*% Xmat
    XCX.inv <- chol2inv( chol(XCX) )

    Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

    # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
    loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(data)%*%Ptemp%*%data)) )

    if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
    return(loglikelihood)

  }
}

#===================================

#===================================
# Next, models with kappa
# Estimates: nugget, variance, kappa

#ROxygen comments ----
#' Constructor functions for global parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of
#' global variance parameters tausq, sigmasq, and nu.
#'
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#' @param cov.model String; the covariance model.
#' @param Scalemat Matrix; contains the scaling quantities from the
#' covariance function.
#' @param Distmat Matrix; contains the scaled distances.
#' @param nugg2.var Fixed values for the covariance of the second nugget term.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_global_loglik1_kappa( data, Xmat, cov.model, Scalemat, Distmat, nugg2.var )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial

make_global_loglik1_kappa <- function( data, Xmat, cov.model, Scalemat, Distmat, nugg2.var ){

  fixed <- c(FALSE, FALSE, FALSE)
  params <- fixed
  function(p){
    params[!fixed] <- p

    # Parameters
    tausq <- params[1]
    sigmasq <- params[2]
    kapp <- params[3]

    Unscl.corr <- cov.spatial( Distmat, cov.model = cov.model,
                               cov.pars = c(1,1), kappa = kapp )
    Corr <- Scalemat*Unscl.corr

    N <- dim(data)[1]
    p <- dim(data)[2]

    Edecomp <- eigen(sigmasq*Corr + diag(rep(tausq,N)) + nugg2.var)
    Dmat.temp <- diag(Edecomp$values)
    Vmat.temp <- Edecomp$vectors

    Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

    XCX <- t(Xmat) %*% Cov.inv %*% Xmat
    XCX.inv <- chol2inv( chol(XCX) )

    Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

    # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
    loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(data)%*%Ptemp%*%data)) )

    if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
    return(loglikelihood)

  }
}


# Estimates: variance, kappa
#ROxygen comments ----
#' Constructor functions for global parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of
#' global variance parameters sigmasq and nu. The nugget variance is
#' taken to be spatially-varying.
#'
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#' @param cov.model String; the covariance model.
#' @param Scalemat Matrix; contains the scaling quantities from the
#' covariance function.
#' @param Distmat Matrix; contains the scaled distances.
#' @param obs.nuggets A vector containing the spatially-varying nuggets
#' corresponding to each data location.
#' @param nugg2.var Fixed values for the covariance of the second nugget term.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_global_loglik2_kappa( data, Xmat, cov.model, Scalemat, Distmat, obs.nuggets, nugg2.var )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial

make_global_loglik2_kappa <- function( data, Xmat, cov.model, Scalemat, Distmat, obs.nuggets, nugg2.var ){

  fixed <- c(FALSE, FALSE)
  params <- fixed
  function(p){
    params[!fixed] <- p

    # Parameters
    sigmasq <- params[1]
    kapp <- params[2]

    Unscl.corr <- cov.spatial( Distmat, cov.model = cov.model,
                               cov.pars = c(1,1), kappa = kapp )
    Corr <- Scalemat*Unscl.corr

    N <- dim(data)[1]
    p <- dim(data)[2]

    Edecomp <- eigen(sigmasq * Corr + diag(obs.nuggets) + nugg2.var)
    Dmat.temp <- diag(Edecomp$values)
    Vmat.temp <- Edecomp$vectors

    Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

    XCX <- t(Xmat) %*% Cov.inv %*% Xmat
    XCX.inv <- chol2inv( chol(XCX) )

    Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

    # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
    loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(data)%*%Ptemp%*%data)) )

    if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
    return(loglikelihood)

  }
}

# Estimates: nugget, kappa
#ROxygen comments ----
#' Constructor functions for global parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of
#' global variance parameters tausq and nu. The process variance is
#' taken to be spatially-varying.
#'
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#' @param cov.model String; the covariance model.
#' @param Scalemat Matrix; contains the scaling quantities from the
#' covariance function.
#' @param Distmat Matrix; contains the scaled distances.
#' @param obs.variance A vector containing the spatially-varying variance
#' corresponding to each data location.
#' @param nugg2.var Fixed values for the covariance of the second nugget term.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_global_loglik3_kappa( data, Xmat, cov.model, Scalemat, Distmat, obs.variance, nugg2.var )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial

make_global_loglik3_kappa <- function( data, Xmat, cov.model, Scalemat, Distmat, obs.variance, nugg2.var ){

  fixed <- c(FALSE, FALSE)
  params <- fixed
  function(p){
    params[!fixed] <- p

    # Parameters
    tausq <- params[1]
    kapp <- params[2]

    Unscl.corr <- cov.spatial( Distmat, cov.model = cov.model,
                               cov.pars = c(1,1), kappa = kapp )
    Corr <- Scalemat*Unscl.corr

    N <- dim(data)[1]
    p <- dim(data)[2]

    Edecomp <- eigen(diag(sqrt(obs.variance)) %*% Corr %*% diag(sqrt(obs.variance))
                     + diag(rep(tausq,N)) + nugg2.var )
    Dmat.temp <- diag(Edecomp$values)
    Vmat.temp <- Edecomp$vectors

    Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

    XCX <- t(Xmat) %*% Cov.inv %*% Xmat
    XCX.inv <- chol2inv( chol(XCX) )

    Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

    # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
    loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(data)%*%Ptemp%*%data)) )

    if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
    return(loglikelihood)

  }
}


# Estimates: kappa
#ROxygen comments ----
#' Constructor functions for global parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of
#' global variance parameters nu. The process variance
#' and nugget variance are taken to be spatially-varying.
#'
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#' @param cov.model String; the covariance model.
#' @param Scalemat Matrix; contains the scaling quantities from the
#' covariance function.
#' @param Distmat Matrix; contains the scaled distances.
#' @param obs.variance A vector containing the spatially-varying variance
#' corresponding to each data location.
#' @param obs.nuggets A vector containing the spatially-varying nuggets
#' corresponding to each data location.
#' @param nugg2.var Fixed values for the covariance of the second nugget term.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_global_loglik4_kappa( data, Xmat, cov.model, Scalemat, Distmat,
#' obs.variance, obs.nuggets, nugg2.var )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial

make_global_loglik4_kappa <- function( data, Xmat, cov.model, Scalemat, Distmat, obs.variance,
                                       obs.nuggets, nugg2.var ){

  fixed <- c(FALSE)
  params <- fixed
  function(p){
    params[!fixed] <- p

    # Parameters
    kapp <- params[1]

    Unscl.corr <- cov.spatial( Distmat, cov.model = cov.model,
                               cov.pars = c(1,1), kappa = kapp )
    Corr <- Scalemat*Unscl.corr

    N <- dim(data)[1]
    p <- dim(data)[2]

    Edecomp <- eigen(diag( sqrt(obs.variance) ) %*% Corr %*% diag( sqrt(obs.variance) ) + diag(obs.nuggets) + nugg2.var)
    Dmat.temp <- diag(Edecomp$values)
    Vmat.temp <- Edecomp$vectors

    Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

    XCX <- t(Xmat) %*% Cov.inv %*% Xmat
    XCX.inv <- chol2inv( chol(XCX) )

    Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

    # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
    loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(data)%*%Ptemp%*%data)) )

    if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
    return(loglikelihood)

  }
}

#======================================================================================

#======================================================================================
# Function to estimate the local models
#======================================================================================
# Using a subset of the data, calculate MLEs of covariance and mean
# parameters. Options for
# (1) method: ml vs. reml
# (2) smoothness parameter: no kappa vs. estimate kappa vs. fixed kappa
# (3) locally isotropic vs. locally anisotropic
# (4) fixed tausq vs. estimate tausq
#======================================================================================

#ROxygen comments ----
#' Constructor functions for local parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of covariance (and possibly mean) parameters.
#' The function includes options for
#' (1) maximum likelihood (\code{"ml"}) vs. restricted maximum likelihood
#'     (\code{"reml"}),
#' (2) smoothness (\code{kappa}): models without smoothness vs. estimating the
#'     smoothness vs. using fixed smoothness,
#' (3) locally isotropic vs. locally anisotropic, and
#' (4) fixed nugget variance (\code{tausq}): fixed vs. estimated.
#'
#' @param locations A matrix of locations.
#' @param cov.model String; the covariance model.
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#' @param nugg2.var Fixed values for the variance/covariance of the second nugget term; defaults
#' to a matrix of zeros.
#' @param tausq Scalar; fixed value for the nugget variance (when
#' \code{fix.tausq = TRUE}).
#' @param kappa Scalar; fixed value for the smoothness (when \code{fix.kappa = TRUE}).
#' @param fixed Logical vector of \code{FALSE} values; length corresponds to the number
#' of parameters to be estimated.
#' @param method Indicates the estimation method, either maximum likelihood (\code{"ml"})
#' or restricted maximum likelihood (\code{"reml"}).
#' @param local.aniso Logical; indicates if the local covariance should be
#' anisotropic (\code{TRUE}) or isotropic (\code{FALSE}). Defaults to \code{TRUE}.
#' @param fix.kappa Logical; indicates if the kappa parameter should be
#' fixed (\code{TRUE}) or estimated (\code{FALSE}). Defaults to \code{FALSE}
#' (only valid for \code{cov.model = "matern"} and \code{cov.model = "cauchy"}).
#' @param fix.tausq Logical; indicates whether the default nugget term
#' (tau^2) should be fixed (\code{TRUE}) or estimated (\code{FALSE}). Defaults to
#' \code{FALSE}.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_local_lik( locations, cov.model, data, Xmat )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial
#' @importFrom StatMatch mahalanobis.dist

make_local_lik <- function( locations, cov.model, data, Xmat,
                            nugg2.var = matrix(0, nrow(locations), nrow(locations)),
                            tausq = 0,
                            kappa = 0.5, fixed = rep(FALSE, 6),
                            method = "reml", local.aniso = TRUE,
                            fix.tausq = FALSE, fix.kappa = FALSE ){

  params <- fixed

  function(p) {
    params[!fixed] <- p

    if( local.aniso ){ # Locally anisotropic
      lam1 <- params[1]
      lam2 <- params[2]
      eta <- params[3]

      if( cov.model == "matern" || cov.model == "cauchy" ){ # For covariance models with kappa
        if( !fix.kappa ){ # Estimate kappa
          if( !fix.tausq ){ # Estimate tausq
            tausq <- params[4]
            sigmasq <- params[5]
            kappa <- params[6]
            if( method == "ml" ) beta <- params[-c(1:6)]
          }
          if( fix.tausq ){ # Fix tausq
            tausq <- tausq
            sigmasq <- params[4]
            kappa <- params[5]
            if( method == "ml" ) beta <- params[-c(1:5)]
          }
        }
        if( fix.kappa ){ # Fix kappa
          if( !fix.tausq ){ # Estimate tausq
            tausq <- params[4]
            sigmasq <- params[5]
            kappa <- kappa
            if( method == "ml" ) beta <- params[-c(1:5)]
          }
          if( fix.tausq ){ # Fix tausq
            tausq <- tausq
            sigmasq <- params[4]
            kappa <- kappa
            if( method == "ml" ) beta <- params[-c(1:4)]
          }
        }
      }
      if( cov.model != "matern" & cov.model != "cauchy" ){ # For covariance models without kappa
        if( !fix.tausq ){ # Estimate tausq
          tausq <- params[4]
          sigmasq <- params[5]
          if( method == "ml" ) beta <- params[-c(1:5)]
        }
        if( fix.tausq ){ # Fix tausq
          tausq <- tausq
          sigmasq <- params[4]
          if( method == "ml" ) beta <- params[-c(1:4)]
        }
      }

    }
    if( !local.aniso ){ # Locally isotropic
      lam1 <- params[1]
      lam2 <- params[1]
      eta <- 0

      if( cov.model == "matern" || cov.model == "cauchy" ){ # For covariance models with kappa
        if( !fix.kappa ){ # Estimate kappa
          if( !fix.tausq ){ # Estimate tausq
            tausq <- params[2]
            sigmasq <- params[3]
            kappa <- params[4]
            if( method == "ml" ) beta <- params[-c(1:4)]
          }
          if( fix.tausq ){ # Fix tausq
            tausq <- tausq
            sigmasq <- params[2]
            kappa <- params[3]
            if( method == "ml" ) beta <- params[-c(1:3)]
          }
        }
        if( fix.kappa ){ # Fix kappa
          if( !fix.tausq ){ # Estimate tausq
            tausq <- params[2]
            sigmasq <- params[3]
            kappa <- kappa
            if( method == "ml" ) beta <- params[-c(1:3)]
          }
          if( fix.tausq ){ # Fix tausq
            tausq <- tausq
            sigmasq <- params[2]
            kappa <- kappa
            if( method == "ml" ) beta <- params[-c(1:2)]
          }
        }
      }
      if( cov.model != "matern" & cov.model != "cauchy" ){ # For covariance models without kappa
        if( !fix.tausq ){ # Estimate tausq
          tausq <- params[2]
          sigmasq <- params[3]
          if( method == "ml" ) beta <- params[-c(1:3)]
        }
        if( fix.tausq ){ # Fix tausq
          tausq <- tausq
          sigmasq <- params[2]
          if( method == "ml" ) beta <- params[-c(1:2)]
        }
      }
    }

    # Setup and covariance calculation
    N <- dim(locations)[1]
    m <- dim(data)[2]

    Pmat <- matrix(c(cos(eta), -sin(eta), sin(eta), cos(eta)), nrow = 2, byrow = T)
    Dmat <- diag(c(lam1, lam2))
    Sigma <- Pmat %*% Dmat %*% t(Pmat)
    distances <- mahalanobis.dist(data.x = locations, vc = Sigma)
    NS.cov <- nugg2.var + tausq * diag(N) + sigmasq * cov.spatial(distances, cov.model = cov.model,
                                                                        cov.pars = c(1, 1), kappa = kappa)
    cov.chol <- chol(NS.cov)

    # Likelihood calculation
    if( method == "ml" ){
      tmp1 <- backsolve(cov.chol, data - Xmat %*% beta, transpose = TRUE)
      ResCinvRes <- t(tmp1) %*% tmp1
      loglikelihood <- m * sum(log(diag(cov.chol))) + 0.5 * sum(diag(ResCinvRes))
    }
    if( method == "reml" ){
      tmp1 <- backsolve(cov.chol, Xmat, transpose = TRUE)
      XCinvX <- t(tmp1) %*% tmp1
      xcx.chol <- chol(XCinvX)
      tmp2 <- backsolve(cov.chol, data, transpose = TRUE)
      ZCinvZ <- t(tmp2) %*% tmp2
      XCinvZ <- t(Xmat) %*% backsolve(cov.chol, tmp2)
      tmp3 <- backsolve(xcx.chol, XCinvZ, transpose = TRUE)
      qf <- t(tmp3) %*% tmp3
      loglikelihood <- m * sum(log(diag(cov.chol))) + m * sum(log(diag(xcx.chol))) + 0.5 * sum(diag(ZCinvZ)) - 0.5 * sum(diag(qf))
    }

    if (abs(loglikelihood) == Inf) { loglikelihood <- 1e+06 } # Make sure not Inf
    return(loglikelihood)
  }

}


