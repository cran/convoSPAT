#======================================================================================
# Local likelihood estimation for covariance functions with spatially-varying
# parameters: the convoSPAT() package for R
# Mark D. Risser / The Ohio State University / 2014-2015
#======================================================================================

#======================================================================================
# Summary/Plotting Functions
#======================================================================================


#======================================================================================
# Function to print parameter estimates from NSconvo.fit()
#======================================================================================
#ROxygen comments ----
#' Summarize the nonstationary model fit.
#'
#' \code{summary_NS} prints relevant output from the model fitting
#' procedure.
#'
#' @param fit.model An object which is the output of
#' \code{NSconvo.fit}.
#'
#' @return Text containing the model fitting results.
#'
#' @examples
#' \dontrun{
#' summary_NS( fit.model )
#' }
#'
#' @export

summary_NS <- function( fit.model ){

  cat("Locally estimated mean and variance parameters: \n")
  print( round(fit.model$MLEs.save,5) )
  cat("\n")

  cat( "Estimate of the mean coefficients: ", round(fit.model$beta.GLS,4), "\n" )
  cat("\n")
  cat( "Regression table for mean coefficients: \n")
  print( round(fit.model$Mean.coefs, 5) )
  cat("\n")

  #===================================
  # Print the nugget estimate, for:

  # 1. The nonstationary model
  if(is.null(fit.model$ns.nugget) == FALSE){
    if( fit.model$ns.nugget == FALSE ){
      cat( "Estimate of the nugget variance: ", round( fit.model$tausq.est, 4), "\n" )
    }
    if( fit.model$ns.nugget == TRUE ){
      cat( "Spatially-varying nugget variance. \n Average nugget variance: ",
           round( mean(fit.model$tausq.est), 4), "\n")
    }
  }
  # 2. The stationary model
  if(is.null(fit.model$ns.nugget) == TRUE){
    cat( "Estimate of the nugget variance: ", round( fit.model$tausq.est, 4), "\n" )
  }

  #===================================
  # Print the process variance estimate, for:

  # 1. The nonstationary model
  if(is.null(fit.model$ns.variance) == FALSE){
    if( fit.model$ns.variance == FALSE ){
      cat( "Estimate of the process variance: ", round(fit.model$sigmasq.est,4), "\n" )
    }
    if( fit.model$ns.variance == TRUE ){
      cat( "Spatially-varying process variance. \n Average process variance: ",
           round(mean(fit.model$sigmasq.est), 4), "\n")
    }
  }
  # 2. The stationary model
  if(is.null(fit.model$ns.variance) == TRUE){
    cat( "Estimate of the process variance: ", round(fit.model$sigmasq.est,4), "\n" )
  }
  cat( "Estimate of the smoothness: ", round(fit.model$kappa.MLE, 4) )

}


#======================================================================================
# Function to calculate the sample size for each mixture component location
# for a particular mixture component grid and fit radius
#======================================================================================
#ROxygen comments ----
#' Calculate local sample sizes.
#'
#' \code{mc_N} calculates the number of observations (sample size) that
#' fall within a certain fit radius for each mixture component location.
#'
#' @param coords A matrix of observation locations.
#' @param mc.locations A matrix of the mixture component locations to
#' use in the model fitting.
#' @param fit.radius Scalar; defines the fitting radius for local likelihood
#' estimation.
#'
#' @return A vector \code{mc.N.fit}, which summarizes the number of
#' observation locations in \code{coords} that fall within the fit radius
#' for each mixture component location.
#'
#' @examples
#' \dontrun{
#' mc_N( coords = simdata$sim.locations, mc.locations = simdata$mc.locations,
#' fit.radius = 1 )
#' }
#'
#' @export

mc_N <- function( coords, mc.locations, fit.radius ){


  K <- dim(mc.locations)[1]
  mc.N.fit <- rep(NA,K)

  for( k in 1:K ){

    temp.locs <- coords[ abs(coords[,1]-mc.locations[k,1]) <= fit.radius
                         & (abs(coords[,2] - mc.locations[k,2]) <= fit.radius), ]

    # Isolate the data/locations to be used for calculating the local kernel
    distances <- rep(NA,dim(temp.locs)[1])

    for(i in 1:dim(temp.locs)[1]){
      distances[i] <- sqrt(sum((temp.locs[i,] - mc.locations[k,])^2))
    }

    temp.locations <- temp.locs[distances <= fit.radius,]
    n.fit <- dim(temp.locations)[1]

    mc.N.fit[k] <- n.fit

  }

  return(mc.N.fit)


}


#======================================================================================
# Evaluation criteria
#======================================================================================
#ROxygen comments ----
#' Evaluation criteria
#'
#' Calculate two evaluation criteria -- continuous rank probability score
#' (CRPS) and mean squared prediction error (MSPE) -- comparing hold-out
#' data and predictions.
#'
#' @param holdout.data Observed/true data that has been held out for model
#' comparison.
#' @param pred.mean Predicted mean values corresponding to the hold-out
#' locations.
#' @param pred.SDs Predicted standard errors corresponding to the hold-out
#' locations.
#'
#' @return A list with the following components:
#' \item{CRPS}{The CRPS averaged over all hold-out locations.}
#' \item{MSPE}{The mean squared prediction error.}
#'
#' @examples
#' \dontrun{
#' evaluate_CV( holdout.data = simdata$sim.data[holdout.index],
#' pred.mean = pred.NS$pred.means, pred.SDs = pred.NS$pred.SDs )
#' }
#'
#' @export
#' @importFrom stats dnorm
#' @importFrom stats pnorm

evaluate_CV <- function( holdout.data, pred.mean, pred.SDs ){

  M <- length(holdout.data)

  # CRPS
  CRPS_out <- rep(NA,M)

  for(m in 1:M){
    sigma <- pred.SDs[m]
    zscore <- (holdout.data[m] - pred.mean[m])/sigma

    CRPS_out[m] <- sigma*( (1/sqrt(pi)) - 2*dnorm(zscore) - zscore*( 2*pnorm(zscore) - 1 ) )
  }

  # MSPE
  MSPE <- mean( (holdout.data - pred.mean)^2 )

  # Output
  output <- list( CRPS = mean(CRPS_out),
                  MSPE = MSPE )
  return(output)
}


#======================================================================================
# Plots from the nonstationary model
#======================================================================================
#ROxygen comments ----
#' Plot estimated anisotropy ellipses.
#'
#' This function plots a subset of the estimated anisotropy ellipses
#' over the spatial domain. This function works only when the locations
#' are on a grid.
#'
#' @param locations A N x 2 matrix of observed locations.
#' @param kernels A N x 2 x 2 array of estimated anisotropy ellipses
#' corresponding to the observed locations.
#' @param length Scalar; how many ellipses should be plotted in
#' each direction.
#'
#' @return A plot of estimated ellipses is printed.
#'
#' @examples
#' \dontrun{
#' plot_ellipses( locations = simdata$coords,
#' kernels = fit.model$kernel.ellipses, length = 4 )
#' }
#'
#' @export
#' @importFrom ellipse ellipse
#' @importFrom graphics lines
#' @importFrom graphics plot

plot_ellipses <- function( locations, kernels, length ){

  # Set up a grid to plot a subset of the ellipses
  xmin <- min(locations[,1])
  xmax <- max(locations[,1])
  ymin <- min(locations[,2])
  ymax <- max(locations[,2])

  N <- dim(locations)[1]

  x_plot <- seq(from = xmin, to = xmax, length = length )
  y_plot <- seq(from = ymin, to = ymax, length = length )
  plot_locs <- expand.grid(x_plot, y_plot)
  plot_locs <- matrix(c(plot_locs[,1], plot_locs[,2]), ncol=2, byrow=F)
  K <- dim(plot_locs)[1]

  distances <- matrix(NA, N, K)
  for(n in 1:N){
    for(k in 1:K){
      distances[n,k] <- sqrt(sum((locations[n,] - plot_locs[k,])^2))
    }
  }

  plot_seq <- rep(NA,K)
  for(k in 1:K){
    plot_seq[k] <- which.min(distances[,k])
  }

  # Plot the ellipses
  plot( ellipse( kernels[,,plot_seq[1]], centre = locations[plot_seq[1],], level = 0.5),
        type="l", xlab="Longitude", ylab="Latitude", main="Estimated kernels", asp=1,
        xlim=c(xmin, xmax), ylim=c(ymin, ymax))
  for(k in 2:K){
    lines( ellipse( kernels[,,plot_seq[k]], centre = locations[plot_seq[k],], level = 0.5) )
  }

}

#ROxygen comments ----
#' Plot kriging predictions and standard errors.
#'
#' This function uses \code{\link[fields]{image.plot}} and
#' \code{\link[fields]{quilt.plot}} to quickly plot predictions. If the
#' predicted values are on a rectangular grid, \code{image.plot} is
#' used to improve the appearance.
#'
#' @param locations A N x 2 matrix of prediction locations.
#' @param predmeans Vector; predicted mean values corresponding to the
#' prediction locations.
#' @param predSDs Vector; predicted standard errors corresponding to the
#' prediction locations.
#' @param main1,main2 String; specifies plot titles for panels 1 and 2.
#' @param grid.locations Logical; indicates if the prediction locations
#' are on a rectangular grid (\code{TRUE}) or not (\code{FALSE}).
#'
#' @return A plot with two panels is printed, showing the predictions and
#' standard errors.
#'
#' @examples
#' \dontrun{
#' plot_preds( locations = pred.locations, predmeans = pred.NS$pred.means,
#' predSDs = pred.NS$pred.SDs, grid.locations = TRUE )
#' }
#'
#' @export
#' @importFrom fields image.plot
#' @importFrom fields quilt.plot
#' @importFrom graphics par
#'
plot_preds <- function( locations, predmeans, predSDs, main1 = "Predictions",
                        main2 = "Prediction standard errors", grid.locations = TRUE ){

  if( grid.locations == TRUE ){
    grid_x <- unique( locations[,1] )
    grid_y <- unique( locations[,2] )
    par(mfrow=c(1,2))
    image.plot(grid_x, grid_y, matrix(predmeans, length(grid_x), length(grid_y)), main=main1,
               xlab="Longitude", ylab="Latitude", asp=1)
    image.plot(grid_x, grid_y, matrix(predSDs, length(grid_x), length(grid_y)), main=main2,
               xlab="Longitude", ylab="Latitude", asp=1)
  }
  if( grid.locations == FALSE ){
    par(mfrow=c(1,2))
    quilt.plot( locations, predmeans, xlab="Longitude", ylab="Latitude", main = main1, asp=1 )
    quilt.plot( locations, predSDs, xlab="Longitude", ylab="Latitude", main = main2, asp=1 )
  }

}

#ROxygen comments ----
#' Plot estimated anisotropy ellipses.
#'
#' This function plots the estimated anisotropy ellipses for each
#' of the mixture component locations.
#'
#' @param mc.locations A matrix of mixture component locations.
#' @param mc.kernels An array of kernel matrices corresponding to each
#' of the mixture component locations.
#' @param fit.radius Scalar; defines the fit radius used for the local
#' likelihood estimation.
#' @param aniso.mat 2 x 2 matrix; contains the estimated anisotropy
#' ellipse from the stationary model (for comparison).
#'
#' @return A plot of the estimated ellipses is printed.
#'
#' @examples
#' \dontrun{
#' plot.mc( mc.locations = simdata$mc.locations,
#' mc.kernels = fit.model$mc.kernels, fit.radius = 1 )
#' }
#'
#' @export
#' @importFrom ellipse ellipse
#' @importFrom plotrix draw.circle
#' @importFrom graphics lines
#' @importFrom graphics plot

plot_mc <- function( mc.locations, mc.kernels, fit.radius, aniso.mat = NULL ){

  xmin <- min(mc.locations[,1])
  xmax <- max(mc.locations[,1])
  ymin <- min(mc.locations[,2])
  ymax <- max(mc.locations[,2])

  K <- dim(mc.locations)[1]

  # Plot the ellipses
  plot( ellipse( mc.kernels[,,1], centre = mc.locations[1,], level = 0.5),
        type="l", xlab="Longitude", ylab="Latitude", main="Estimated ellipses for NS (red) \n and S (blue) models", asp=1,
        xlim=c(xmin - 0.25*(xmax-xmin), xmax + 0.25*(xmax-xmin)),
        ylim=c(ymin - 0.25*(ymax-ymin), ymax + 0.25*(ymax-ymin)),
        col = 2 )
  if( is.null(aniso.mat) == FALSE ){
    lines( ellipse(aniso.mat, centre = mc.locations[1,], level = 0.5), lty="dotted",
           col = 4 )
  }
  draw.circle(mc.locations[1,1], mc.locations[1,2], fit.radius, lty="dashed")

  for(k in 2:K){
    lines( ellipse( mc.kernels[,,k], centre = mc.locations[k,], level = 0.5), col=2 )
    plotrix::draw.circle(mc.locations[k,1], mc.locations[k,2], fit.radius, lty="dashed")
    if( is.null(aniso.mat) == FALSE ){
      lines( ellipse(aniso.mat, centre = mc.locations[k,], level = 0.5), lty="dotted",
             col = 4 )
    }

  }

}


#ROxygen comments ----
#' Plot estimated correlations.
#'
#' This function plots the estimated correlation over the spatial
#' domain between a reference point and all other prediction locations.
#'
#' @param model.fit.obj The object which is the output of either
#' \code{NSconvo.fit} or \code{Aniso.fit}.
#' @param ref.loc Vector of length 2; the reference location.
#' @param all.pred.locs A matrix of all prediction locations.
#' @param NS.model Logical; indicates if \code{model.fit.obj} is
#' from \code{NSconvo.fit} (\code{TRUE}) or \code{Aniso.fit}
#' (\code{FALSE}).
#' @param grid Logical; indicates if the \code{all.pred.locs}
#' are on a rectangular grid (\code{TRUE}) or not (\code{FALSE}).
#'
#' @return A plot of the estimated correlation is printed.
#'
#' @examples
#' \dontrun{
#' plot_correlation( model.fit.obj = fit.model, ref.loc = c(1.5, 3.5),
#' all.pred.locs = pred.locs, NS.model = TRUE )
#' }
#'
#' @export
#' @importFrom fields image.plot
#' @importFrom fields quilt.plot

plot_correlation <- function( model.fit.obj, ref.loc, all.pred.locs, NS.model = TRUE, grid = TRUE ){

  M <- dim(all.pred.locs)[1]

  if( NS.model == TRUE ){

    # mc kernels and locations
    mc.kern <- model.fit.obj$mc.kernels
    mc.loc <- model.fit.obj$mc.locations
    K <- dim(mc.loc)[1]
    lambda.w <- model.fit.obj$lambda.w

    # Smoothness
    kappa <- model.fit.obj$kappa.MLE

    cov.model <- model.fit.obj$cov.model

    all.pred.weights <- matrix(NA, M, K)
    for(m in 1:M){
      for(k in 1:K){
        all.pred.weights[m,k] <- exp(-sum((all.pred.locs[m,] - mc.loc[k,])^2)/(2*lambda.w))
      }
      # Normalize the weights
      all.pred.weights[m,] <- all.pred.weights[m,]/sum(all.pred.weights[m,])
    }
    pred.weight <- rep(NA, K)
    for(k in 1:K){
      pred.weight[k] <- exp(-sum((ref.loc - mc.loc[k,])^2)/(2*lambda.w))
    }
    # Normalize the weights
    pred.weight <- pred.weight/sum(pred.weight)


    #=======================================
    # Calculate the prediction kernel ellipses
    pred.kernel.ellipses <- array(0, dim=c(2,2,M))

    for(m in 1:M){
      for(k in 1:K){
        pred.kernel.ellipses[,,m] <- pred.kernel.ellipses[,,m] + all.pred.weights[m,k]*mc.kern[,,k]
      }
    }
    pred.kernel.ellipse <- matrix( rep(0,4), nrow=2, ncol=2 )
    for(k in 1:K){
      pred.kernel.ellipse <- pred.kernel.ellipse + pred.weight[k]*mc.kern[,,k]
    }

    Scale.cross <- rep(NA, M)
    Dist.cross <- rep(NA, M)

    # Calculate the elements of the Mx1 cross-correlation matrix.

    # Diagonal elements
    Kerneli <- pred.kernel.ellipse
    det_i <- Kerneli[1,1]*Kerneli[2,2] - Kerneli[1,2]*Kerneli[2,1]
    Ui <- chol(Kerneli)
    for(j in 1:M){

      Kernelj <- pred.kernel.ellipses[,,j]
      det_j <- Kernelj[1,1]*Kernelj[2,2] - Kernelj[1,2]*Kernelj[2,1]

      avg_ij <- 0.5 * (Kerneli + Kernelj)
      Uij <- chol(avg_ij)
      det_ij <- avg_ij[1,1]*avg_ij[2,2] - avg_ij[1,2]*avg_ij[2,1]
      vec_ij <- backsolve(Uij, (ref.loc - all.pred.locs[j,]), transpose = TRUE)

      Scale.cross[j] <- sqrt( sqrt(det_i*det_j) / det_ij )
      Dist.cross[j] <- sqrt(sum(vec_ij^2))

    }

    Unscl.cross <- geoR::cov.spatial( Dist.cross, cov.model = cov.model,
                                cov.pars = c(1,1), kappa = kappa )
    Cov <- matrix(Scale.cross * Unscl.cross, ncol=1)


  }

  #==========================
  # Stationary model
  #==========================
  if( NS.model == FALSE ){

    # Kernel parameters
    lam1 <- model.fit.obj$aniso.pars[1]
    lam2 <- model.fit.obj$aniso.pars[2]
    eta <- model.fit.obj$aniso.pars[3]

    # Smoothness
    kappa <- model.fit.obj$kappa.MLE

    cov.model <- model.fit.obj$cov.model

    # Calculate the correlations
    Pmat <- matrix(c(cos(eta),-sin(eta),sin(eta),cos(eta)),nrow=2,byrow=T)
    Dmat <- diag(c(lam1,lam2))

    Sigma <- Pmat %*% Dmat %*% t(Pmat)

    distances <- StatMatch::mahalanobis.dist( data.x = all.pred.locs,
                                   data.y = t(ref.loc),
                                   vc = Sigma )
    Cov <- cov.spatial(distances, cov.model = cov.model, cov.pars = c(1,1), kappa = kappa)
    Cov <- matrix(Cov, ncol=1)
  }


  # Plot
  if( grid == TRUE ){
  grid_x <- unique( all.pred.locs[,1] )
  grid_y <- unique( all.pred.locs[,2] )
  image.plot(grid_x, grid_y, matrix(Cov, length(grid_x), length(grid_y)),
             xlab="Longitude", ylab="Latitude", asp=1,
             main = "Estimated correlations")
  }
  if( grid == FALSE){
  	quilt.plot(all.pred.locs, c(Cov), asp=1)
  }


}

