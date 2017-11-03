#======================================================================================
# Local likelihood estimation for covariance functions with spatially-varying
# parameters: the convoSPAT() package for R
# Mark D. Risser / The Ohio State University / 2014-2015
#======================================================================================

#======================================================================================
# Fit/Predict Functions
#======================================================================================

#======================================================================================
# Fit the nonstationary model
#======================================================================================
# The NSconvo.fit() function estimates the parameters of the nonstationary
# convolution-based spatial model. Required inputs are the observed data and
# locations (a geoR object with $coords and $data).
# Optional inputs include mixture component locations (if not provided, the number of mixture component
# locations are required), the fit radius, the covariance model (exponential is
# the default), and whether or not the nugget and process variance
# will be spatially-varying.
# The output of the model is the mixture component locations, mixture component kernels, estimates of
# mu (mean), tausq (nugget variance), sigmasq (process variance), local MLEs,
# the covariance model, and the MLE covariance matrix.
#======================================================================================
#ROxygen comments ----
#' Fit the nonstationary spatial model
#'
#' \code{NSconvo_fit} estimates the parameters of the nonstationary
#' convolution-based spatial model. Required inputs are the observed data and
#' locations (a geoR object with $coords and $data).
#' Optional inputs include mixture component locations (if not provided,
#' the number of mixture component locations are required), the fit radius,
#' the covariance model (exponential is the default), and whether or not the
#' nugget and process variance will be spatially-varying.
#'
#' @param geodata A list containing elements \code{coords} and \code{data} as
#' described next. Typically an object of the class "\code{geodata}", although
#' a geodata object only allows \code{data} to be a vector (no replicates).
#' If not provided, the arguments \code{coords} and \code{data} must be
#' provided instead.
#' @param sp.SPDF A "\code{SpatialPointsDataFrame}" object, which contains the
#' spatial coordinates and additional attribute variables corresponding to the
#' spatoal coordinates
#' @param coords An N x 2 matrix where each row has the two-dimensional
#' coordinates of the N data locations. By default, it takes the \code{coords}
#' component of the argument \code{geodata}, if provided.
#' @param data A vector or matrix with N rows, containing the data values.
#' Inputting a vector corresponds to a single replicate of data, while
#' inputting a matrix corresponds to replicates. In the case of replicates,
#' the model assumes the replicates are independent and identically
#' distributed.
#' @param cov.model A string specifying the model for the correlation
#' function; following \code{geoR}, defaults to \code{"exponential"}.
#' Options available in this package are: "\code{exponential}",
#' \code{"cauchy"}, \code{"matern"}, \code{"circular"}, \code{"cubic"},
#' \code{"gaussian"}, \code{"spherical"}, and \code{"wave"}. For further
#' details, see documentation for \code{\link[geoR]{cov.spatial}}.
#' @param mean.model An object of class \code{\link[stats]{formula}},
#' specifying the mean model to be used. Defaults to an intercept only.
#' @param mc.locations Optional; matrix of mixture component locations.
#' @param N.mc Optional; if \code{mc.locations} is not specified, the
#' function will create a rectangular grid of size \code{N.mc} over the
#' spatial domain.
#' @param lambda.w Scalar; tuning parameter for the weight function.
#' Defaults to be the square of one-half of the minimum distance between
#' mixture component locations.
#' @param fixed.nugg2.var Optional; describes the variance/covariance for
#' a fixed (second) nugget term (represents a known error term). Either
#' a vector of length N containing a station-specific variances (implying
#' independent error) or an NxN covariance matrix (implying dependent error).
#' Defaults to zero.
#' @param mean.model.df Optional data frame; refers to the variables used
#' in \code{mean.model}. Important when using categorical variables in
#' \code{mean.model}, as a subset of the full design matrix will likely
#' be rank deficient. Specifying \code{mean.model.df} allows \code{NSconvo_fit}
#' to calculate a design matrix specific to the points used to fit each
#' local model.
#' @param mc.kernels Optional specification of mixture component kernel
#' matrices (based on expert opinion, etc.).
#' @param fit.radius Scalar; specifies the fit radius or neighborhood size
#' for the local likelihood estimation.
#' @param ns.nugget Logical; indicates if the nugget variance (tausq) should
#' be spatially-varying (\code{TRUE}) or constant (\code{FALSE}).
#' @param ns.variance Logical; indicates if the process variance (sigmasq)
#' should be spatially-varying (\code{TRUE}) or constant (\code{FALSE}).
#' @param ns.mean Logical; indicates if the mean coefficeints (beta)
#' should be spatially-varying (\code{TRUE}) or constant (\code{FALSE}).
#' @param local.aniso Logical; indicates if the local covariance should be
#' anisotropic (\code{TRUE}) or isotropic (\code{FALSE}). Defaults to \code{TRUE}.
#' In the case of a locally isotropic model, the bounds and initial values
#' for lam will default to the first element of \code{local.pars.LB},
#' \code{local.pars.UB}, and \code{local.ini.pars} (while still required, the
#' second and third elements of these vectors will be ignored.)
#' @param fix.kappa Logical; indicates if the kappa parameter should be
#' fixed (\code{TRUE}) or estimated (\code{FALSE}). Defaults to \code{FALSE}
#' (only valid for \code{cov.model = "matern"} and \code{cov.model = "cauchy"}).
#' @param kappa Scalar; value of the kappa parameter. Only used if
#' \code{fix.kappa = TRUE}.
#' @param fix.tausq Logical; indicates whether the default nugget term
#' (tau^2) should be fixed (\code{TRUE}) or estimated (\code{FALSE}). Defaults to
#' \code{FALSE}.
#' @param tausq Scalar; fixed value for the nugget variance (when
#' \code{fix.tausq = TRUE}).
#' @param method Indicates the estimation method, either maximum likelihood
#' (\code{"ml"}) or restricted maximum likelihood (\code{"reml"}).
#' @param print.progress Logical; if \code{TRUE}, text indicating the progress
#' of local model fitting in real time.
#'
#' @param local.pars.LB,local.pars.UB Optional vectors of lower and upper
#' bounds, respectively, used by the \code{"L-BFGS-B"} method option in the
#' \code{\link[stats]{optim}} function for the local parameter estimation.
#' Each vector must be of length five,
#' containing values for lam1, lam2, tausq, sigmasq, and nu. Default for
#' \code{local.pars.LB} is \code{rep(1e-05,5)}; default for
#' \code{local.pars.UB} is \code{c(max.distance/2, max.distance/2, 4*resid.var, 4*resid.var, 100)},
#' where \code{max.distance} is the maximum interpoint distance of the
#' observed data and \code{resid.var} is the residual variance from using
#' \code{\link[stats]{lm}} with \code{mean.model}.
#'
#' @param global.pars.LB,global.pars.UB Optional vectors of lower and upper
#' bounds, respectively, used by the \code{"L-BFGS-B"} method option in the
#' \code{\link[stats]{optim}} function for the global parameter estimation.
#' Each vector must be of length three,
#' containing values for tausq, sigmasq, and nu. Default for
#' \code{global.pars.LB} is \code{rep(1e-05,3)}; default for
#' \code{global.pars.UB} is \code{c(4*resid.var, 4*resid.var, 100)},
#' where \code{resid.var} is the residual variance from using
#' \code{\link[stats]{lm}} with \code{mean.model}.
#'
#' @param local.ini.pars Optional vector of initial values used by the
#' \code{"L-BFGS-B"} method option in the \code{\link[stats]{optim}}
#' function for the local parameter estimation. The vector must be of length
#' five, containing values for lam1, lam2, tausq, sigmasq, and nu. Defaults
#' to \code{c(max.distance/10, max.distance/10, 0.1*resid.var, 0.9*resid.var, 1)},
#' where \code{max.distance} is the maximum interpoint distance of the
#' observed data and \code{resid.var} is the residual variance from using
#' \code{\link[stats]{lm}} with \code{mean.model}.
#'
#' @param global.ini.pars Optional vector of initial values used by the
#' \code{"L-BFGS-B"} method option in the \code{\link[stats]{optim}}
#' function for the global parameter estimation. The vector must be of length
#' three, containing values for tausq, sigmasq, and nu. Defaults to
#' \code{c(0.1*resid.var, 0.9*resid.var, 1)}, where \code{resid.var} is the
#' residual variance from using \code{\link[stats]{lm}} with \code{mean.model}.
#'
#' @return A "NSconvo" object, with the following components:
#' \item{mc.locations}{Mixture component locations used for the simulated
#' data.}
#' \item{mc.kernels}{Mixture component kernel matrices used for the simulated
#' data.}
#' \item{MLEs.save}{Table of local maximum likelihood estimates for each
#' mixture component location.}
#' \item{kernel.ellipses}{\code{N.obs} x 2 x 2 array, containing the kernel
#' matrices corresponding to each of the simulated values.}
#' \item{data}{Observed data values.}
#' \item{beta.GLS}{Generalized least squares estimates of beta,
#' the mean coefficients. For \code{ns.mean = FALSE}, this is a vector
#' (containing the global mean coefficients); for \code{ns.mean = TRUE},
#' this is a matrix (one column for each mixture component location).}
#' \item{beta.cov}{Covariance matrix of the generalized least squares
#' estimate of beta. For \code{ns.mean = FALSE}, this is a matrix
#' (containing the covariance of theglobal mean coefficients); for
#' \code{ns.mean = TRUE}, this is an array (one matrix for each mixture
#' component location).}
#' \item{Mean.coefs}{"Regression table" for the mean coefficient estimates,
#' listing the estimate, standard error, and t-value (for \code{ns.mean =
#' FALSE} only).}
#' \item{tausq.est}{Estimate of tausq (nugget variance), either scalar (when
#' \code{ns.nugget = "FALSE"}) or a vector of length N (when
#' \code{ns.nugget = "TRUE"}), which contains the estimated nugget variance
#' for each observation location.}
#' \item{sigmasq.est}{Estimate of sigmasq (process variance), either scalar
#' (when \code{ns.variance = "FALSE"}) or a vector of length N (when
#' \code{ns.variance = "TRUE"}), which contains the estimated process
#' variance for each observation location.}
#' \item{beta.est}{Estimate of beta (mean coefficients), either a vector
#' (when \code{ns.mean = "FALSE"}) or a matrix with N rows (when
#' \code{ns.mean = "TRUE"}), each row of which contains the estimated
#' (smoothed) mean coefficients for each observation location.}
#' \item{kappa.MLE}{Scalar maximum likelihood estimate for kappa (when
#' applicable).}
#' \item{Cov.mat}{Estimated covariance matrix (\code{N.obs} x \code{N.obs})
#' using all relevant parameter estimates.}
#' \item{Cov.mat.chol}{Cholesky of \code{Cov.mat} (i.e., \code{chol(Cov.mat)}),
#' the estimated covariance matrix (\code{N.obs} x \code{N.obs}).}
#' \item{cov.model}{String; the correlation model used for estimation.}
#' \item{ns.nugget}{Logical, indicating if the nugget variance was estimated
#' as spatially-varing (\code{TRUE}) or constant (\code{FALSE}).}
#' \item{ns.variance}{Logical, indicating if the process variance was
#' estimated as spatially-varying (\code{TRUE}) or constant (\code{FALSE}).}
#' \item{fixed.nugg2.var}{N x N matrix with the fixed
#' variance/covariance for the second (measurement error) nugget term (defaults
#' to zero).}
#' \item{coords}{N x 2 matrix of observation locations.}
#' \item{global.loglik}{Scalar value of the maximized likelihood from the
#' global optimization (if available).}
#' \item{Xmat}{Design matrix, obtained from using \code{\link[stats]{lm}}
#' with \code{mean.model}.}
#' \item{lambda.w}{Tuning parameter for the weight function.}
#' \item{fix.kappa}{Logical, indicating if kappa was fixed (\code{TRUE}) or
#' estimated (\code{FALSE}).}
#' \item{kappa}{Scalar; fixed value of kappa.}
#'
#' @examples
#' # Using white noise data
#' fit.model <- NSconvo_fit( coords = cbind( runif(100), runif(100)),
#' data = rnorm(100), fit.radius = 0.4, N.mc = 4 )
#'
#' @export
#' @importFrom geoR cov.spatial
#' @importFrom stats lm
#' @importFrom stats optim
#' @importFrom stats dist

NSconvo_fit <- function( geodata = NULL, sp.SPDF = NULL,
                         coords = geodata$coords, data = geodata$data,
                         cov.model = "exponential", mean.model = data ~ 1,
                         mc.locations = NULL, N.mc = NULL, lambda.w = NULL,
                         fixed.nugg2.var = NULL, mean.model.df = NULL,
                         mc.kernels = NULL, fit.radius = NULL,
                         ns.nugget = FALSE, ns.variance = FALSE,
                         ns.mean = FALSE, local.aniso = TRUE,
                         fix.tausq = FALSE, tausq = 0,
                         fix.kappa = FALSE, kappa = 0.5,
                         method = "reml", print.progress = TRUE,
                         local.pars.LB = NULL, local.pars.UB = NULL,
                         global.pars.LB = NULL, global.pars.UB = NULL,
                         local.ini.pars = NULL, global.ini.pars = NULL ){

  if( is.null(fit.radius) ){cat("\nPlease specify a fitting radius.\n")}
  #===========================================================================
  # Formatting for coordinates/data
  #===========================================================================
  if( is.null(geodata) == FALSE ){
    if( class(geodata) != "geodata" ){
      stop("Please use a geodata object for the 'geodata = ' input.")
    }
    coords <- geodata$coords
    data <- geodata$data
  }
  if( is.null(sp.SPDF) == FALSE ){
    if( class(sp.SPDF) != "SpatialPointsDataFrame" ){
      stop("Please use a SpatialPointsDataFrame object for the 'sp.SPDF = ' input.")
    }
    geodata <- geoR::as.geodata( sp.SPDF )
    coords <- geodata$coords
    data <- geodata$data
  }

  coords <- as.matrix(coords)
  N <- dim(coords)[1]
  data <- as.matrix(data, nrow=N)
  p <- dim(data)[2]

  # Make sure cov.model is one of the permissible options
  if( cov.model != "cauchy" & cov.model != "matern" & cov.model != "circular" &
        cov.model != "cubic" & cov.model != "gaussian" & cov.model != "exponential" &
        cov.model != "spherical" & cov.model != "wave" ){
    stop("Please specify a valid covariance model (cauchy, matern,\ncircular, cubic, gaussian,
          exponential, spherical, or wave).")
  }

  # Check that ns.mean = TRUE is only used where applicable
  if( ns.mean == TRUE ){
    if( ns.nugget == FALSE || ns.variance == FALSE ){
      stop("Cannot use ns.mean = TRUE and either ns.nugget = FALSE or ns.variance = FALSE (currently unsupported).")
    }
  }

  # Check that fix.tausq = TRUE is only used where applicable
  if( fix.tausq == TRUE ){
    if( ns.nugget == FALSE || ns.variance == FALSE ){
      stop("Cannot use fix.tausq == TRUE and either ns.nugget = FALSE or ns.variance = FALSE (currently unsupported).")
    }
  }

  #===========================================================================
  # Calculate the mixture component locations if not user-specified
  #===========================================================================
  if( is.null(mc.locations) == TRUE ){ # Calculate the mixture component locations

    if( is.null(N.mc) == TRUE ){
      cat("Please enter the desired number of mixture component locations. \n")
    }

    lon_min <- min(coords[,1])
    lon_max <- max(coords[,1])
    lat_min <- min(coords[,2])
    lat_max <- max(coords[,2])

    #=======================================
    # mixture component knot locations
    #=======================================
    mc_x <- seq(from = lon_min + 0.5*(lon_max - lon_min)/floor(sqrt(N.mc)),
                   to = lon_max - 0.5*(lon_max - lon_min)/floor(sqrt(N.mc)),
                   length = floor(sqrt(N.mc)) )
    mc_y <- seq(from = lat_min + 0.5*(lat_max - lat_min)/floor(sqrt(N.mc)),
                   to = lat_max - 0.5*(lat_max - lat_min)/floor(sqrt(N.mc)),
                   length = floor(sqrt(N.mc)) )
    mc.locations <- expand.grid( mc_x, mc_y )
    mc.locations <- matrix(c(mc.locations[,1], mc.locations[,2]), ncol=2, byrow=F)
  }

  K <- dim(mc.locations)[1]

  #===========================================================================
  # Check the mixture component locations
  #===========================================================================
  check.mc.locs <- mc_N( coords, mc.locations, fit.radius )
  cat("\n-----------------------------------------------------------\n")
  cat(paste("Fitting the nonstationary model: ", K, " local models with\nlocal sample sizes ranging between ", min(check.mc.locs),
            " and ", max(check.mc.locs), ".", sep = "" ))
  if( ns.nugget == FALSE & ns.variance == FALSE ){
    cat("\nConstant nugget and constant variance.")
  }
  if( ns.nugget == FALSE & ns.variance == TRUE ){
    cat("\nConstant nugget and spatially-varying variance.")
  }
  if( ns.nugget == TRUE & ns.variance == FALSE ){
    cat("\nSpatially-varying nugget and constant variance.")
  }
  if( ns.nugget == TRUE & ns.variance == TRUE ){
    cat("\nSpatially-varying nugget and spatially-varying variance.")
  }
  if( local.aniso ){
    cat(paste("\nLocally anisotropic ", cov.model, " covariance.", sep = ""))
  }
  if( !local.aniso ){
    cat(paste("\nLocally isotropic ", cov.model, " covariance.", sep = ""))
  }

  if( min(check.mc.locs) < 5 ){cat("\nWARNING: at least one of the mc locations has too few data points.\n")}
  cat("\n-----------------------------------------------------------\n")

  #===========================================================================
  # Set the tuning parameter, if not specified
  #===========================================================================
  if( is.null(lambda.w) == TRUE ){
    lambda.w <- ( 0.5*min(dist(mc.locations)) )^2
  }

  #===========================================================================
  # Set the second nugget variance, if not specified
  #===========================================================================
  if( is.null(fixed.nugg2.var) == TRUE ){
    fixed.nugg2.var <- matrix(0,N,N)
  } else{
    if( !is.matrix(fixed.nugg2.var) ){ # Convert to a covariance matrix if specified as a vector
      fixed.nugg2.var <- diag(fixed.nugg2.var)
    }
  }

  #===========================================================================
  # Calculate the design matrix
  #===========================================================================
  if( is.null(mean.model.df) == TRUE ){
    OLS.model <- lm( mean.model, x=TRUE )
    Xmat <- matrix( unname( OLS.model$x ), nrow=N )
    beta.names <- colnames(OLS.model$x)
  }
  if( is.null(mean.model.df) == FALSE ){
    OLS.model <- lm( mean.model, x=TRUE, data = mean.model.df )
    Xmat <- matrix( unname( OLS.model$x ), nrow=N )
    beta.names <- colnames(OLS.model$x)
  }

  q <- ncol(Xmat)

  #===========================================================================
  # Specify lower, upper, and initial parameter values for optim()
  #===========================================================================
  lon_min <- min(coords[,1])
  lon_max <- max(coords[,1])
  lat_min <- min(coords[,2])
  lat_max <- max(coords[,2])

  max.distance <- sqrt( sum((c(lon_min,lat_min) - c(lon_max,lat_max))^2))

  if( p > 1 ){
    ols.sigma <- NULL
    for(i in 1:length(names(summary(OLS.model)))){
      ols.sigma <- c( ols.sigma, summary(OLS.model)[[i]]$sigma )
    }
    resid.var <- (max(ols.sigma))^2
  }
  if( p == 1 ){
    resid.var <- summary(OLS.model)$sigma^2
  }

  #=================================
  # Lower limits for optim()
  #=================================
  if( is.null(local.pars.LB) == TRUE ){
    lam1.LB <- 1e-05
    lam2.LB <- 1e-05
    tausq.local.LB <- 1e-05
    sigmasq.local.LB <- 1e-05
    kappa.local.LB <- 1e-05
  }
  if( is.null(local.pars.LB) == FALSE ){
    lam1.LB <- local.pars.LB[1]
    lam2.LB <- local.pars.LB[2]
    tausq.local.LB <- local.pars.LB[3]
    sigmasq.local.LB <- local.pars.LB[4]
    kappa.local.LB <- local.pars.LB[5]
  }
  if( is.null(global.pars.LB) == TRUE ){
    tausq.global.LB <- 1e-05
    sigmasq.global.LB <- 1e-05
    kappa.global.LB <- 1e-05
  }
  if( is.null(global.pars.LB) == FALSE ){
    tausq.global.LB <- global.pars.LB[1]
    sigmasq.global.LB <- global.pars.LB[2]
    kappa.global.LB <- global.pars.LB[3]
  }

  #=================================
  # Upper limits for optim()
  #=================================
  if( is.null(local.pars.UB) == TRUE ){
    lam1.UB <- max.distance/4
    lam2.UB <- max.distance/4
    tausq.local.UB <- 4*resid.var
    sigmasq.local.UB <- 4*resid.var
    kappa.local.UB <- 30
  }
  if( is.null(local.pars.UB) == FALSE ){
    lam1.UB <- local.pars.UB[1]
    lam2.UB <- local.pars.UB[2]
    tausq.local.UB <- local.pars.UB[3]
    sigmasq.local.UB <- local.pars.UB[4]
    kappa.local.UB <- local.pars.UB[5]
  }
  if( is.null(global.pars.UB) == TRUE ){
    tausq.global.UB <- 4*resid.var
    sigmasq.global.UB <- 4*resid.var
    kappa.global.UB <- 30
  }
  if( is.null(global.pars.UB) == FALSE ){
    tausq.global.UB <- global.pars.UB[1]
    sigmasq.global.UB <- global.pars.UB[2]
    kappa.global.UB <- global.pars.UB[3]
  }

  #=================================
  # Initial values for optim()
  #=================================

  # Local estimation
  if( is.null(local.ini.pars) == TRUE ){
    lam1.init <- max.distance/10
    lam2.init <- max.distance/10
    tausq.local.init <- 0.1*resid.var
    sigmasq.local.init <- 0.9*resid.var
    kappa.local.init <- 1
  }
  if( is.null(local.ini.pars) == FALSE ){
    lam1.init <- local.ini.pars[1]
    lam2.init <- local.ini.pars[2]
    tausq.local.init <- local.ini.pars[3]
    sigmasq.local.init <- local.ini.pars[4]
    kappa.local.init <- local.ini.pars[5]
  }

  # Global estimation
  if( is.null(global.ini.pars) == TRUE ){
    tausq.global.init <- 0.1*resid.var
    sigmasq.global.init <- 0.9*resid.var
    kappa.global.init <- 1
  }
  if( is.null(global.ini.pars) == FALSE ){
    tausq.global.init <- global.ini.pars[1]
    sigmasq.global.init <- global.ini.pars[2]
    kappa.global.init <- global.ini.pars[3]
  }

  #===========================================================================
  # Calculate the mixture component kernels (and other parameters) if not
  # user-specified
  #===========================================================================
  if( is.null(mc.kernels) == TRUE ){

    # Storage for the mixture component kernels
    mc.kernels <- array(NA, dim=c(2, 2, K))
    MLEs.save <- matrix(NA, K, 7)
    beta.GLS.save <- matrix(NA, K, ncol(Xmat))
    beta.cov.save <- array(NA, dim = c(ncol(Xmat), ncol(Xmat), K))
    beta.coefs.save <- list()

    # Estimate the kernel function for each mixture component location,
    # completely specified by the kernel covariance matrix
    for( k in 1:K ){

      if( is.null(mean.model.df) == TRUE ){
        # Select coordinates in a square
        coords.sub <- (abs(coords[,1] - mc.locations[k,1]) <= fit.radius
                       & (abs(coords[,2] - mc.locations[k,2]) <= fit.radius))
        # Subset
        temp.locs <- coords[coords.sub,]
        temp.dat <- data[coords.sub,]
        temp.n2.var <- fixed.nugg2.var[coords.sub, coords.sub]
        X.tem <- as.matrix(Xmat[coords.sub,])

        # Isolate the data/locations to be used for calculating the local kernel
        distances <- rep(NA,dim(temp.locs)[1])

        for(i in 1:dim(temp.locs)[1]){
          distances[i] <- sqrt(sum((temp.locs[i,] - mc.locations[k,])^2))
        }

        temp.locations <- temp.locs[distances <= fit.radius,]
        Xtemp <- X.tem[distances <= fit.radius,]
        n.fit <- dim(temp.locations)[1]
        temp.dat <- as.matrix( temp.dat, nrow=n.fit)
        temp.data <- temp.dat[distances <= fit.radius,]
        temp.data <- as.matrix(temp.data, nrow=n.fit)
        temp.nugg2.var <- temp.n2.var[distances <= fit.radius, distances <= fit.radius]
      }
      if( is.null(mean.model.df) == FALSE ){

        # Select coordinates in a square
        coords.sub <- (abs(coords[,1] - mc.locations[k,1]) <= fit.radius
                       & (abs(coords[,2] - mc.locations[k,2]) <= fit.radius))

        temp.locs <- coords[coords.sub,]
        temp.dat <- data[coords.sub,]
        temp.n2.var <- fixed.nugg2.var[coords.sub, coords.sub]
        temp.mmdf <- mean.model.df[coords.sub, ]

        # Isolate the data/locations to be used for calculating the local kernel
        distances <- rep(NA,dim(temp.locs)[1])

        for(i in 1:dim(temp.locs)[1]){
          distances[i] <- sqrt(sum((temp.locs[i,] - mc.locations[k,])^2))
        }

        temp.locations <- temp.locs[distances <= fit.radius,]
        n.fit <- dim(temp.locations)[1]
        temp.dat <- as.matrix( temp.dat, nrow=n.fit)
        temp.data <- temp.dat[distances <= fit.radius,]
        temp.data <- as.matrix(temp.data, nrow=n.fit)
        temp.nugg2.var <- temp.n2.var[distances <= fit.radius,distances <= fit.radius]
        temp.mmdf <- temp.mmdf[distances <= fit.radius,]
        Xtemp <- matrix( unname( lm( mean.model, x=TRUE, data = temp.mmdf )$x ), nrow=n.fit )
      }

      if( print.progress ){
        if(k == 1){
          cat("Calculating the parameter set for:\n")
        }
        cat("mc location ", k,", using ", n.fit," observations...\n", sep="")
      }

      #####################################################
      # Local estimation
      if( local.aniso ){ # Locally anisotropic
        if( cov.model == "matern" || cov.model == "cauchy" ){ # For covariance models with kappa
          if( !fix.kappa ){ # Estimate kappa
            if( !fix.tausq ){ # Estimate tausq
              if( method == "ml" ){
                f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                            data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                            fixed = rep(FALSE, 6 + q), method = method,
                                            local.aniso = local.aniso, fix.tausq = fix.tausq,
                                            fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
                MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, tausq.local.init,
                                       sigmasq.local.init, kappa.local.init, rep(0, q) ),
                               fn = f_loglik, method = "L-BFGS-B",
                               lower = c( lam1.LB, lam2.LB, 0, tausq.local.LB, sigmasq.local.LB,
                                          kappa.local.LB, rep(-Inf, q) ),
                               upper = c( lam1.UB, lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB,
                                          kappa.local.UB, rep(Inf, q) ) )
                covMLE <- MLEs$par[1:6]
                betaMLE <- MLEs$par[-c(1:6)]
              }
              if( method == "reml" ){
                f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                            data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                            fixed = rep(FALSE, 6), method = method,
                                            local.aniso = local.aniso, fix.tausq = fix.tausq,
                                            fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
                MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, tausq.local.init,
                                       sigmasq.local.init, kappa.local.init ),
                               fn = f_loglik, method = "L-BFGS-B",
                               lower = c( lam1.LB, lam2.LB, 0, tausq.local.LB, sigmasq.local.LB,
                                          kappa.local.LB ),
                               upper = c( lam1.UB, lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB,
                                          kappa.local.UB ) )
                covMLE <- MLEs$par
              }
            }
            if( fix.tausq ){ # Fix tausq
              if( method == "ml" ){
                f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                            data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                            fixed = rep(FALSE, 5 + q), method = method,
                                            local.aniso = local.aniso, fix.tausq = fix.tausq,
                                            fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
                MLEs <- optim( par = c(lam1.init, lam2.init, pi/4,
                                       sigmasq.local.init, kappa.local.init, rep(0, q) ),
                               fn = f_loglik, method = "L-BFGS-B",
                               lower = c( lam1.LB, lam2.LB, 0, sigmasq.local.LB,
                                          kappa.local.LB, rep(-Inf, q) ),
                               upper = c( lam1.UB, lam2.UB, pi/2, sigmasq.local.UB,
                                          kappa.local.UB, rep(Inf, q) ) )
                covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4:5])
                betaMLE <- MLEs$par[-c(1:5)]
              }
              if( method == "reml" ){
                f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                            data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                            fixed = rep(FALSE, 5), method = method,
                                            local.aniso = local.aniso, fix.tausq = fix.tausq,
                                            fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
                MLEs <- optim( par = c(lam1.init, lam2.init, pi/4,
                                       sigmasq.local.init, kappa.local.init ),
                               fn = f_loglik, method = "L-BFGS-B",
                               lower = c( lam1.LB, lam2.LB, 0, sigmasq.local.LB,
                                          kappa.local.LB ),
                               upper = c( lam1.UB, lam2.UB, pi/2, sigmasq.local.UB,
                                          kappa.local.UB ) )
                covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4:5])
              }
            }
          }
          if( fix.kappa ){ # Fix kappa
            if( !fix.tausq ){ # Estimate tausq
              if( method == "ml" ){
                f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                            data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                            fixed = rep(FALSE, 5 + q), method = method,
                                            local.aniso = local.aniso, fix.tausq = fix.tausq,
                                            fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
                MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, tausq.local.init,
                                       sigmasq.local.init, rep(0, q) ),
                               fn = f_loglik, method = "L-BFGS-B",
                               lower = c( lam1.LB, lam2.LB, 0, tausq.local.LB, sigmasq.local.LB,
                                          rep(-Inf, q) ),
                               upper = c( lam1.UB, lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB,
                                          rep(Inf, q) ) )
                covMLE <- c(MLEs$par[1:5], kappa)
                betaMLE <- MLEs$par[-c(1:5)]
              }
              if( method == "reml" ){
                f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                            data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                            fixed = rep(FALSE, 5), method = method,
                                            local.aniso = local.aniso, fix.tausq = fix.tausq,
                                            fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
                MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, tausq.local.init,
                                       sigmasq.local.init ),
                               fn = f_loglik, method = "L-BFGS-B",
                               lower = c( lam1.LB, lam2.LB, 0, tausq.local.LB, sigmasq.local.LB ),
                               upper = c( lam1.UB, lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB ) )
                covMLE <- c(MLEs$par[1:5], kappa)
              }
            }
            if( fix.tausq ){ # Fix tausq
              if( method == "ml" ){
                f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                            data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                            fixed = rep(FALSE, 4 + q), method = method,
                                            local.aniso = local.aniso, fix.tausq = fix.tausq,
                                            fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
                MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, sigmasq.local.init, rep(0, q) ),
                               fn = f_loglik, method = "L-BFGS-B",
                               lower = c( lam1.LB, lam2.LB, 0, sigmasq.local.LB, rep(-Inf, q) ),
                               upper = c( lam1.UB, lam2.UB, pi/2, sigmasq.local.UB, rep(Inf, q) ) )
                covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4], kappa)
                betaMLE <- MLEs$par[-c(1:4)]
              }
              if( method == "reml" ){
                f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                            data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                            fixed = rep(FALSE, 4), method = method,
                                            local.aniso = local.aniso, fix.tausq = fix.tausq,
                                            fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
                MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, sigmasq.local.init ),
                               fn = f_loglik, method = "L-BFGS-B",
                               lower = c( lam1.LB, lam2.LB, 0, sigmasq.local.LB ),
                               upper = c( lam1.UB, lam2.UB, pi/2, sigmasq.local.UB ) )
                covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4], kappa)
              }
            }
          }
        }
        if( cov.model != "matern" & cov.model != "cauchy" ){ # For covariance models without kappa
          if( !fix.tausq ){ # Estimate tausq
            if( method == "ml" ){
              f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                          data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                          fixed = rep(FALSE, 5 + q), method = method,
                                          local.aniso = local.aniso, fix.tausq = fix.tausq,
                                          fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
              MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, tausq.local.init,
                                     sigmasq.local.init, rep(0, q) ),
                             fn = f_loglik, method = "L-BFGS-B",
                             lower = c( lam1.LB, lam2.LB, 0, tausq.local.LB, sigmasq.local.LB,
                                        rep(-Inf, q) ),
                             upper = c( lam1.UB, lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB,
                                        rep(Inf, q) ) )
              covMLE <- c(MLEs$par[1:5], NA)
              betaMLE <- MLEs$par[-c(1:5)]
            }
            if( method == "reml" ){
              f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                          data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                          fixed = rep(FALSE, 5), method = method,
                                          local.aniso = local.aniso, fix.tausq = fix.tausq,
                                          fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
              MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, tausq.local.init,
                                     sigmasq.local.init ),
                             fn = f_loglik, method = "L-BFGS-B",
                             lower = c( lam1.LB, lam2.LB, 0, tausq.local.LB, sigmasq.local.LB ),
                             upper = c( lam1.UB, lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB ) )
              covMLE <- c(MLEs$par[1:5], NA)
            }







          }
          if( fix.tausq ){ # Fix tausq
            if( method == "ml" ){
              f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                          data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                          fixed = rep(FALSE, 4 + q), method = method,
                                          local.aniso = local.aniso, fix.tausq = fix.tausq,
                                          fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
              MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, sigmasq.local.init, rep(0, q) ),
                             fn = f_loglik, method = "L-BFGS-B",
                             lower = c( lam1.LB, lam2.LB, 0, sigmasq.local.LB, rep(-Inf, q) ),
                             upper = c( lam1.UB, lam2.UB, pi/2, sigmasq.local.UB, rep(Inf, q) ) )
              covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4], NA)
              betaMLE <- MLEs$par[-c(1:4)]
            }
            if( method == "reml" ){
              f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                          data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                          fixed = rep(FALSE, 4), method = method,
                                          local.aniso = local.aniso, fix.tausq = fix.tausq,
                                          fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
              MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, sigmasq.local.init ),
                             fn = f_loglik, method = "L-BFGS-B",
                             lower = c( lam1.LB, lam2.LB, 0, sigmasq.local.LB ),
                             upper = c( lam1.UB, lam2.UB, pi/2, sigmasq.local.UB ) )
              covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4], NA)
            }
          }
        }

      }
      if( !local.aniso ){ # Locally isotropic
        if( cov.model == "matern" || cov.model == "cauchy" ){ # For covariance models with kappa
          if( !fix.kappa ){ # Estimate kappa
            if( !fix.tausq ){ # Estimate tausq
              if( method == "ml" ){
                f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                            data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                            fixed = rep(FALSE, 4 + q), method = method,
                                            local.aniso = local.aniso, fix.tausq = fix.tausq,
                                            fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
                MLEs <- optim( par = c(lam1.init, tausq.local.init,
                                       sigmasq.local.init, kappa.local.init, rep(0, q) ),
                               fn = f_loglik, method = "L-BFGS-B",
                               lower = c( lam1.LB, tausq.local.LB, sigmasq.local.LB,
                                          kappa.local.LB, rep(-Inf, q) ),
                               upper = c( lam1.UB, tausq.local.UB, sigmasq.local.UB,
                                          kappa.local.UB, rep(Inf, q) ) )
                covMLE <- c(MLEs$par[1], MLEs$par[1], 0, MLEs$par[2:4] )
                betaMLE <- MLEs$par[-c(1:4)]
              }
              if( method == "reml" ){
                f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                            data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                            fixed = rep(FALSE, 4), method = method,
                                            local.aniso = local.aniso, fix.tausq = fix.tausq,
                                            fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
                MLEs <- optim( par = c(lam1.init, tausq.local.init,
                                       sigmasq.local.init, kappa.local.init ),
                               fn = f_loglik, method = "L-BFGS-B",
                               lower = c( lam1.LB, tausq.local.LB, sigmasq.local.LB,
                                          kappa.local.LB ),
                               upper = c( lam1.UB, tausq.local.UB, sigmasq.local.UB,
                                          kappa.local.UB ) )
                covMLE <- c(MLEs$par[1], MLEs$par[1], 0, MLEs$par[2:4] )
              }
            }
            if( fix.tausq ){ # Fix tausq
              if( method == "ml" ){
                f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                            data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                            fixed = rep(FALSE, 3 + q), method = method,
                                            local.aniso = local.aniso, fix.tausq = fix.tausq,
                                            fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
                MLEs <- optim( par = c(lam1.init, sigmasq.local.init, kappa.local.init, rep(0, q) ),
                               fn = f_loglik, method = "L-BFGS-B",
                               lower = c( lam1.LB, sigmasq.local.LB, kappa.local.LB, rep(-Inf, q) ),
                               upper = c( lam1.UB, sigmasq.local.UB, kappa.local.UB, rep(Inf, q) ) )
                covMLE <- c(MLEs$par[1], MLEs$par[1], 0, tausq, MLEs$par[2:3] )
                betaMLE <- MLEs$par[-c(1:3)]
              }
              if( method == "reml" ){
                f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                            data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                            fixed = rep(FALSE, 3), method = method,
                                            local.aniso = local.aniso, fix.tausq = fix.tausq,
                                            fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
                MLEs <- optim( par = c(lam1.init, sigmasq.local.init, kappa.local.init ),
                               fn = f_loglik, method = "L-BFGS-B",
                               lower = c( lam1.LB, sigmasq.local.LB, kappa.local.LB ),
                               upper = c( lam1.UB, sigmasq.local.UB, kappa.local.UB ) )
                covMLE <- c(MLEs$par[1], MLEs$par[1], 0, tausq, MLEs$par[2:3] )
              }
            }
          }
          if( fix.kappa ){ # Fix kappa
            if( !fix.tausq ){ # Estimate tausq
              if( method == "ml" ){
                f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                            data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                            fixed = rep(FALSE, 3 + q), method = method,
                                            local.aniso = local.aniso, fix.tausq = fix.tausq,
                                            fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
                MLEs <- optim( par = c(lam1.init, tausq.local.init,
                                       sigmasq.local.init, rep(0, q) ),
                               fn = f_loglik, method = "L-BFGS-B",
                               lower = c( lam1.LB, tausq.local.LB, sigmasq.local.LB,
                                          rep(-Inf, q) ),
                               upper = c( lam1.UB, tausq.local.UB, sigmasq.local.UB,
                                          rep(Inf, q) ) )
                covMLE <- c(MLEs$par[1], MLEs$par[1], 0, MLEs$par[2:3], kappa )
                betaMLE <- MLEs$par[-c(1:3)]
              }
              if( method == "reml" ){
                f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                            data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                            fixed = rep(FALSE, 3), method = method,
                                            local.aniso = local.aniso, fix.tausq = fix.tausq,
                                            fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
                MLEs <- optim( par = c(lam1.init, tausq.local.init,
                                       sigmasq.local.init ),
                               fn = f_loglik, method = "L-BFGS-B",
                               lower = c( lam1.LB, tausq.local.LB, sigmasq.local.LB ),
                               upper = c( lam1.UB, tausq.local.UB, sigmasq.local.UB ) )
                covMLE <- c(MLEs$par[1], MLEs$par[1], 0, MLEs$par[2:3], kappa )
              }
            }
            if( fix.tausq ){ # Fix tausq
              if( method == "ml" ){
                f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                            data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                            fixed = rep(FALSE, 2 + q), method = method,
                                            local.aniso = local.aniso, fix.tausq = fix.tausq,
                                            fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
                MLEs <- optim( par = c(lam1.init, sigmasq.local.init, rep(0, q) ),
                               fn = f_loglik, method = "L-BFGS-B",
                               lower = c( lam1.LB, sigmasq.local.LB, rep(-Inf, q) ),
                               upper = c( lam1.UB, sigmasq.local.UB, rep(Inf, q) ) )
                covMLE <- c(MLEs$par[1], MLEs$par[1], 0, tausq, MLEs$par[2], kappa )
                betaMLE <- MLEs$par[-c(1:2)]
              }
              if( method == "reml" ){
                f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                            data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                            fixed = rep(FALSE, 2), method = method,
                                            local.aniso = local.aniso, fix.tausq = fix.tausq,
                                            fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
                MLEs <- optim( par = c(lam1.init, sigmasq.local.init ),
                               fn = f_loglik, method = "L-BFGS-B",
                               lower = c( lam1.LB, sigmasq.local.LB ),
                               upper = c( lam1.UB, sigmasq.local.UB ) )
                covMLE <- c(MLEs$par[1], MLEs$par[1], 0, tausq, MLEs$par[2], kappa )
              }
            }
          }
        }
        if( cov.model != "matern" & cov.model != "cauchy" ){ # For covariance models without kappa
          if( !fix.tausq ){ # Estimate tausq
            if( method == "ml" ){
              f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                          data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                          fixed = rep(FALSE, 3 + q), method = method,
                                          local.aniso = local.aniso, fix.tausq = fix.tausq,
                                          fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
              MLEs <- optim( par = c(lam1.init, tausq.local.init,
                                     sigmasq.local.init, rep(0, q) ),
                             fn = f_loglik, method = "L-BFGS-B",
                             lower = c( lam1.LB, tausq.local.LB, sigmasq.local.LB,
                                        rep(-Inf, q) ),
                             upper = c( lam1.UB, tausq.local.UB, sigmasq.local.UB,
                                        rep(Inf, q) ) )
              covMLE <- c(MLEs$par[1], MLEs$par[1], 0, MLEs$par[2:3], NA )
              betaMLE <- MLEs$par[-c(1:3)]
            }
            if( method == "reml" ){
              f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                          data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                          fixed = rep(FALSE, 3), method = method,
                                          local.aniso = local.aniso, fix.tausq = fix.tausq,
                                          fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
              MLEs <- optim( par = c(lam1.init, tausq.local.init,
                                     sigmasq.local.init ),
                             fn = f_loglik, method = "L-BFGS-B",
                             lower = c( lam1.LB, tausq.local.LB, sigmasq.local.LB ),
                             upper = c( lam1.UB, tausq.local.UB, sigmasq.local.UB ) )
              covMLE <- c(MLEs$par[1], MLEs$par[1], 0, MLEs$par[2:3], NA )
            }
          }
          if( fix.tausq ){ # Fix tausq
            if( method == "ml" ){
              f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                          data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                          fixed = rep(FALSE, 2 + q), method = method,
                                          local.aniso = local.aniso, fix.tausq = fix.tausq,
                                          fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
              MLEs <- optim( par = c(lam1.init, sigmasq.local.init, rep(0, q) ),
                             fn = f_loglik, method = "L-BFGS-B",
                             lower = c( lam1.LB, sigmasq.local.LB, rep(-Inf, q) ),
                             upper = c( lam1.UB, sigmasq.local.UB, rep(Inf, q) ) )
              covMLE <- c(MLEs$par[1], MLEs$par[1], 0, tausq, MLEs$par[2], NA )
              betaMLE <- MLEs$par[-c(1:2)]
            }
            if( method == "reml" ){
              f_loglik <- make_local_lik( locations = temp.locations, cov.model = cov.model,
                                          data = temp.data, Xmat = Xtemp, nugg2.var = temp.nugg2.var,
                                          fixed = rep(FALSE, 2), method = method,
                                          local.aniso = local.aniso, fix.tausq = fix.tausq,
                                          fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
              MLEs <- optim( par = c(lam1.init, sigmasq.local.init ),
                             fn = f_loglik, method = "L-BFGS-B",
                             lower = c( lam1.LB, sigmasq.local.LB ),
                             upper = c( lam1.UB, sigmasq.local.UB ) )
              covMLE <- c(MLEs$par[1], MLEs$par[1], 0, tausq, MLEs$par[2], NA )
            }
          }
        }
      }
      #####################################################

      # Warnings?
      if( MLEs$convergence != 0 ){
        if( MLEs$convergence == 52 ){
          cat( paste("  There was a NON-FATAL error with optim(): \n  ",
                     MLEs$convergence, "  ", MLEs$message, "\n", sep = "") )
        }
        else{
          cat( paste("  There was an error with optim(): \n  ",
                     MLEs$convergence, "  ", MLEs$message, "\n", sep = "") )
        }
      }

      # Save the kernel matrix
      mc.kernels[,,k] <- kernel_cov( covMLE[1:3] )

      # Calculate spatially-varying mean coefficients if needed
      if( ns.mean ){
        if( method == "ml" ){
          beta.GLS.save[k,] <- betaMLE
          beta.cov.save[,,k] <- matrix(NA, q, q)
          beta.coefs.save[[k]] <- NA
        }
        if( method == "reml" ){
          dist.k <- mahalanobis.dist( data.x = temp.locations, vc = mc.kernels[,,k] )
          NS.cov.k <- covMLE[5]*cov.spatial(dist.k, cov.model = cov.model,
                                            cov.pars = c(1,1), kappa = covMLE[6])
          Data.cov.k <- NS.cov.k + diag(rep(covMLE[4],n.fit)) + temp.nugg2.var
          Data.chol.k <- chol(Data.cov.k)
          tX.Cinv.k <- t(backsolve(Data.chol.k, backsolve(Data.chol.k, Xtemp, transpose = TRUE)))
          beta.cov.k <- chol2inv( chol( tX.Cinv.k%*%Xtemp) )/p
          beta.GLS.k <- (p*beta.cov.k %*% tX.Cinv.k %*% temp.data)/p

          beta.GLS.save[k,] <- beta.GLS.k
          beta.cov.save[,,k] <- beta.cov.k
          beta.coefs.save[[k]] <- data.frame( Estimate = beta.GLS.k,
                                              Std.Error = sqrt(diag(beta.cov.k)),
                                              t.val = beta.GLS.k/sqrt(diag(beta.cov.k)) )
        }
      }

      # Save all MLEs
      # Parameter order: n, lam1, lam2, eta, tausq, sigmasq, kappa
      MLEs.save[k,] <- c( n.fit, covMLE )
    }

    # Put the MLEs into a data frame
    MLEs.save <- data.frame( MLEs.save )
    names(MLEs.save) <- c("n","lam1", "lam2", "eta", "tausq", "sigmasq", "kappa" )
  }

  #===========================================================================
  # Calculate the weights for each observation location
  #===========================================================================
  weights <- matrix(NA, N, K)
  for(n in 1:N){for(k in 1:K){
    weights[n,k] <- exp(-sum((coords[n,] - mc.locations[k,])^2)/(2*lambda.w))
  }
    # Normalize the weights
    weights[n,] <- weights[n,]/sum(weights[n,])
  }

  #===========================================================================
  # Calculate the kernel ellipses and other spatially-varying quantities
  #===========================================================================
  kernel.ellipses <- array(0, dim=c(2,2,N))
  for(n in 1:N){
    for(k in 1:K){
      kernel.ellipses[,,n] <- kernel.ellipses[,,n] + weights[n,k]*mc.kernels[,,k]
    }
  }

  # If specified: calculate the spatially-varying nugget and variance
  if( ns.nugget == TRUE ){
    mc.nuggets <- as.numeric(MLEs.save$tausq)
    obs.nuggets <- rep(0,N)
    for(n in 1:N){
      for(k in 1:K){
        obs.nuggets[n] <- obs.nuggets[n] + weights[n,k]*mc.nuggets[k]
      }
    }
    # obs.nuggets <- obs.nuggets + fixed.nugg2.var
  }

  if( ns.variance == TRUE ){
    mc.variance <- as.numeric(MLEs.save$sigmasq)
    obs.variance <- rep(0,N)
    for(n in 1:N){
      for(k in 1:K){
        obs.variance[n] <- obs.variance[n] + weights[n,k]*mc.variance[k]
      }
    }
  }

  if( ns.mean == TRUE ){
    obs.beta <- matrix(0, N, ncol(Xmat))
    for( t in 1:ncol(Xmat)){
      for(n in 1:N){
        for(k in 1:K){
          obs.beta[n,t] <- obs.beta[n,t] + weights[n,k]*beta.GLS.save[k,t]
        }
      }
    }
  }
  cat("-----------------------------------------------------------\n")

  #===========================================================================
  # Global parameter estimation
  #===========================================================================
  # First, for models without kappa.
  if( cov.model != "matern" & cov.model != "cauchy" ){

    cat("Calculating the nonstationary correlation matrix.\n")
    KAPPA <- NULL
    #====================================
    # Calculate the correlation matrix

    Scale.mat <- matrix(rep(NA, N^2), nrow=N)
    Dist.mat <- matrix(rep(NA, N^2), nrow=N)
    # Calculate the elements of the observed correlations matrix.
    for(i in 1:N){

      # Diagonal elements
      Kerneli <- kernel.ellipses[,,i]
      det_i <- Kerneli[1,1]*Kerneli[2,2] - Kerneli[1,2]*Kerneli[2,1]

      Scale.mat[i,i] <- 1
      Dist.mat[i,i] <- 0

      # Ui <- chol(Kerneli)

      if(i < N){
        for(j in (i+1):N){ # Off-diagonal elements

          Kernelj <- kernel.ellipses[,,j]
          det_j <- Kernelj[1,1]*Kernelj[2,2] - Kernelj[1,2]*Kernelj[2,1]

          avg_ij <- 0.5 * (Kerneli + Kernelj)
          det_ij <- avg_ij[1,1]*avg_ij[2,2] - avg_ij[1,2]*avg_ij[2,1]

          Scale.mat[i,j] <- sqrt( sqrt(det_i*det_j) / det_ij )
          Dist.mat[i,j] <- sqrt( t(coords[i,]-coords[j,]) %*% solve(avg_ij) %*% (coords[i,]-coords[j,]) )

          Scale.mat[j,i] <- Scale.mat[i,j]
          Dist.mat[j,i] <- Dist.mat[i,j]

        }
      }
    }
    Unscl.corr <- cov.spatial( Dist.mat, cov.model = cov.model,
                               cov.pars = c(1,1), kappa = KAPPA )
    NS.corr <- Scale.mat*Unscl.corr

    #====================================
    # Global parameter estimation
    if( ns.nugget == FALSE & ns.variance == FALSE ){

      if(print.progress){
        cat("Calculating the variance parameter MLEs. \n")
      }

      overall.lik1 <- make_global_loglik1( data = data, Xmat = Xmat,
                                           Corr = NS.corr,
                                           nugg2.var = fixed.nugg2.var )

      overall.MLEs <- optim( c( tausq.global.init, sigmasq.global.init ),
                             overall.lik1,
                             method = "L-BFGS-B",
                             lower=c( tausq.global.LB, sigmasq.global.LB ),
                             upper=c( tausq.global.UB, sigmasq.global.UB ) )

      if(print.progress){
        if( overall.MLEs$convergence != 0 ){
          if( overall.MLEs$convergence == 52 ){
            cat( paste("  There was a NON-FATAL error with optim(): \n  ",
                       overall.MLEs$convergence, "  ", overall.MLEs$message, "\n", sep = "") )
          }
          else{
            cat( paste("  There was an error with optim(): \n  ",
                       overall.MLEs$convergence, "  ", overall.MLEs$message, "\n", sep = "") )
          }
        }
      }

      tausq.MLE <- overall.MLEs$par[1]
      sigmasq.MLE <- overall.MLEs$par[2]

      global.lik <- overall.MLEs$value

      ObsNuggMat <- diag(rep(tausq.MLE,N)) + fixed.nugg2.var
      ObsCov <- sigmasq.MLE * NS.corr

      obs.variance <- rep(sigmasq.MLE, N)
    }

    if( ns.nugget == TRUE & ns.variance == FALSE ){
      if(print.progress){
        cat("Calculating the variance parameter MLEs. \n")
      }

      overall.lik2 <- make_global_loglik2( data = data,
                                           Xmat = Xmat,
                                           Corr = NS.corr,
                                           obs.nuggets = obs.nuggets,
                                           nugg2.var = fixed.nugg2.var )

      overall.MLEs <- optim( sigmasq.global.init, overall.lik2, method = "L-BFGS-B",
                             lower=c( sigmasq.global.LB ),
                             upper=c( sigmasq.global.UB ) )
      if(print.progress){
        if( overall.MLEs$convergence != 0 ){
          if( overall.MLEs$convergence == 52 ){
            cat( paste("  There was a NON-FATAL error with optim(): \n  ",
                       overall.MLEs$convergence, "  ", overall.MLEs$message, "\n", sep = "") )
          }
          else{
            cat( paste("  There was an error with optim(): \n  ",
                       overall.MLEs$convergence, "  ", overall.MLEs$message, "\n", sep = "") )
          }
        }
      }
      sigmasq.MLE <- overall.MLEs$par[1]
      global.lik <- overall.MLEs$value

      ObsNuggMat <- diag(obs.nuggets) + fixed.nugg2.var
      ObsCov <- sigmasq.MLE * NS.corr

      obs.variance <- rep(sigmasq.MLE, N)

    }

    if( ns.nugget == FALSE & ns.variance == TRUE ){
      if(print.progress){
        cat("Calculating the variance parameter MLEs. \n")
      }

      overall.lik3 <- make_global_loglik3( data = data,
                                           Xmat = Xmat,
                                           Corr = NS.corr,
                                           obs.variance = obs.variance,
                                           nugg2.var = fixed.nugg2.var )

      overall.MLEs <- optim( tausq.global.init, overall.lik3, method = "L-BFGS-B",
                             lower=c( tausq.global.LB ),
                             upper=c( tausq.global.UB ) )

      if(print.progress){
        if( overall.MLEs$convergence != 0 ){
          if( overall.MLEs$convergence == 52 ){
            cat( paste("  There was a NON-FATAL error with optim(): \n  ",
                       overall.MLEs$convergence, "  ", overall.MLEs$message, "\n", sep = "") )
          }
          else{
            cat( paste("  There was an error with optim(): \n  ",
                       overall.MLEs$convergence, "  ", overall.MLEs$message, "\n", sep = "") )
          }
        }
      }
      tausq.MLE <- overall.MLEs$par[1]
      global.lik <- overall.MLEs$value

      ObsNuggMat <- diag(rep(tausq.MLE,N)) + fixed.nugg2.var
      ObsCov <- diag(sqrt(obs.variance)) %*% NS.corr %*% diag(sqrt(obs.variance))

    }

    if( ns.nugget == TRUE & ns.variance == TRUE ){

      Cov <- diag( sqrt(obs.variance) ) %*% NS.corr %*% diag( sqrt(obs.variance) )

      global.lik <- NA
      ObsNuggMat <- diag(obs.nuggets) + fixed.nugg2.var
      ObsCov <- Cov

    }
    kappa.MLE <- NA
  }

  # Next, for models with kappa.
  if( cov.model == "matern" || cov.model == "cauchy" ){
    cat("Calculating the nonstationary correlation matrix.\n")

    #====================================
    # Calculate the correlation matrix

    Scale.mat <- matrix(rep(NA, N^2), nrow=N)
    Dist.mat <- matrix(rep(NA, N^2), nrow=N)
    # Calculate the elements of the observed correlations matrix.
    for(i in 1:N){

      # Diagonal elements
      Kerneli <- kernel.ellipses[,,i]
      det_i <- Kerneli[1,1]*Kerneli[2,2] - Kerneli[1,2]*Kerneli[2,1]

      Scale.mat[i,i] <- 1
      Dist.mat[i,i] <- 0

      #Ui <- chol(Kerneli)

      if(i < N){
        for(j in (i+1):N){ # Off-diagonal elements

          Kernelj <- kernel.ellipses[,,j]
          det_j <- Kernelj[1,1]*Kernelj[2,2] - Kernelj[1,2]*Kernelj[2,1]

          avg_ij <- 0.5 * (Kerneli + Kernelj)
          det_ij <- avg_ij[1,1]*avg_ij[2,2] - avg_ij[1,2]*avg_ij[2,1]

          Scale.mat[i,j] <- sqrt( sqrt(det_i*det_j) / det_ij )
          Dist.mat[i,j] <- sqrt( t(coords[i,]-coords[j,]) %*% solve(avg_ij) %*% (coords[i,]-coords[j,]) )

          Scale.mat[j,i] <- Scale.mat[i,j]
          Dist.mat[j,i] <- Dist.mat[i,j]

        }
      }
    }

    #====================================
    # Global parameter estimation
    if( !fix.kappa ){
      if( ns.nugget == FALSE & ns.variance == FALSE ){
        if(print.progress){
          cat("Calculating the variance and smoothness parameter MLEs. \n")
        }

        overall.lik1.kappa <- make_global_loglik1_kappa( data = data, Xmat = Xmat,
                                                         cov.model = cov.model,
                                                         Scalemat = Scale.mat,
                                                         Distmat = Dist.mat,
                                                         nugg2.var = fixed.nugg2.var )

        overall.MLEs <- optim(c( tausq.global.init, sigmasq.global.init, kappa.global.init ),
                              overall.lik1.kappa,
                              method = "L-BFGS-B",
                              lower=c( tausq.global.LB, sigmasq.global.LB, kappa.global.LB ),
                              upper=c( tausq.global.UB, sigmasq.global.UB, kappa.global.UB ) )
        if(print.progress){
          if( overall.MLEs$convergence != 0 ){
            if( overall.MLEs$convergence == 52 ){
              cat( paste("  There was a NON-FATAL error with optim(): \n  ",
                         overall.MLEs$convergence, "  ", overall.MLEs$message, "\n", sep = "") )
            }
            else{
              cat( paste("  There was an error with optim(): \n  ",
                         overall.MLEs$convergence, "  ", overall.MLEs$message, "\n", sep = "") )
            }
          }
        }
        tausq.MLE <- overall.MLEs$par[1]
        sigmasq.MLE <- overall.MLEs$par[2]
        kappa.MLE <- overall.MLEs$par[3]

        global.lik <- overall.MLEs$value
        ObsNuggMat <- diag(rep(tausq.MLE,N)) + fixed.nugg2.var
        ObsCov <- sigmasq.MLE * Scale.mat * cov.spatial( Dist.mat, cov.model = cov.model,
                                                         cov.pars = c(1,1), kappa = kappa.MLE )
        obs.variance <- rep(sigmasq.MLE, N)
      }

      if( ns.nugget == TRUE & ns.variance == FALSE ){
        if(print.progress){
          cat("Calculating the variance and smoothness parameter MLEs. \n")
        }

        overall.lik2.kappa <- make_global_loglik2_kappa( data = data, Xmat = Xmat,
                                                         cov.model = cov.model,
                                                         Scalemat = Scale.mat,
                                                         Distmat = Dist.mat,
                                                         obs.nuggets = obs.nuggets,
                                                         nugg2.var = fixed.nugg2.var )

        overall.MLEs <- optim(c( sigmasq.global.init, kappa.global.init ),
                              overall.lik2.kappa,
                              method = "L-BFGS-B",
                              lower=c( sigmasq.global.LB, kappa.global.LB ),
                              upper=c( sigmasq.global.UB, kappa.global.UB ) )
        if(print.progress){
          if( overall.MLEs$convergence != 0 ){
            if( overall.MLEs$convergence == 52 ){
              cat( paste("  There was a NON-FATAL error with optim(): \n  ",
                         overall.MLEs$convergence, "  ", overall.MLEs$message, "\n", sep = "") )
            }
            else{
              cat( paste("  There was an error with optim(): \n  ",
                         overall.MLEs$convergence, "  ", overall.MLEs$message, "\n", sep = "") )
            }
          }
        }
        sigmasq.MLE <- overall.MLEs$par[1]
        kappa.MLE <- overall.MLEs$par[2]

        global.lik <- overall.MLEs$value
        ObsNuggMat <- diag(obs.nuggets) + fixed.nugg2.var
        ObsCov <- sigmasq.MLE * Scale.mat * cov.spatial( Dist.mat, cov.model = cov.model,
                                                         cov.pars = c(1,1), kappa = kappa.MLE )

        obs.variance <- rep(sigmasq.MLE, N)

      }

      if( ns.nugget == FALSE & ns.variance == TRUE ){
        if(print.progress){
          cat("Calculating the variance and smoothness parameter MLEs. \n")
        }

        overall.lik3.kappa <- make_global_loglik3_kappa( data = data, Xmat = Xmat,
                                                         cov.model = cov.model,
                                                         Scalemat = Scale.mat,
                                                         Distmat = Dist.mat,
                                                         obs.variance = obs.variance,
                                                         nugg2.var = fixed.nugg2.var )

        overall.MLEs <- optim(c( tausq.global.init, kappa.global.init ),
                              overall.lik3.kappa,
                              method = "L-BFGS-B",
                              lower=c( tausq.global.LB, kappa.global.LB ),
                              upper=c( tausq.global.UB, kappa.global.UB ) )
        if(print.progress){
          if( overall.MLEs$convergence != 0 ){
            if( overall.MLEs$convergence == 52 ){
              cat( paste("  There was a NON-FATAL error with optim(): \n  ",
                         overall.MLEs$convergence, "  ", overall.MLEs$message, "\n", sep = "") )
            }
            else{
              cat( paste("  There was an error with optim(): \n  ",
                         overall.MLEs$convergence, "  ", overall.MLEs$message, "\n", sep = "") )
            }
          }
        }
        tausq.MLE <- overall.MLEs$par[1]
        kappa.MLE <- overall.MLEs$par[2]

        global.lik <- overall.MLEs$value
        ObsNuggMat <- diag(rep(tausq.MLE,N)) + fixed.nugg2.var
        ObsCov <- diag(sqrt(obs.variance)) %*% (Scale.mat * cov.spatial( Dist.mat, cov.model = cov.model,
                                                                         cov.pars = c(1,1), kappa = kappa.MLE )) %*% diag(sqrt(obs.variance))

      }

      if( ns.nugget == TRUE & ns.variance == TRUE ){
        if(print.progress){
          cat("Calculating the smoothness parameter MLE. \n")
        }

        overall.lik4.kappa <- make_global_loglik4_kappa( data = data, Xmat = Xmat,
                                                         cov.model = cov.model,
                                                         Scalemat = Scale.mat,
                                                         Distmat = Dist.mat,
                                                         obs.nuggets = obs.nuggets,
                                                         obs.variance = obs.variance,
                                                         nugg2.var = fixed.nugg2.var )

        overall.MLEs <- optim( kappa.global.init, overall.lik4.kappa, method = "L-BFGS-B",
                               lower=c( kappa.global.LB ),
                               upper=c( kappa.global.UB ) )
        if(print.progress){
          if( overall.MLEs$convergence != 0 ){
            if( overall.MLEs$convergence == 52 ){
              cat( paste("  There was a NON-FATAL error with optim(): \n  ",
                         overall.MLEs$convergence, "  ", overall.MLEs$message, "\n", sep = "") )
            }
            else{
              cat( paste("  There was an error with optim(): \n  ",
                         overall.MLEs$convergence, "  ", overall.MLEs$message, "\n", sep = "") )
            }
          }
        }

        kappa.MLE <- overall.MLEs$par[1]

        global.lik <- overall.MLEs$value

        ObsNuggMat <- diag(obs.nuggets) + fixed.nugg2.var
        ObsCov <- diag(sqrt(obs.variance)) %*% (Scale.mat * cov.spatial( Dist.mat, cov.model = cov.model,
                                                                         cov.pars = c(1,1), kappa = kappa.MLE )) %*% diag(sqrt(obs.variance))

      }
    } else{

      kappa.MLE <- kappa

      global.lik <- NA

      ObsNuggMat <- diag(obs.nuggets) + fixed.nugg2.var
      ObsCov <- diag(sqrt(obs.variance)) %*% (Scale.mat * cov.spatial( Dist.mat, cov.model = cov.model,
                                                                       cov.pars = c(1,1), kappa = kappa.MLE )) %*% diag(sqrt(obs.variance))


    }


  }

  Data.Cov <- ObsNuggMat + ObsCov
  #===========================================================================
  # Calculate the GLS estimate of beta
  #===========================================================================
  Data.Cov.chol <- chol(Data.Cov)
  if( ns.mean == FALSE ){
    tX.Cinv <- t(backsolve(Data.Cov.chol, backsolve(Data.Cov.chol, Xmat, transpose = TRUE)))
    beta.cov <- chol2inv( chol( tX.Cinv%*%Xmat) )/p
    beta.GLS <- (p*beta.cov %*% tX.Cinv %*% data)/p
    Mean.coefs <- data.frame( Estimate = beta.GLS,
                              Std.Error = sqrt(diag(beta.cov)),
                              t.val = beta.GLS/sqrt(diag(beta.cov)) )
  }
  if( ns.mean == TRUE ){
    beta.cov <- beta.cov.save
    beta.GLS <- beta.GLS.save
    Mean.coefs <- beta.coefs.save
  }

  #===========================================================================
  # Output
  #===========================================================================
  if( ns.nugget == TRUE ){
    tausq.out <- obs.nuggets # This is the variance of nugget1 only
  }
  if( ns.nugget == FALSE ){
    tausq.out <- tausq.MLE
  }

  if( ns.variance == TRUE ){
    sigmasq.out <- obs.variance
  }
  if( ns.variance == FALSE ){
    sigmasq.out <- sigmasq.MLE
  }

  if( ns.mean == TRUE ){
    beta.out <- obs.beta
    names(beta.out) <- beta.names
  }
  if( ns.mean == FALSE ){
    beta.out <- beta.GLS
  }

  output <- list( mc.kernels = mc.kernels,
                  mc.locations = mc.locations,
                  MLEs.save = MLEs.save,
                  kernel.ellipses = kernel.ellipses,
                  data = data,
                  beta.GLS = beta.GLS,
                  beta.cov = beta.cov,
                  Mean.coefs = Mean.coefs,
                  tausq.est = tausq.out,
                  sigmasq.est = sigmasq.out,
                  beta.est = beta.out,
                  kappa.MLE = kappa.MLE,
                  Cov.mat = Data.Cov,
                  Cov.mat.chol = Data.Cov.chol,
                  cov.model = cov.model,
                  ns.nugget = ns.nugget,
                  ns.variance = ns.variance,
                  ns.mean = ns.mean,
                  fixed.nugg2.var = fixed.nugg2.var,
                  coords = coords,
                  global.loglik = global.lik,
                  Xmat = Xmat,
                  lambda.w = lambda.w,
                  fix.kappa = fix.kappa, kappa = kappa )

  class(output) <- "NSconvo"
  cat("Done.")
  cat("\n-----------------------------------------------------------\n")

  return(output)

}

#======================================================================================
# Calculate predictions using the output of NSconvo_fit()
#======================================================================================
# Using the output from NSconvo_fit(), calculate the kriging predictors
# and kriging standard errors for prediction locations of interest.
#======================================================================================
#ROxygen comments ----
#' Obtain predictions at unobserved locations for the nonstationary
#' spatial model.
#'
#' \code{predict.NSconvo} calculates the kriging predictor and corresponding
#' standard errors at unmonitored sites.
#'
#' @param object A "NSconvo" object, from \code{NSconvo_fit}.
#' @param pred.coords Matrix of locations where predictions are required.
#' @param pred.covariates Matrix of covariates for the prediction locations,
#' NOT including an intercept. The number of columns for this matrix must
#' match the design matrix from \code{mean.model} in \code{\link{NSconvo_fit}}.
#' Defaults to an intercept only.
#' @param pred.fixed.nugg2.var An optional vector or matrix describing the
#' the variance/covariance a fixed second nugget term (corresponds to
#' \code{fixed.nugg2.var} in \code{NSconvo_fit}; often useful if conducting
#' prediction for held-out data). Defaults to zero.
#' @param ... additional arguments affecting the predictions produced.
#'
#' @return A list with the following components:
#' \item{pred.means}{Vector of the kriging predictor, for each location in
#' \code{pred.coords}.}
#' \item{pred.SDs}{Vector of the kriging standard errors, for each location
#' in \code{pred.coords}.}
#'
#' @examples
#' \dontrun{
#' pred.NS <- predict( NSconvo.obj,
#' pred.coords = matrix(c(1,1), ncol=2),
#' pred.covariates = matrix(c(1,1), ncol=2) )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial
#'

predict.NSconvo <- function(object, pred.coords, pred.covariates = NULL,
                            pred.fixed.nugg2.var = NULL, ... )
{
  if( !inherits(object, "NSconvo") ){
    warning("Object is not of type NSconvo.")
  }
  else{
    {
      pred.coords <- as.matrix(pred.coords)
      M <- dim(pred.coords)[1]
      mc.locations <- object$mc.locations
      K <- dim(mc.locations)[1]
      mc.kernels <- object$mc.kernels
      kernel.ellipses <- object$kernel.ellipses
      ns.nugget <- object$ns.nugget
      ns.variance <- object$ns.variance
      ns.mean <- object$ns.mean
      beta.MLE <- as.matrix(object$beta.GLS)
      beta.est <- object$beta.est
      tausq.est <- object$tausq.est
      sigmasq.est <- object$sigmasq.est
      if( object$fix.kappa ){
        kappa.MLE <- object$kappa
      } else{
        kappa.MLE <- object$kappa.MLE
      }
      mc.MLEs <- object$MLEs.save
      Cov.mat.chol <- object$Cov.mat.chol
      data <- object$data
      N <- length(object$data)
      coords <- object$coords
      cov.model <- object$cov.model
      Xmat <- object$Xmat
      lambda.w <- object$lambda.w
      fixed.nugg2.var <- object$fixed.nugg2.var

      if (is.null(pred.covariates) == TRUE) {
        Xpred <- rep(1, M)
      }
      if (is.null(pred.covariates) == FALSE) {
        Xpred <- cbind(rep(1, M), pred.covariates)
      }

      data <- matrix(data, nrow = N)
      p <- dim(data)[2]
      pred.weights <- matrix(NA, M, K)
      for (m in 1:M) {
        for (k in 1:K) {
          pred.weights[m, k] <- exp(-sum((pred.coords[m, ] -
                                            mc.locations[k, ])^2)/(2 * lambda.w))
        }
        pred.weights[m, ] <- pred.weights[m, ]/sum(pred.weights[m,])
      }
      pred.kernel.ellipses <- array(0, dim = c(2, 2, M))
      for (m in 1:M) {
        for (k in 1:K) {
          pred.kernel.ellipses[, , m] <- pred.kernel.ellipses[,, m] + pred.weights[m, k] * mc.kernels[, , k]
        }
      }
      if (ns.nugget == TRUE) {
        mc.nuggets <- mc.MLEs$tausq
        pred.nuggets <- rep(0, M)
        for (m in 1:M) {
          for (k in 1:K) {
            pred.nuggets[m] <- pred.nuggets[m] + pred.weights[m,
                                                              k] * mc.nuggets[k]
          }
        }
      }
      if (ns.nugget == FALSE) {
        pred.nuggets <- rep(tausq.est, M)
      }
      if (ns.variance == TRUE) {
        obs.variance <- sigmasq.est
        mc.variance <- mc.MLEs$sigmasq
        pred.variance <- rep(0, M)
        for (m in 1:M) {
          for (k in 1:K) {
            pred.variance[m] <- pred.variance[m] + pred.weights[m,
                                                                k] * mc.variance[k]
          }
        }
      }
      if (ns.variance == FALSE) {
        obs.variance <- rep(sigmasq.est, N)
        pred.variance <- rep(sigmasq.est, M)
      }

      if( is.null(pred.fixed.nugg2.var) == TRUE ){
        pred.fixed.nugg2.var <- matrix(0,M,M)
      } else{
        if( !is.matrix(pred.fixed.nugg2.var) ){ # Convert to a covariance matrix if specified as a vector
          pred.fixed.nugg2.var <- diag(pred.fixed.nugg2.var)
        }
      }

      if( ns.mean == TRUE ){
        pred.beta <- matrix(0, M, ncol(Xmat))
        for( t in 1:ncol(Xmat)){
          for(m in 1:M){
            for(k in 1:K){
              pred.beta[m,t] <- pred.beta[m,t] + pred.weights[m,k]*beta.MLE[k,t]
            }
          }
        }
      }

      cat("-----------------------------------------------------------\n")
      cat(paste("Calculating the ", M, " by ", N, " cross-correlation matrix.\n", sep=""))
      if( max(c(M,N)) > 1000 ){
        cat("NOTE: this step can be VERY time consuming.\n")
      }
      Scale.cross <- matrix(NA, M, N)
      Dist.cross <- matrix(NA, M, N)
      cat("Progress: \n")
      cat("|--------|---------|---------|---------|")
      cat("\n|")
      for (i in 1:N) {
        Kerneli <- kernel.ellipses[, , i]
        det_i <- Kerneli[1, 1] * Kerneli[2, 2] - Kerneli[1, 2] *
          Kerneli[2, 1]
        for (j in 1:M) {
          Kernelj <- pred.kernel.ellipses[, , j]
          det_j <- Kernelj[1, 1] * Kernelj[2, 2] - Kernelj[1,
                                                           2] * Kernelj[2, 1]
          avg_ij <- 0.5 * (Kerneli + Kernelj)
          Uij <- chol(avg_ij)
          det_ij <- avg_ij[1, 1] * avg_ij[2, 2] - avg_ij[1,
                                                         2] * avg_ij[2, 1]
          vec_ij <- backsolve(Uij, (coords[i, ] - pred.coords[j,
                                                              ]), transpose = TRUE)
          Scale.cross[j, i] <- sqrt(sqrt(det_i * det_j)/det_ij)
          Dist.cross[j, i] <- sqrt(sum(vec_ij^2))
        }
        if (i%%floor(N/38) == 0) {
          cat("-")
        }
      }
      cat("|\n")
      Unscl.cross <- cov.spatial(Dist.cross, cov.model = cov.model,
                                 cov.pars = c(1, 1), kappa = kappa.MLE)
      NS.cross.corr <- Scale.cross * Unscl.cross
      CrossCov <- diag(sqrt(pred.variance)) %*% NS.cross.corr %*%
        diag(sqrt(obs.variance))

      crscov.Cinv <- t(backsolve(Cov.mat.chol,
                                 backsolve(Cov.mat.chol, t(CrossCov), transpose = TRUE)))

      cat("Calculating kriging means and standard errors...\n")
      pred.means <- matrix(NA, M, p)
      if( ns.mean == FALSE ){
        for (i in 1:p) {
          pred.means[, i] <- Xpred %*% beta.MLE + crscov.Cinv %*% (object$data[, i] - Xmat %*% beta.MLE)
        }
      }
      if( ns.mean == TRUE ){
        pred.Xbeta.hat <- rep(0,M)
        obs.Xbeta.hat <- rep(0,N)
        for( t in 1:ncol(Xmat) ){
          pred.Xbeta.hat <- pred.Xbeta.hat + Xpred[,t]*pred.beta[,t]
          obs.Xbeta.hat <- obs.Xbeta.hat + Xmat[,t]*beta.est[,t]
        }
        for (i in 1:p) {
          pred.means[, i] <- pred.Xbeta.hat + crscov.Cinv %*% (object$data[, i] - obs.Xbeta.hat)
        }
      }

      pred.SDs <- sqrt((pred.variance + pred.nuggets + diag(pred.fixed.nugg2.var)) - diag(crscov.Cinv %*% t(CrossCov)))
      output <- list(pred.means = pred.means, pred.SDs = pred.SDs)
      cat("Done.")
      cat("\n-----------------------------------------------------------\n")
      return(output)
    }
  }
}

#======================================================================================
# Fit the anisotropic model
#======================================================================================
# The Aniso_fit() function estimates the parameters of the anisotropic spatial
# model. Required inputs are the observed data and locations (a geoR object with
# $coords and $data). Optional inputs include the covariance model (exponential is
# the default) and the mean model.
#======================================================================================
#ROxygen comments ----
#' Fit the stationary spatial model
#'
#' \code{Aniso_fit} estimates the parameters of the stationary spatial model.
#' Required inputs are the observed data and locations (a geoR object
#' with $coords and $data). Optional inputs include the covariance model
#' (exponential is the default).
#'
#' @param geodata A list containing elements \code{coords} and \code{data} as
#' described next. Typically an object of the class "\code{geodata}", although
#' a geodata object only allows \code{data} to be a vector (no replicates).
#' If not provided, the arguments \code{coords} and \code{data} must be
#' provided instead.
#' @param sp.SPDF A "\code{SpatialPointsDataFrame}" object, which contains the
#' spatial coordinates and additional attribute variables corresponding to the
#' spatoal coordinates
#' @param coords An N x 2 matrix where each row has the two-dimensional
#' coordinates of the N data locations. By default, it takes the \code{coords}
#' component of the argument \code{geodata}, if provided.
#' @param data A vector or matrix with N rows, containing the data values.
#' Inputting a vector corresponds to a single replicate of data, while
#' inputting a matrix corresponds to replicates. In the case of replicates,
#' the model assumes the replicates are independent and identically
#' distributed.
#' @param cov.model A string specifying the model for the correlation
#' function; following \code{geoR}, defaults to \code{"exponential"}.
#' Options available in this package are: "\code{exponential}",
#' \code{"cauchy"}, \code{"matern"}, \code{"circular"}, \code{"cubic"},
#' \code{"gaussian"}, \code{"spherical"}, and \code{"wave"}. For further
#' details, see documentation for \code{\link[geoR]{cov.spatial}}.
#' @param mean.model An object of class \code{\link[stats]{formula}},
#' specifying the mean model to be used. Defaults to an intercept only.
#' @param fixed.nugg2.var Optional; describes the variance/covariance for
#' a fixed (second) nugget term (represents a known error term). Either
#' a vector of length N containing a station-specific variances (implying
#' independent error) or an NxN covariance matrix (implying dependent error).
#' Defaults to zero.
#' @param fix.tausq Logical; indicates whether the default nugget term
#' (tau^2) should be fixed (\code{TRUE}) or estimated (\code{FALSE}). Defaults to
#' \code{FALSE}.
#' @param tausq Scalar; fixed value for the nugget variance (when
#' \code{fix.tausq = TRUE}).
#'
#' @param fix.kappa Logical; indicates if the kappa parameter should be
#' fixed (\code{TRUE}) or estimated (\code{FALSE}). Defaults to \code{FALSE}
#' (only valid for \code{cov.model = "matern"} and \code{cov.model = "cauchy"}).
#' @param kappa Scalar; value of the kappa parameter. Only used if
#' \code{fix.kappa = TRUE}.
#' @param method Indicates the estimation method, either maximum likelihood
#' (\code{"ml"}) or restricted maximum likelihood (\code{"reml"}).
#'
#' @param local.pars.LB,local.pars.UB Optional vectors of lower and upper
#' bounds, respectively, used by the \code{"L-BFGS-B"} method option in the
#' \code{\link[stats]{optim}} function for the local parameter estimation.
#' Each vector must be of length five,
#' containing values for lam1, lam2, tausq, sigmasq, and nu. Default for
#' \code{local.pars.LB} is \code{rep(1e-05,5)}; default for
#' \code{local.pars.UB} is \code{c(max.distance/2, max.distance/2, 4*resid.var, 4*resid.var, 100)},
#' where \code{max.distance} is the maximum interpoint distance of the
#' observed data and \code{resid.var} is the residual variance from using
#' \code{\link[stats]{lm}} with \code{mean.model}.
#'
#' @param local.ini.pars Optional vector of initial values used by the
#' \code{"L-BFGS-B"} method option in the \code{\link[stats]{optim}}
#' function for the local parameter estimation. The vector must be of length
#' five, containing values for lam1, lam2, tausq, sigmasq, and nu. Defaults
#' to \code{c(max.distance/10, max.distance/10, 0.1*resid.var, 0.9*resid.var, 1)},
#' where \code{max.distance} is the maximum interpoint distance of the
#' observed data and \code{resid.var} is the residual variance from using
#' \code{\link[stats]{lm}} with \code{mean.model}.
#'
#' @return A list with the following components:
#' \item{MLEs.save}{Table of local maximum likelihood estimates for each
#' mixture component location.}
#' \item{data}{Observed data values.}
#' \item{beta.GLS}{Vector of generalized least squares estimates of beta,
#' the mean coefficients.}
#' \item{beta.cov}{Covariance matrix of the generalized least squares
#' estimate of beta.}
#' \item{Mean.coefs}{"Regression table" for the mean coefficient estimates,
#' listing the estimate, standard error, and t-value.}
#' \item{Cov.mat}{Estimated covariance matrix (\code{N.obs} x \code{N.obs})
#' using all relevant parameter estimates.}
#' \item{Cov.mat.chol}{Cholesky of \code{Cov.mat} (i.e., \code{chol(Cov.mat)}),
#' the estimated covariance matrix (\code{N.obs} x \code{N.obs}).}
#' \item{aniso.pars}{Vector of MLEs for the anisotropy parameters lam1,
#' lam2, eta.}
#' \item{aniso.mat}{2 x 2 anisotropy matrix, calculated from
#' \code{aniso.pars}.}
#' \item{tausq.est}{Scalar maximum likelihood estimate of tausq (nugget
#' variance).}
#' \item{sigmasq.est}{Scalar maximum likelihood estimate of sigmasq
#' (process variance).}
#' \item{kappa.MLE}{Scalar maximum likelihood estimate for kappa (when
#' applicable).}
#' \item{fixed.nugg2.var}{N x N matrix with the fixed
#' variance/covariance for the second (measurement error) nugget term (defaults
#' to zero).}
#' \item{cov.model}{String; the correlation model used for estimation.}
#' \item{coords}{N x 2 matrix of observation locations.}
#' \item{global.loglik}{Scalar value of the maximized likelihood from the
#' global optimization (if available).}
#' \item{Xmat}{Design matrix, obtained from using \code{\link[stats]{lm}}
#' with \code{mean.model}.}
#' \item{fix.kappa}{Logical, indicating if kappa was fixed (\code{TRUE}) or
#' estimated (\code{FALSE}).}
#' \item{kappa}{Scalar; fixed value of kappa.}
#'
#' @examples
#' # Using iid standard Gaussian data
#' aniso.fit <- Aniso_fit( coords = cbind(runif(100), runif(100)),
#' data = rnorm(100) )
#'
#'
#' @export
#' @importFrom geoR cov.spatial
#' @importFrom StatMatch mahalanobis.dist
#' @importFrom stats lm
#' @importFrom stats optim

Aniso_fit <- function( geodata = NULL, sp.SPDF = NULL,
                       coords = geodata$coords, data = geodata$data,
                       cov.model = "exponential", mean.model = data ~ 1,
                       fixed.nugg2.var = NULL, method = "reml",
                       fix.tausq = FALSE, tausq = 0,
                       fix.kappa = FALSE, kappa = 0.5,
                       local.pars.LB = NULL, local.pars.UB = NULL,
                       local.ini.pars = NULL ){

  #===========================================================================
  # Formatting for coordinates/data
  #===========================================================================
  if( is.null(geodata) == FALSE ){
    if( class(geodata) != "geodata" ){
      cat("\nPlease use a geodata object for the 'geodata = ' input.\n")
    }
    coords <- geodata$coords
    data <- geodata$data
  }
  if( is.null(sp.SPDF) == FALSE ){
    if( class(sp.SPDF) != "SpatialPointsDataFrame" ){
      cat("\nPlease use a SpatialPointsDataFrame object for the 'sp.SPDF = ' input.\n")
    }
    geodata <- geoR::as.geodata( sp.SPDF )
    coords <- geodata$coords
    data <- geodata$data
  }

  N <- dim(coords)[1]
  data <- as.matrix(data)
  p <- dim(data)[2]

  # Make sure cov.model is one of the permissible options
  if( cov.model != "cauchy" & cov.model != "matern" & cov.model != "circular" &
      cov.model != "cubic" & cov.model != "gaussian" & cov.model != "exponential" &
      cov.model != "spherical" & cov.model != "wave" ){
    cat("Please specify a valid covariance model (cauchy, matern, circular, cubic, gaussian, exponential, spherical, or wave).")
  }

  #===========================================================================
  # Calculate the design matrix
  #===========================================================================
  OLS.model <- lm( mean.model, x=TRUE )

  Xmat <- matrix( unname( OLS.model$x ), nrow=N )
  q <- ncol(Xmat)

  #===========================================================================
  # Set the second nugget variance, if not specified
  #===========================================================================
  if( is.null(fixed.nugg2.var) == TRUE ){
    fixed.nugg2.var <- matrix(0,N,N)
  } else{
    if( !is.matrix(fixed.nugg2.var) ){ # Convert to a covariance matrix if specified as a vector
      fixed.nugg2.var <- diag(fixed.nugg2.var)
    }
  }

  #===========================================================================
  # Specify lower, upper, and initial parameter values for optim()
  #===========================================================================
  lon_min <- min(coords[,1])
  lon_max <- max(coords[,1])
  lat_min <- min(coords[,2])
  lat_max <- max(coords[,2])

  max.distance <- sqrt( sum((c(lon_min,lat_min) - c(lon_max,lat_max))^2))

  if( p > 1 ){
    ols.sigma <- NULL
    for(i in 1:length(names(summary(OLS.model)))){
      ols.sigma <- c( ols.sigma, summary(OLS.model)[[i]]$sigma )
    }
    resid.var <- (max(ols.sigma))^2
  }
  if( p == 1 ){
    resid.var <- summary(OLS.model)$sigma^2
  }

  #=================================
  # Lower limits for optim()
  #=================================
  if( is.null(local.pars.LB) == TRUE ){
    lam1.LB <- 1e-05
    lam2.LB <- 1e-05
    tausq.local.LB <- 1e-05
    sigmasq.local.LB <- 1e-05
    kappa.local.LB <- 1e-05
  }
  if( is.null(local.pars.LB) == FALSE ){
    lam1.LB <- local.pars.LB[1]
    lam2.LB <- local.pars.LB[2]
    tausq.local.LB <- local.pars.LB[3]
    sigmasq.local.LB <- local.pars.LB[4]
    kappa.local.LB <- local.pars.LB[5]
  }

  #=================================
  # Upper limits for optim()
  #=================================
  if( is.null(local.pars.UB) == TRUE ){
    lam1.UB <- max.distance/4
    lam2.UB <- max.distance/4
    tausq.local.UB <- 4*resid.var
    sigmasq.local.UB <- 4*resid.var
    kappa.local.UB <- 30
  }
  if( is.null(local.pars.UB) == FALSE ){
    lam1.UB <- local.pars.UB[1]
    lam2.UB <- local.pars.UB[2]
    tausq.local.UB <- local.pars.UB[3]
    sigmasq.local.UB <- local.pars.UB[4]
    kappa.local.UB <- local.pars.UB[5]
  }

  #=================================
  # Initial values for optim()
  #=================================
  if( is.null(local.ini.pars) == TRUE ){
    lam1.init <- max.distance/10
    lam2.init <- max.distance/10
    tausq.local.init <- 0.1*resid.var
    sigmasq.local.init <- 0.9*resid.var
    kappa.local.init <- 1
  }
  if( is.null(local.ini.pars) == FALSE ){
    lam1.init <- local.ini.pars[1]
    lam2.init <- local.ini.pars[2]
    tausq.local.init <- local.ini.pars[3]
    sigmasq.local.init <- local.ini.pars[4]
    kappa.local.init <- local.ini.pars[5]
  }

  #===========================================================================
  # MLEs
  #===========================================================================

  cat("\n-------------------------------------------------------\n")
  cat("Fitting the stationary (anisotropic) model.")
  cat("\n-------------------------------------------------------\n")

  cat("Estimating the variance/covariance parameters. \n")

  #####################################################
  local.aniso <- TRUE
  if( cov.model == "matern" || cov.model == "cauchy" ){ # For covariance models with kappa
    if( !fix.kappa ){ # Estimate kappa
      if( !fix.tausq ){ # Estimate tausq
        if( method == "ml" ){
          f_loglik <- make_local_lik( locations = coords, cov.model = cov.model,
                                      data = data, Xmat = Xmat, nugg2.var = fixed.nugg2.var,
                                      fixed = rep(FALSE, 6 + q), method = method,
                                      local.aniso = local.aniso, fix.tausq = fix.tausq,
                                      fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
          MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, tausq.local.init,
                                 sigmasq.local.init, kappa.local.init, rep(0, q) ),
                         fn = f_loglik, method = "L-BFGS-B",
                         lower = c( lam1.LB, lam2.LB, 0, tausq.local.LB, sigmasq.local.LB,
                                    kappa.local.LB, rep(-Inf, q) ),
                         upper = c( lam1.UB, lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB,
                                    kappa.local.UB, rep(Inf, q) ) )
          covMLE <- MLEs$par[1:6]
          betaMLE <- MLEs$par[-c(1:6)]
        }
        if( method == "reml" ){
          f_loglik <- make_local_lik( locations = coords, cov.model = cov.model,
                                      data = data, Xmat = Xmat, nugg2.var = fixed.nugg2.var,
                                      fixed = rep(FALSE, 6), method = method,
                                      local.aniso = local.aniso, fix.tausq = fix.tausq,
                                      fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
          MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, tausq.local.init,
                                 sigmasq.local.init, kappa.local.init ),
                         fn = f_loglik, method = "L-BFGS-B",
                         lower = c( lam1.LB, lam2.LB, 0, tausq.local.LB, sigmasq.local.LB,
                                    kappa.local.LB ),
                         upper = c( lam1.UB, lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB,
                                    kappa.local.UB ) )
          covMLE <- MLEs$par
        }
      }
      if( fix.tausq ){ # Fix tausq
        if( method == "ml" ){
          f_loglik <- make_local_lik( locations = coords, cov.model = cov.model,
                                      data = data, Xmat = Xmat, nugg2.var = fixed.nugg2.var,
                                      fixed = rep(FALSE, 5 + q), method = method,
                                      local.aniso = local.aniso, fix.tausq = fix.tausq,
                                      fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
          MLEs <- optim( par = c(lam1.init, lam2.init, pi/4,
                                 sigmasq.local.init, kappa.local.init, rep(0, q) ),
                         fn = f_loglik, method = "L-BFGS-B",
                         lower = c( lam1.LB, lam2.LB, 0, sigmasq.local.LB,
                                    kappa.local.LB, rep(-Inf, q) ),
                         upper = c( lam1.UB, lam2.UB, pi/2, sigmasq.local.UB,
                                    kappa.local.UB, rep(Inf, q) ) )
          covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4:5])
          betaMLE <- MLEs$par[-c(1:5)]
        }
        if( method == "reml" ){
          f_loglik <- make_local_lik( locations = coords, cov.model = cov.model,
                                      data = data, Xmat = Xmat, nugg2.var = fixed.nugg2.var,
                                      fixed = rep(FALSE, 5), method = method,
                                      local.aniso = local.aniso, fix.tausq = fix.tausq,
                                      fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
          MLEs <- optim( par = c(lam1.init, lam2.init, pi/4,
                                 sigmasq.local.init, kappa.local.init ),
                         fn = f_loglik, method = "L-BFGS-B",
                         lower = c( lam1.LB, lam2.LB, 0, sigmasq.local.LB,
                                    kappa.local.LB ),
                         upper = c( lam1.UB, lam2.UB, pi/2, sigmasq.local.UB,
                                    kappa.local.UB ) )
          covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4:5])
        }
      }
    }
    if( fix.kappa ){ # Fix kappa
      if( !fix.tausq ){ # Estimate tausq
        if( method == "ml" ){
          f_loglik <- make_local_lik( locations = coords, cov.model = cov.model,
                                      data = data, Xmat = Xmat, nugg2.var = fixed.nugg2.var,
                                      fixed = rep(FALSE, 5 + q), method = method,
                                      local.aniso = local.aniso, fix.tausq = fix.tausq,
                                      fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
          MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, tausq.local.init,
                                 sigmasq.local.init, rep(0, q) ),
                         fn = f_loglik, method = "L-BFGS-B",
                         lower = c( lam1.LB, lam2.LB, 0, tausq.local.LB, sigmasq.local.LB,
                                    rep(-Inf, q) ),
                         upper = c( lam1.UB, lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB,
                                    rep(Inf, q) ) )
          covMLE <- c(MLEs$par[1:5], kappa)
          betaMLE <- MLEs$par[-c(1:5)]
        }
        if( method == "reml" ){
          f_loglik <- make_local_lik( locations = coords, cov.model = cov.model,
                                      data = data, Xmat = Xmat, nugg2.var = fixed.nugg2.var,
                                      fixed = rep(FALSE, 5), method = method,
                                      local.aniso = local.aniso, fix.tausq = fix.tausq,
                                      fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
          MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, tausq.local.init,
                                 sigmasq.local.init ),
                         fn = f_loglik, method = "L-BFGS-B",
                         lower = c( lam1.LB, lam2.LB, 0, tausq.local.LB, sigmasq.local.LB ),
                         upper = c( lam1.UB, lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB ) )
          covMLE <- c(MLEs$par[1:5], kappa)
        }
      }
      if( fix.tausq ){ # Fix tausq
        if( method == "ml" ){
          f_loglik <- make_local_lik( locations = coords, cov.model = cov.model,
                                      data = data, Xmat = Xmat, nugg2.var = fixed.nugg2.var,
                                      fixed = rep(FALSE, 4 + q), method = method,
                                      local.aniso = local.aniso, fix.tausq = fix.tausq,
                                      fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
          MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, sigmasq.local.init, rep(0, q) ),
                         fn = f_loglik, method = "L-BFGS-B",
                         lower = c( lam1.LB, lam2.LB, 0, sigmasq.local.LB, rep(-Inf, q) ),
                         upper = c( lam1.UB, lam2.UB, pi/2, sigmasq.local.UB, rep(Inf, q) ) )
          covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4], kappa)
          betaMLE <- MLEs$par[-c(1:4)]
        }
        if( method == "reml" ){
          f_loglik <- make_local_lik( locations = coords, cov.model = cov.model,
                                      data = data, Xmat = Xmat, nugg2.var = fixed.nugg2.var,
                                      fixed = rep(FALSE, 4), method = method,
                                      local.aniso = local.aniso, fix.tausq = fix.tausq,
                                      fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
          MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, sigmasq.local.init ),
                         fn = f_loglik, method = "L-BFGS-B",
                         lower = c( lam1.LB, lam2.LB, 0, sigmasq.local.LB ),
                         upper = c( lam1.UB, lam2.UB, pi/2, sigmasq.local.UB ) )
          covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4], kappa)
        }
      }
    }
  }
  if( cov.model != "matern" & cov.model != "cauchy" ){ # For covariance models without kappa
    if( !fix.tausq ){ # Estimate tausq
      if( method == "ml" ){
        f_loglik <- make_local_lik( locations = coords, cov.model = cov.model,
                                    data = data, Xmat = Xmat, nugg2.var = fixed.nugg2.var,
                                    fixed = rep(FALSE, 5 + q), method = method,
                                    local.aniso = local.aniso, fix.tausq = fix.tausq,
                                    fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
        MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, tausq.local.init,
                               sigmasq.local.init, rep(0, q) ),
                       fn = f_loglik, method = "L-BFGS-B",
                       lower = c( lam1.LB, lam2.LB, 0, tausq.local.LB, sigmasq.local.LB,
                                  rep(-Inf, q) ),
                       upper = c( lam1.UB, lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB,
                                  rep(Inf, q) ) )
        covMLE <- c(MLEs$par[1:5], NA)
        betaMLE <- MLEs$par[-c(1:5)]
      }
      if( method == "reml" ){
        f_loglik <- make_local_lik( locations = coords, cov.model = cov.model,
                                    data = data, Xmat = Xmat, nugg2.var = fixed.nugg2.var,
                                    fixed = rep(FALSE, 5), method = method,
                                    local.aniso = local.aniso, fix.tausq = fix.tausq,
                                    fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
        MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, tausq.local.init,
                               sigmasq.local.init ),
                       fn = f_loglik, method = "L-BFGS-B",
                       lower = c( lam1.LB, lam2.LB, 0, tausq.local.LB, sigmasq.local.LB ),
                       upper = c( lam1.UB, lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB ) )
        covMLE <- c(MLEs$par[1:5], NA)
      }
    }
    if( fix.tausq ){ # Fix tausq
      if( method == "ml" ){
        f_loglik <- make_local_lik( locations = coords, cov.model = cov.model,
                                    data = data, Xmat = Xmat, nugg2.var = fixed.nugg2.var,
                                    fixed = rep(FALSE, 4 + q), method = method,
                                    local.aniso = local.aniso, fix.tausq = fix.tausq,
                                    fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
        MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, sigmasq.local.init, rep(0, q) ),
                       fn = f_loglik, method = "L-BFGS-B",
                       lower = c( lam1.LB, lam2.LB, 0, sigmasq.local.LB, rep(-Inf, q) ),
                       upper = c( lam1.UB, lam2.UB, pi/2, sigmasq.local.UB, rep(Inf, q) ) )
        covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4], NA)
        betaMLE <- MLEs$par[-c(1:4)]
      }
      if( method == "reml" ){
        f_loglik <- make_local_lik( locations = coords, cov.model = cov.model,
                                    data = data, Xmat = Xmat, nugg2.var = fixed.nugg2.var,
                                    fixed = rep(FALSE, 4), method = method,
                                    local.aniso = local.aniso, fix.tausq = fix.tausq,
                                    fix.kappa = fix.kappa, tausq = tausq, kappa = kappa )
        MLEs <- optim( par = c(lam1.init, lam2.init, pi/4, sigmasq.local.init ),
                       fn = f_loglik, method = "L-BFGS-B",
                       lower = c( lam1.LB, lam2.LB, 0, sigmasq.local.LB ),
                       upper = c( lam1.UB, lam2.UB, pi/2, sigmasq.local.UB ) )
        covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4], NA)
      }
    }
  }
  #####################################################

  # Warnings?
  if( MLEs$convergence != 0 ){
    if( MLEs$convergence == 52 ){
      cat( paste("  There was a NON-FATAL error with optim(): \n  ",
                 MLEs$convergence, "  ", MLEs$message, "\n", sep = "") )
    }
    else{
      cat( paste("  There was an error with optim(): \n  ",
                 MLEs$convergence, "  ", MLEs$message, "\n", sep = "") )
    }
  }

  # Save all MLEs
  # Parameter order: lam1, lam2, eta, tausq, sigmasq, mu, kappa
  MLEs.save <- covMLE

  # Put the MLEs into a data frame
  names(MLEs.save) <- c( "lam1", "lam2", "eta", "tausq", "sigmasq", "kappa" )
  MLEs.save <- data.frame( t(MLEs.save) )

  aniso.mat <- kernel_cov( covMLE[1:3] )

  global.lik <-  -MLEs$value

  #===========================================================================
  # Calculate the covariance matrix using the MLEs
  #===========================================================================
  distances <- mahalanobis.dist( data.x = coords, vc = aniso.mat )
  NS.cov <- covMLE[5] * cov.spatial( distances, cov.model = cov.model,
                                       cov.pars = c(1,1), kappa = covMLE[6] )
  Data.Cov <- NS.cov + diag(rep(covMLE[4], N)) + fixed.nugg2.var

  #===========================================================================
  # Calculate the GLS estimate of beta
  #===========================================================================
  Data.Cov.chol <- chol(Data.Cov)

  if( method == "reml" ){
    tX.Cinv <- t(backsolve(Data.Cov.chol, backsolve(Data.Cov.chol, Xmat, transpose = TRUE)))
    beta.cov <- chol2inv( chol( tX.Cinv%*%Xmat) )/p
    beta.GLS <- (p*beta.cov %*% tX.Cinv %*% data)/p
    Mean.coefs <- data.frame( Estimate = beta.GLS,
                              Std.Error = sqrt(diag(beta.cov)),
                              t.val = beta.GLS/sqrt(diag(beta.cov)) )
  }
  if( method == "ml" ){
    beta.cov <- matrix(NA, q, q)
    beta.GLS <- betaMLE
    Mean.coefs <- data.frame( Estimate = beta.GLS,
                              Std.Error = sqrt(diag(beta.cov)),
                              t.val = beta.GLS/sqrt(diag(beta.cov)) )
  }


  #===========================================================================
  # Output
  #===========================================================================
  output <- list( MLEs.save = MLEs.save,
                  data = data,
                  beta.GLS = beta.GLS,
                  beta.cov = beta.cov,
                  Mean.coefs = Mean.coefs,
                  Cov.mat = Data.Cov,
                  Cov.mat.chol = Data.Cov.chol,
                  aniso.pars =  covMLE[1:3],
                  aniso.mat = aniso.mat,
                  tausq.est = covMLE[4],
                  sigmasq.est = covMLE[5],
                  kappa.MLE = covMLE[6],
                  fixed.nugg2.var = fixed.nugg2.var,
                  cov.model = cov.model,
                  coords = coords,
                  Xmat = Xmat,
                  global.loglik = global.lik,
                  fix.kappa = fix.kappa, kappa = kappa )

  class(output) <- "Aniso"

  return(output)

}

#======================================================================================
# Calculate predictions using the output of Aniso.fit()
#======================================================================================
# Using the output from Aniso.fit(), calculate the kriging predictors
# and kriging standard errors for prediction locations of interest.
#======================================================================================
#ROxygen comments ----
#' Obtain predictions at unobserved locations for the stationary
#' spatial model.
#'
#' \code{predict.Aniso} calculates the kriging predictor and corresponding
#' standard errors at unmonitored sites.
#'
#' @param object An "Aniso" object, from \code{Aniso_fit}.
#' @param pred.coords Matrix of locations where predictions are required.
#' @param pred.covariates Matrix of covariates for the prediction locations,
#' NOT including an intercept. The number of columns for this matrix must
#' match the design matrix from \code{mean.model} in \code{\link{NSconvo_fit}}.
#' Defaults to an intercept only.
#' @param pred.fixed.nugg2.var An optional vector or matrix describing the
#' the variance/covariance a fixed second nugget term (corresponds to
#' \code{fixed.nugg2.var} in \code{Aniso_fit}; often useful if conducting
#' prediction for held-out data). Defaults to zero.
#'
#' @return A list with the following components:
#' \item{pred.means}{Vector of the kriging predictor, for each location in
#' \code{pred.coords}.}
#' \item{pred.SDs}{Vector of the kriging standard errors, for each location
#' in \code{pred.coords}.}
#' @param ... additional arguments affecting the predictions produced.
#'
#' @examples
#' \dontrun{
#' pred.S <- predict( Aniso.obj,
#' pred.coords = cbind(runif(300),runif(300)) )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial
#' @importFrom StatMatch mahalanobis.dist

predict.Aniso <- function(object, pred.coords, pred.covariates = NULL,
                          pred.fixed.nugg2.var = NULL, ... )
{
  if( !inherits(object, "Aniso") ){
    warning("Object is not of type Aniso")
  }
  else{
    {
      M <- dim(pred.coords)[1]
      beta.MLE <- as.matrix(object$beta.GLS)
      tausq.est <- object$tausq.est
      sigmasq.est <- object$sigmasq.est
      if( object$fix.kappa ){
        kappa.MLE <- object$kappa
      } else{
        kappa.MLE <- object$kappa.MLE
      }
      Cov.mat.chol <- object$Cov.mat.chol
      data <- object$data
      N <- length(object$data)
      coords <- object$coords
      cov.model <- object$cov.model
      Xmat <- object$Xmat
      aniso.mat <- object$aniso.mat
      if (is.null(pred.covariates) == TRUE) {
        Xpred <- rep(1, M)
      }
      if (is.null(pred.covariates) == FALSE) {
        Xpred <- cbind(rep(1, M), pred.covariates)
      }
      if( is.null(pred.fixed.nugg2.var) == TRUE ){
        pred.fixed.nugg2.var <- matrix(0,N,N)
      } else{
        if( !is.matrix(pred.fixed.nugg2.var) ){ # Convert to a covariance matrix if specified as a vector
          pred.fixed.nugg2.var <- diag(pred.fixed.nugg2.var)
        }
      }

      data <- matrix(data, nrow = N)
      p <- dim(data)[2]
      CC.distances <- mahalanobis.dist(data.x = pred.coords, data.y = coords,
                                       vc = aniso.mat)
      CrossCov <- sigmasq.est * cov.spatial(CC.distances, cov.model = cov.model,
                                            cov.pars = c(1, 1), kappa = kappa.MLE)
      crscov.Cinv <- t(backsolve(Cov.mat.chol,
                                 backsolve(Cov.mat.chol, t(CrossCov), transpose = TRUE)))

      pred.means <- matrix(NA, M, p)
      for (i in 1:p) {
        pred.means[, i] <- Xpred %*% beta.MLE + crscov.Cinv %*% (object$data[, i] - Xmat %*% beta.MLE)
      }
      pred.SDs <- sqrt(rep(tausq.est + sigmasq.est, M) + diag(pred.fixed.nugg2.var) - diag(crscov.Cinv %*% t(CrossCov)))

      output <- list(pred.means = pred.means, pred.SDs = pred.SDs)
      return(output)
    }
  }
}
