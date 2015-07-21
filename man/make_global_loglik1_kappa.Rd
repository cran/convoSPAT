% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/convoSPAT_paramEst.R
\name{make_global_loglik1_kappa}
\alias{make_global_loglik1_kappa}
\title{Constructor functions for global parameter estimation.}
\usage{
make_global_loglik1_kappa(data, Xmat, cov.model, Scalemat, Distmat)
}
\arguments{
\item{data}{A vector or matrix of data to use in the likelihood
calculation.}

\item{Xmat}{The design matrix for the mean model.}

\item{cov.model}{String; the covariance model.}

\item{Scalemat}{Matrix; contains the scaling quantities from the
covariance function.}

\item{Distmat}{Matrix; contains the scaled distances.}
}
\value{
This function returns another function for use in \code{optim}.
}
\description{
This function generates another function to be used within \code{optim} to
obtain maximum likelihood estimates of
global variance parameters tausq, sigmasq, and nu.
}
\examples{
\dontrun{
make_global_loglik1_kappa( data, Xmat, cov.model, Scalemat, Distmat )
}
}
