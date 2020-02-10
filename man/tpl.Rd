\name{tpl}
\alias{tpl}
\title{Shrinkage estimation of a 3PL model}
\description{
Estimation of a 3PL model using maximum likelihood estimation, penalized maximum likelihood
estimation and the bias reduction methodology.
}
\usage{
tpl(data, penalized = FALSE, BR = FALSE, lambda = 0, cv = FALSE, K = 5, nq = 21, 
    empirical = TRUE, M = 5000, ncores = NULL, trace = TRUE, init_BR = FALSE)
}
\arguments{
 \item{data}{a dataframe.}
 \item{penalized}{logical. If \code{TRUE} ridge-type penalization is applied.}
 \item{BR}{logical. If \code{TRUE} the bias reduction methodology is applied.}
 \item{lambda}{vector of tuning parameters to perform ridge-type penalization.}
 \item{cv}{logical. If \code{TRUE} K-fold cross-validation is performed to select \eqn{\lambda}.}
 \item{K}{number of folds for cross-validation.}
 \item{nq}{number of quadrature points.}
 \item{empirical}{logical. If \code{TRUE} the empirical approximation is used to compute the expected values when applying bias reduction methodology. If \code{FALSE} the Monte Carlo approximation is used.}
 \item{M}{Number of replications used for Monte Carlo approximation when applying bias reduction methodology if argument \code{empirical} is \code{FALSE}.}
 \item{ncores}{Number of cores used parallelize computations when applying bias reduction methodology. In Windows OS it is set to 1.}
 \item{trace}{logical. If \code{TRUE} tracing information is produced.}
 \item{init_BR}{logical. If \code{TRUE} the MLE is used as starting point for bias reduced estimation.}
}
\value{
An object of class \code{tpl} with components:
\item{est_MLE}{vector of maximum likelihood estimates.}
\item{logLik}{log-likelihood value at maximum likelihood estimates.}
\item{conv_MLE}{convergence of maximum likelihood estimation.}
\item{est_pen}{matrix of penalized maximum likelihood estimates obtained for each \eqn{\lambda}.}
\item{conv_pen}{vector convergence of maximum likelihood estimation; one for each \eqn{\lambda}.}
\item{est_pen_sel}{vector of penalized maximum likelihood estimates at the selected value of \eqn{\lambda}.}
\item{cv_est}{list of matrixes of item parameter estimates for each value of \eqn{\lambda}. Each matrix contains the estimates for each of the K folds.}
\item{cv_lik}{list of vectors of log-likelihood values computed on the validation set for each value of \eqn{\lambda}.}
\item{lambda_sel}{selected value of \eqn{\lambda}.}
\item{lambda}{vector of values for \eqn{\lambda}.}
\item{est_br}{vector of estimates obtained with the bias reduction methodology.}
\item{conv_br}{convergence of the bias reduction method.}
\item{itemnames}{column names of \code{data}.}
%\item{}{}
}
%\details{
%}
%\references{
%}
\examples{
# simulate the data
set.seed(1)
dataset<-simtpl(n = 1000, diff = rnorm(30), discr = rbeta(30, 3, 3)*0.6+0.7, 
  guess = rbeta(30, 3, 5)*0.4)

# estimate the parameters using penalized maximum likelihood
mod<-tpl(data = dataset, penalized = TRUE, lambda = seq(0.05, 0.5, by = 0.05), 
  cv = TRUE)
mod
print(mod, IRTparam = FALSE)

# estimate the parameters using bias reduction methodology
modBR<-tpl(data = dataset, BR = TRUE, ncores = 2)
modBR
print(modBR, IRTparam = FALSE)

}
\author{Michela Battauz, Ruggero Bellio}

