\name{S3PL-package}
\alias{S3PL-package}
\alias{S3PL}
\docType{package}
\title{
Shrinkage estimation of the three-parameter logistic model
}
\description{
Estimation of the Three Parameter Logistic (3PL) model using penalized likelihood estimation and bias reduction methodology.
}
%\details{
%}
\author{
\packageAuthor{S3PL}

Maintainer: \packageMaintainer{S3PL}
}
\references{
Battauz, M. and Bellio, R. (2021). Shrinkage estimation of the three-parameter logistic model. Br J Math Stat Psychol. https://doi.org/10.1111/bmsp.12241
}
%\keyword{ package }
\examples{
# simulate the data
set.seed(1)
dataset<-simtpl(n = 1000, diff = rnorm(30), discr = rbeta(30, 3, 3) * 0.6 + 0.7, 
  guess = rbeta(30, 3, 5) * 0.4)

# estimate the parameters using penalized maximum likelihood
mod<-tpl(data = dataset, penalized = TRUE, lambda = seq(0.05, 0.5, by = 0.05), cv = TRUE)
mod
print(mod, IRTparam = FALSE)

# plot the regularization paths for the penalized estimates
plot(mod)
plot(mod, IRTparam = FALSE)
plot(mod, onlyGuess = FALSE)
plot(mod, IRTparam = FALSE, onlyGuess = FALSE)

# estimate the parameters using bias reduction methodology
modBR<-tpl(data = dataset, BR = TRUE)
modBR
print(modBR, IRTparam = FALSE)
}
