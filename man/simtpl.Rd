\name{simtpl}
\alias{simtpl}
\title{Simulation of dataset}
\description{
Simulation of a dataset from a 3PL model
}
\usage{
simtpl(n, diff, discr, guess)
}
\arguments{
 \item{n}{sample size.}
 \item{diff}{vector of difficulty parameters.}
 \item{discr}{vector of discrimination parameters.}
 \item{guess}{vector of guessing parameters.}
}
\details{
The three vectors of item parameters should have the same lenght,
which is equal to the number of items.
The number of items is denoted by \eqn{J}.
}
\value{
\code{data.frame} with dimension \eqn{n X J} containing the binary responses.
}
\examples{
dataset <- simtpl(n = 1000, diff = rnorm(30), discr = rbeta(30, 3, 3) * 0.6 + 0.7, 
  guess = rbeta(30, 3, 5) * 0.4)
head(dataset)
}


