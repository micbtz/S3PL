\name{print.tpl}
\alias{print.tpl}
\title{Print estimates of the 3PL model}
\description{
Print method for the output of function \code{\link{tpl}}.
}
\usage{
\method{print}{tpl}(x, IRTparam = TRUE, ...)
}
\arguments{
 \item{x}{an object of class \code{tpl} obtained from function \code{\link{tpl}}.}
 \item{IRTparam}{logical. If \code{TRUE} the usual IRT parameterization is applied. See details for more information.}
 \item{\dots}{further arguments passed to or from other methods.}
}
\details{
The parameter of the 3PL model are estimated using the following parameterization:
\deqn{\pi_{ij} = F(\beta_{3j}) + (1 - F(\beta_{3j})) 
\frac{\exp(\beta_{1i} + \beta_{2i} z)}{1 + 
\exp(\beta_{1j} + \beta_{2j} \theta_i)}, }{ \pi_{ij} = 
F(\beta_{3j}) + (1 - F(\beta_{3j})) * \{exp(\beta_{1j} + \beta_{2j} \theta_i)\}/
\{1 + exp(\beta_{1j} + \beta_{2j} \theta_i)\}, } 
where \eqn{\pi_{ij}} denotes the conditional probability of responding correctly to the \eqn{j}th item given \eqn{\theta_i}, and 
\deqn{F(\beta_{3j}) = \exp(F(\beta_{3j}))/(1-\exp(F(\beta_{3j})))}.
If argument \code{IRTparam} is \code{TRUE} the estimates are transformed to
the usual IRT parameterization, that is, 
\deqn{\pi_{ij} = c_j + (1 - c_j) \frac{\exp[a_i (\theta_i - b_j)]}{1 + 
\exp[a_j (\theta_i - b_j)]}, }{ \pi_i = c_j + (1 - c_j) * 
exp\{a_i (\theta_i - b_j)\} / [1 + exp\{a_j (\theta - b_j)\}], }
where \eqn{a_j = \beta_{2j}} and \eqn{b_j = -\beta_{1j}/\beta_{2j}}.
}
%\references{
%}
\seealso{
	\code{\link{tpl}}
}


