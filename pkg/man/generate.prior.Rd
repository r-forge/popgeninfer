\name{generate.prior}
\Rdversion{0.8}
\alias{generate.prior}
\title{Generate Prior Values for the Parameters theta and M}
\description{
Generate a set of prior values for the parameters theta and M
}
\usage{
  generate.prior(number.of.simulations,prior.theta,min.theta,max.theta,prior.M,
    min.M,max.M)
}
\arguments{
  \item{number.of.simulations}{total number of simulations}
  \item{prior.theta}{a character string indicating the prior distribution for theta. Possible values are "UNI" (uniform prior) and"LOG" (log-unifomr prior)}
  \item{min.theta}{Minimum value for theta}
  \item{max.theta}{Maximum value for theta}
  \item{prior.M}{a character string indicating the prior distribution for M. Possible values are "UNI" (uniform prior) and"LOG" (log-unifomr prior)}
  \item{min.M}{Minimum value for M}
  \item{max.M}{Maximum value for M}
}
\details{
Generate a set of prior values for the parameters theta and M, which will be used by the coalescent-based simulator to simulate artificial data.
}
\value{
	a matrix of prior values.
}
\examples{
## This is to generate a prior distribution of the model parameters.

prior <- generate.prior(number.of.simulations = 1e3,prior.theta = 'UNI',min.theta = 0.1,
  max.theta = 5,prior.M = 'UNI',min.M = 0.1,max.M = 5)

prior
}
