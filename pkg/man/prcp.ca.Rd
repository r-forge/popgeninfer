\name{prcp.ca}
\Rdversion{0.8}
\alias{prcp.ca}
\title{Perform a Principal Component Analysis}
\description{
Perform a principal component analysis (PCA) based on the simulated summary statistics, together with those measured on the observed data.
}
\usage{
	prcp.ca(stats,target,...)
}
\arguments{
  \item{stats}{a data frame containing a set of summary statistics measured on simulated data}
  \item{target}{a data frame containing a set of summary statistics measured on the observed data}
  \item{...}{options for plotting}
}
\details{
Perform a principal component analysis (PCA) based on the simulated summary statistics, together with those measured on the observed data.
}
\value{
	an object of the class pca dudi, and a graphical representation of the PCA.
}
\examples{
	
## This is to generate a pseudo-observed dataset.

prior <- generate.prior(number.of.simulations = 1,prior.theta = 'UNI',min.theta = 1.0,
  max.theta = 1.0,prior.M = 'UNI',min.M = 1.0,max.M = 1.0)

target <- sim.island.model(number.of.simulations = 1,mutation.model = 'SMM',
  total.number.of.demes = 10,number.of.loci = 20,number.of.sampled.demes = 10,
  sample.sizes = 50)

## This is to generate a prior distribution of the model parameters.

prior <- generate.prior(number.of.simulations = 1e3,prior.theta = 'UNI',min.theta = 0.1,
  max.theta = 5,prior.M = 'UNI',min.M = 0.1,max.M = 5)

## This is to generate summary statistics from simulated data.

stats <- sim.island.model(number.of.simulations = 1e3,mutation.model = 'SMM',
  total.number.of.demes = 10,number.of.loci = 20,number.of.sampled.demes = 10,
  sample.sizes = 50)

## Perform a principal component analysis (PCA) 

prcp.ca(stats,target)

}
