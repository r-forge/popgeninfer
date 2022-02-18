\name{sim.island.model}
\Rdversion{0.8}
\alias{sim.island.model}
\title{Simulate Artificial Data from an Island Model}
\description{
Simulate artificial data using the coalescent, and return summary statistics from each dataset.
}
\usage{
  sim.island.model(number.of.simulations,mutation.model,total.number.of.demes,
    number.of.loci,number.of.sampled.demes,sample.sizes)
}

\arguments{
  \item{number.of.simulations}{total number of simulations}
  \item{mutation.model}{a character string indicating the mutation model. Possible values are "IAM" (infinite island model),"KAM" (K-allele model) and "SMM" (stepwise mutation model)}
  \item{total.number.of.demes}{total number of demes in the island model}
  \item{number.of.loci}{total number of simulated loci}
  \item{number.of.sampled.demes}{number of demes sampled}
  \item{sample.sizes}{number of genes sampled in each deme}
}
\details{
Simulate artificial data using a coalescent-based simulator, and return summary statistics from each dataset.
}
\value{
	a matrix of summary statistics.
}
\examples{
## This is to generate a prior distribution of the model parameters.

prior <- generate.prior(number.of.simulations = 1e3,prior.theta = 'UNI',min.theta = 0.1,
  max.theta = 5,prior.M = 'UNI',min.M = 0.1,max.M = 5)

## This is to generate summary statistics from simulated data.

stats <- sim.island.model(number.of.simulations = 1e3,mutation.model = 'SMM',
  total.number.of.demes = 10,number.of.loci = 20,number.of.sampled.demes = 10,
  sample.sizes = 50)

stats
}
