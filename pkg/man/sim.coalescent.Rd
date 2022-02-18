\name{sim.coalescent}
\Rdversion{0.8}
\alias{sim.coalescent}
\title{Simulate a Sample of Genes Using the Coalescent}
\description{
Simulate a sample of genes using the coalescent.
}
\usage{
  sim.coalescent(total.number.of.demes,effective.size,migration.rate,mutation.rate,
    number.of.sampled.demes,sample.sizes,number.of.loci)
}
\arguments{
  \item{total.number.of.demes}{total number of demes in the island model}
  \item{effective.size}{effective population size of each deme}
  \item{migration.rate}{rate of migration among demes}
  \item{mutation.rate}{mutation rate. The model assumes a 2-alleles model of mutations}
  \item{number.of.sampled.demes}{number of demes sampled}
  \item{sample.sizes}{number of genes sampled in each deme}
  \item{number.of.loci}{number of loci}
}
\details{
Simulate a sample of genes, using a coalescent-based simulator.
}
\value{
	a \sQuote{sample} object
}
\examples{
## This is to simulate a sample of genes.

sample <- sim.coalescent(total.number.of.demes = 100,effective.size = 1000,
  migration.rate = 0.001,mutation.rate = 0.000001,number.of.sampled.demes = 10,
  sample.sizes = 50,number.of.loci = 2)
sample
}
