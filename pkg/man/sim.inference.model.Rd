\name{sim.inference.model}
\Rdversion{0.8}
\alias{sim.inference.model}
\title{Simulate a Sample of Genes}
\description{
Simulate a sample of genes.
}
\usage{
  sim.inference.model(number.of.sampled.demes,sample.sizes,M,pi)
}
\arguments{
  \item{number.of.sampled.demes}{number of demes sampled}
  \item{sample.sizes}{number of genes sampled in each deme}
  \item{M}{twice the product of effective size of each deme and migration rate}
  \item{pi}{vector of allele A frequencies at each locus from the migrant pool}
}
\details{
Simulate a sample of genes, using a diffusion approximation for the probability of gene counts within each deme.
}
\value{
	a \sQuote{sample} object
}
\examples{
## This is to simulate a sample of genes.

sample <- sim.inference.model(number.of.sampled.demes = 10,sample.sizes = 50,M = 2,pi = runif(2))
sample
}
