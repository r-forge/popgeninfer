\name{expected.fst}
\Rdversion{0.8}
\alias{expected.fst}
\title{Compute the Expected \eqn{F_{ST}}}
\description{
Compute the expected \eqn{F_{ST}} in a finite island model.
}
\usage{
	expected.fst(total.number.of.demes,effective.size,migration.rate,mutation.rate)
}
\arguments{
  \item{total.number.of.demes}{total number of demes in the island model}
  \item{effective.size}{effective population size of each deme}
  \item{migration.rate}{rate of migration among demes}
  \item{mutation.rate}{mutation rate. The model assumes a 2-alleles model of mutations}
}
\details{
Computes the expected \eqn{F_{ST}} in a finite island model of population structure, with equal sized demes and isotropic migration across demes. The mutation model is a 2-alleles model. This function is based upon eq.(6) in Rousset (1997), adapted for haploids.
}
\value{
	Expected \eqn{F_{ST}} in a finite island model of population structure, with equal sized demes and isotropic migration across demes}
\examples{
## This is to get the expected level of population differentiation in an
## island model with 100 demes, each of size 1000, that exchange migrants
## at a rate of 0.001. The mutation rate is assumed to be 0.000001.

expected.fst(total.number.of.demes = 100,effective.size = 1000,migration.rate = 0.001,
  mutation.rate = 0.000001)
}
\references{
Rousset, F. (1996) Equilibrium values of measures of population subdivision for stepwise mutation processes, \emph{Genetics} \bold{142}: 1357-1362.
}
