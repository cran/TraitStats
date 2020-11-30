\name{traitstats.genocov}
\alias{traitstats.genocov}
\title{
Genotypic Variance and Covariance
}
\description{
Computation of genotypic variance and covariance matrix.
}
\usage{
traitstats.genocov(Treatment, Replication, DataFile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Treatment}{Treatment column in the DataFile}
	\item{Replication}{Replication column in the DataFile}
	\item{DataFile}{Input Data file name}
}
\details{
Genotypic variance and covariance matrix is computed as per the methodology given by Al-Jibour et al., (1958).
}
\value{
Genotypic variance and co-variance matrix. Diagonal value indicating the genetic variance and upper and lower half of the matrix indicating the genotypic covariance.
}
\references{
Al-Jibour HA, Miller PA and Robinson HP. 1958. Genotypic and environmental variance in upland cotton cross of interspecific origin. Agronomy Journal, 50:633-637.
}
\author{
Nitesh, S.D., Parashuram Patroti and Shilpa Parashuram
}
\examples{
data("Barley")
traitstats.genocov(Treatment=Barley$trt,
                  Replication=Barley$rep,
                  DataFile=Barley)
}
