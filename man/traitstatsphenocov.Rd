\name{traitstats.phenocov}
\alias{traitstats.phenocov}
\title{
Phenotypic Variance and Covariance
}
\description{
Computation of phenotypic variance and covariance matrix.
}
\usage{
traitstats.phenocov(Treatment, Replication, DataFile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Treatment}{Treatment column in the DataFile}
	\item{Replication}{Replication column in the DataFile}
	\item{DataFile}{Input Data file name}
}
\details{
Phenotypic variance and covariance matrix is computed as per the methodology given by Al-Jibour et al., (1958).
}
\value{
Phenotypic variance and co-variance matrix. Diagonal value indicating the phenotypic variance and upper and lower half of the matrix indicating the phenotypic covariance.
}
\references{
Al-Jibour HA, Miller PA and Robinson HP. 1958. Genotypic and environmental variance in upland cotton cross of interspecific origin. Agronomy Journal, 50:633-637.
}
\author{
Nitesh, S.D., Parashuram Patroti and Shilpa Parashuram
}

\examples{
data("Barley")
traitstats.phenocov(Treatment=Barley$trt,
                    Replication=Barley$rep,
                    DataFile=Barley)
}
