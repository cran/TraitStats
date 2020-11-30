\name{traitstats.envcov}
\alias{traitstats.envcov}
\title{
Environmental Variance and Covariance
}
\description{
Computation of environment variance and covariance matrix.
}
\usage{
traitstats.envcov(Treatment, Replication, DataFile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Treatment}{Treatment column in the DataFile}
	\item{Replication}{Replication column in the DataFile}
	\item{DataFile}{Input Data file name}
}
\details{
Environmental variance and covariance matrix is computed as per the methodology given by Al-Jibour et al., (1958).
}
\value{
Environment variance and co-variance matrix. Diagonal value indicating the environmental variance and upper and lower half of the matrix indicating the environmental covariance.
}
\references{
Al-Jibour HA, Miller PA and Robinson HP. 1958. Genotypic and environmental variance in upland cotton cross of interspecific origin. Agronomy Journal, 50:633-637.
}
\author{
Nitesh, S.D., Parashuram Patroti and Shilpa Parashuram
}

\examples{
data("Barley")
traitstats.envcov(Treatment=Barley$trt,
                  Replication=Barley$rep,
                  DataFile=Barley)
}
