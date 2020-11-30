\name{traitstats.envcorl}
\alias{traitstats.envcorl}
\title{
Environmental Correlation Analysis
}
\description{
An extract of environment correlation coefficient for between traits
}
\usage{
traitstats.envcorl(Treatment, Replication, DataFile)
}
\arguments{
  \item{Treatment}{Treatment column in the DataFile}
	\item{Replication}{Replication column in the DataFile}
	\item{DataFile}{Input Data file name}
}
\details{
Environmental Correlation coefficient is calculated according to the Al-Jibouri et al., (1958).
}
\value{
Environmental correlation coefficient matrix.
}
\references{
Al-Jibour HA, Miller PA and Robinson HP. 1958. Genotypic and environmental variance in upland cotton cross of interspecific origin. Agronomy Journal, 50:633-637.
}
\author{
Nitesh, S.D., Parashuram Patroti, and Shilpa Parashuram.
}
\examples{
data("Barley")
traitstats.envcorl(Treatment=Barley$trt,
                  Replication=Barley$rep,
                  DataFile=Barley)
  }

