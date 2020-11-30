\name{traitstats.phenocorl}
\alias{traitstats.phenocorl}
\title{
Phenotypic Correlation Analysis
}
\description{
An extract of phenotypic correlation coefficient for between traits
}
\usage{
traitstats.phenocorl(Treatment, Replication, DataFile)
}
\arguments{
  \item{Treatment}{Treatment column in the DataFile}
	\item{Replication}{Replication column in the DataFile}
	\item{DataFile}{Input Data file name}
}
\details{
Phenotypic Correlation coefficient is calculated according to the Al-Jibouri et al., (1958).
}
\value{
Phenotypic correlation coefficient matrix.
}
\references{
Al-Jibour HA, Miller PA and Robinson HP. 1958. Genotypic and environmental variance in upland cotton cross of interspecific origin. Agronomy Journal, 50:633-637.
}
\author{
Nitesh, S.D., Parashuram Patroti, and Shilpa Parashuram.
}
\examples{
data("Barley")
traitstats.phenocorl(Treatment=Barley$trt,
                  Replication=Barley$rep,
                  DataFile=Barley)
  }
