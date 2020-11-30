\name{traitstats.rbdstats}
\alias{traitstats.rbdstats}
\title{
RCBD Descriptive Statistics
}
\description{
Descriptive Statistics for Trait Mean, Minimum, Maximum, SE(m), SE(d) and CV.
}
\usage{
traitstats.rbdstats(Treatment, Replication, DataFile)
}
\arguments{
  \item{Treatment}{Treatment column in the DataFile}
	\item{Replication}{Replication column in the DataFile}
	\item{DataFile}{Input Data file name}
}
\details{
traitstats.rbdstats function analysis the input DataFile according to the methodology given by Panse and Sukhatme (1984) and generates the RBD descriptive statistics for Trait Mean, Minimum, Maximum, Standard Error(mean), Standard Error(deviation) and Coefficient of Variation.
}
\value{
	\code{Trait Mean}\cr
	\code{Minimum}\cr
	\code{Maximum}\cr
	\code{StdError(m)}\cr
	\code{StdError(d)}\cr
	\code{CV}\cr
}
\references{
Panse and Sukhatme. 1984. Statistical Methods for Agricultural Workers. ICAR, New Delhi.
}
\author{
Nitesh, S.D, Parashuram Patroti and Shilpa Parashuram
}

\seealso{
\code{\link{traitstats.rbd},\link{traitstats.rbdmeantable},\link{traitstats.rbdsummary}}
}
\examples{
data("Barley")
traitstats.rbdstats(Barley$trt,
                    Barley$rep,
                    Barley)
}
