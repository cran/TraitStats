\name{traitstats.rbd}
\alias{traitstats.rbd}
\title{
RCBD Bulk Data Analysis
}
\description{
Bulk RCBD data analysis for ANOVA and Descriptive statistics
}
\usage{
traitstats.rbd(Treatment, Replication, DataFile)
}
\arguments{
  \item{Treatment}{Treatment column in the DataFile}
	\item{Replication}{Replication column in the DataFile}
	\item{DataFile}{Input Data file name}
}
\details{
traitstats.rbd function analysis the input DataFile according to the methodology given by Panse and Sukhatme (1984). Futhermore, it also represent the analysis results along with descriptive statistics parameters like Trait mean, Minimum, Maximum, CV, StdError(m) and StdError(d) in the tabular representation of the  publication ready format.
}
\value{
\code{Number of Replication}\cr
	\code{Number of Treatment}\cr
	\code{Treatment Mean Sum of Square along with Test of Significance}\cr
	\code{Replication Mean Sum of Square along with Test of Significance}\cr
\code{Error Mean Sum of Square}\cr
	\code{Trait Mean}\cr
	\code{Minimum}\cr
	\code{Maximum}\cr
	\code{CV}\cr
	\code{StdError(m)}\cr
	\code{StdError(d)}
}
\references{
Panse and Sukhatme. 1984. Statistical Methods for Agricultural Workers. ICAR, New Delhi.
}
\author{
Nitesh, S.D., Parashuram Patroti and Shilpa Parashuram
}
\seealso{
\code{\link{traitstats.rbdmeantable}, \link{traitstats.rbdstats},\link{traitstats.rbdsummary}}
}
\examples{
data("Barley")
traitstats.rbd (Barley$trt,
                Barley$rep,
                Barley)

}
