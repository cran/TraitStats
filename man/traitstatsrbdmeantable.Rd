\name{traitstats.rbdmeantable}
\alias{traitstats.rbdmeantable}
\title{
RCBD Mean Table
}
\description{
Individual treatment mean table along with Grand trait mean, CV, SE(m) and CD.
}
\usage{
traitstats.rbdmeantable(Treatment, Replication, DataFile)
}
\arguments{
  \item{Treatment}{Treatment column in the DataFile}
	\item{Replication}{Replication column in the DataFile}
	\item{DataFile}{Input Data file name}
}
\details{
traitstats.rbdmeantable function analysis the input DataFile according to the methodology given by Panse and Sukhatme (1984) and generates the data frame consisting of the treatment mean over replication. Furthermore, it also represent the descriptive statistics parameter – Trait mean, StdError(m), CV and CD at 5 percent level of significance in the tabular representation of the publication ready format.
}
\value{
\code{Individual Treatment mean}\cr
\code{Trait Mean}\cr
\code{SE(m)}\cr
\code{CV}\cr
\code{CD}
}
\references{
Panse and Sukhatme. 1984. Statistical Methods for Agricultural Workers. ICAR, New Delhi.
}
\author{
Nitesh, S.D, Parashuram Patroti and Shilpa Parashuram
}
\seealso{
\code{\link{traitstats.rbd}, \link{traitstats.rbdstats},\link{traitstats.rbdsummary}},
}
\examples{
data("Barley")
traitstats.rbdmeantable(Barley$trt,
                        Barley$rep,
                        Barley)
}
