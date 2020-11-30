\name{traitstats.rbdsummary}
\alias{traitstats.rbdsummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RCBD Data Analysis Report
}
\description{
RCBD data analysis report in PDF formate compiling the ANOVA, Descriptive Statistics and Genetic Parameter along with Data distribution graph.
}
\usage{
traitstats.rbdsummary(Treatment, Replication, DataFile)
}
\arguments{
  \item{Treatment}{Treatment column in the DataFile}
	\item{Replication}{Replication column in the DataFile}
	\item{DataFile}{Input Data file name}
}
\details{
traitstats.rbdsummary function analysis the input DataFile according to the methodology given by Panse and Sukhatme (1984) and generates the RCBD data analysis results comprising of Data point distribution, ANOVA table, Descriptive statistics and Genetic parameter arranged in the PDF format output file.
}
\value{
\code{Data distribution graph}\cr
	\code{ANOVA table : Rep, Trt, MStrt, MSrep and MSerror along the test of significance}\cr
	\code{Descriptive statistics table : Trait Mean, Min, Max, SE(m), SE(d), CV and CD.}\cr
	\code{Genetic Parameter : GCV, PCV, h2, GA and GAM along with the category}\cr
}
\references{
Panse and Sukhatme. 1984. Statistical Methods for Agricultural Workers. ICAR, New Delhi.
}
\author{
Nitesh, S.D., Parashuram Patroti and Shilpa Prashuram
}
\seealso{
\code{\link{traitstats.rbd},\link{traitstats.rbdmeantable}, \link{traitstats.rbdstats}}
}
\examples{
data("Barley")
traitstats.rbd (Barley$trt,
                Barley$rep,
                Barley)
}
