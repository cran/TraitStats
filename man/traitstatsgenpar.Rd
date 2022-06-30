\name{traitstats.genpar}
\alias{traitstats.genpar}
\title{
Genetic Parameter Analysis
}
\description{
traitstats.genpar performs genetic variability analysis for all the traits in the DataFile.
}
\usage{
traitstats.genpar(Treatment, Replication, DataFile)
}
\arguments{
  \item{Treatment}{Treatment column in the DataFile}
	\item{Replication}{Replication column in the DataFile}
	\item{DataFile}{Input Data file name}
}
\details{
\code{traitstats.genpar} performs genetic parameter analysis from the input values extracted from the ANOVA analysis for individual traits and computes several variability estimates. .

The phenotypic, genotypic and environmental variance
(\ifelse{html}{\out{<i>&sigma;<sup>2</sup> p </i>}}{\eqn{\sigma^{2}_{p}}},
\ifelse{html}{\out{<i>&sigma;<sup>2</sup> g </i>}}{\eqn{\sigma^{2}_{g}}}
and
\ifelse{html}{\out{<i>&sigma;<sup>2</sup> e </i>}}{\eqn{\sigma^{2}_{e}}}
) are obtained from the ANOVA tables according to the expected value of mean
square described by Federer and Searle (1976) as follows:

\ifelse{html}{\out{<p style="text-align:
center;">&sigma;<sup>2</sup> p  = Mean sum of squares of
test treatments</p>}}{\deqn{\sigma^{2}_{p} = \textrm{Mean sum of squares of
test treatments}}} \ifelse{html}{\out{<p style="text-align:
center;"> &sigma;<sup>2</sup> e   = Mean sum of squares of
residuals</p>}}{\deqn{\sigma^{2}_{e} = \textrm{Mean sum of squares of
residuals}}} \ifelse{html}{\out{<p style="text-align:
center;">&sigma;<sup>2</sup> g  =
&sigma;<sup>2</sup> p  &minus;
&sigma;<sup>2</sup> e </p>}}{\deqn{\sigma^{2}_{g} =
\sigma^{2}_{p} - \sigma^{2}_{e}}}

Phenotypic and genotypic coefficients of variation
(\ifelse{html}{\out{<i>PCV</i>}}{\eqn{PCV}} and
\ifelse{html}{\out{<i>GCV</i>}}{\eqn{GCV}}) are estimated according to Burton
(1951, 1952) as follows:

\ifelse{html}{\out{<p style="text-align: center;">PCV =
 [ <sup>&sigma;<sup>2</sup> p </sup> &frasl;  &radic;
(<span style="text-decoration:overline">x</span>)  ]
&times; 100</p>}}{\deqn{PCV = \frac{\sigma^{2}_{p}}{\sqrt{\overline{x}}}
\times 100}} \ifelse{html}{\out{<p style="text-align: center;"> GCV =
 [ <sup>&sigma;<sup>2</sup> g </sup> &frasl;  &radic;
(<span style="text-decoration:overline">x</span>)  ]
&times; 100</p>}}{\deqn{GCV = \frac{\sigma^{2}_{g}}{\sqrt{\overline{x}}}
\times 100}}

Where \ifelse{html}{\out{<i><span style="text-decoration:
overline;">x</span></i>}}{\eqn{\overline{x}}} is the mean.

The estimates of \ifelse{html}{\out{<i>PCV</i>}}{\eqn{PCV}} and
\ifelse{html}{\out{<i>GCV</i>}}{\eqn{GCV}} are categorised according to
Sivasubramanian and Madhavamenon (1978) as follows:

\tabular{ll}{ \strong{\emph{CV} (\%)} \tab \strong{Category} \cr x
\ifelse{html}{\out{<}}{\eqn{<}} 10 \tab Low \cr 10
\ifelse{html}{\out{&le;}}{\eqn{\le}} x \ifelse{html}{\out{<}}{\eqn{<}} 20 \tab
Medium \cr \ifelse{html}{\out{&ge;}}{\eqn{\ge}} 20 \tab High }

The broad-sense heritability
(\ifelse{html}{\out{<i>H<sup>2</sup></i>}}{\eqn{H^{2}}}) is calculated
according to method of Lush (1940) as follows:

\ifelse{html}{\out{<p style="text-align: center;"> H<sup>2</sup> =
<sup>&sigma;<sup>2</sup> g </sup> &frasl;
 &sigma;<sup>2</sup> p   </p>}}{\deqn{H^{2} =
\frac{\sigma^{2}_{g}}{\sigma^{2}_{p}}}}

The estimates of broad-sense heritability
(\ifelse{html}{\out{<i>H<sup>2</sup></i>}}{\eqn{H^{2}}}) are categorised
according to Robinson (1966) as follows:

\tabular{ll}{ \strong{\ifelse{html}{\out{<i>H<sup>2</sup></i>}}{\eqn{H^{2}}}}
\tab \strong{Category} \cr x \ifelse{html}{\out{<}}{\eqn{<}} 30 \tab Low \cr
30 \ifelse{html}{\out{&le;}}{\eqn{\le}} x \ifelse{html}{\out{<}}{\eqn{<}} 60
\tab Medium \cr \ifelse{html}{\out{&ge;}}{\eqn{\ge}} 60 \tab High }

Genetic advance (\ifelse{html}{\out{<i>GA</i>}}{\eqn{GA}}) is estimated
and categorised according to Johnson et al., (1955) as follows:

\ifelse{html}{\out{<p style="text-align: center;"> GA = k &times;
&sigma; g  &times;  [ <sup>H<sup>2</sup></sup> &frasl;
 100  ]  </p>}}{\deqn{GA = k \times \sigma_{g} \times
\frac{H^{2}}{100}}}

Where the constant \ifelse{html}{\out{<i>k</i>}}{\eqn{k}} is the standardized
selection differential or selection intensity. The value of
\ifelse{html}{\out{<i>k</i>}}{\eqn{k}} at 5\% proportion selected is 2.063.
Values of \ifelse{html}{\out{<i>k</i>}}{\eqn{k}} at other selected proportions
are available in Appendix Table A of Falconer and Mackay (1996).
}
\value{
Genotypic Coefficient of Variation and GCV Category
	\cr Phenotypic Coefficient of Variation and PCV Category
	\cr Heritability (broad sense)and h2 Category
	\cr Genetic Advance and GA Category
	\cr Genetic Advance percent Mean and GAM Category

}
\references{
Lush JH. (1940). Intra-sire correlations or regressions of offspring on dam as a method of estimating heritability of characteristics. Proceedings of the American Society of Animal Nutrition, 1940(1):293-301.

Burtone GW and De Vane GM. (1953). Estimating heritability in tall Fescus (Festuca arundinaceae) from replicated clonal material. Agronomy Journal,45:478-481.

Johnson HW, Robinson HF and Comstock RE. (1955). Estimates of genetic and environmental variability in soybeans. Agronomy Journal, 47:314-318.

Robinson HF. (1966). Quantitative genetics in relation to breeding on centennial of Mendelism. Indian Journal of Genetics and Plant Breeding, 171.

Dudley JW and Moll RH. (1969). Interpretation and Use of estimates of heritability of genetics variance in Plant Breeding. Crop Science, 9:257-262.

Sivasubramaniam S and Madhavamenon P. (1973). Genotypic and Phenotypic variability in rice. The Madras Agricultural Journal, 60:1093-1096.
}
\author{
Nitesh, S.D., Parashuram Patroti and Shilpa Parashuram
}
\note{
Genetic parameter analysis need to be performed only if the sum of squares of treatment: Test is significant
Negative estimates of variance components is computed are not abnormal. Ref. Dudley and Moll (1969).
}

\examples{
data("Barley")
traitstats.genpar(Treatment=Barley$trt,
                  Replication=Barley$rep,
                  DataFile=Barley)

  }
