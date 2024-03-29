% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Zhang2019_CpGs.R
\docType{data}
\name{Zhang2019_CpGs}
\alias{Zhang2019_CpGs}
\title{Zhang 2019 ultra-precision clock elasticnet Cpgs}
\format{
A data frame with 514 rows (CpGs) and 2 variables:
\describe{
  \item{CpG}{The CpG names from Illumina Array IDs used in the current clock}
  \item{coef}{The weighted linear regression beta value of fit for each CpG in the clock.}
}
}
\source{
<https://doi.org/10.1186/s13073-019-0667-1>
}
\usage{
Zhang2019_CpGs
}
\description{
The CpGs and Weights to calculate the Zhang 2019 elasticnet predictor. For the BLUP version, please see the original publication.
}
\keyword{datasets}
