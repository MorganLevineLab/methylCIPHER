#' @title Zhang 2019 ultra-precision clock elasticnet Cpgs
#'
#' @description The CpGs and Weights to calculate the Zhang 2019 elasticnet predictor. For the BLUP version, please see the original publication.
#'
#' @format A data frame with 514 rows (CpGs) and 2 variables:
#' \describe{
#'   \item{CpG}{The CpG names from Illumina Array IDs used in the current clock}
#'   \item{coef}{The weighted linear regression beta value of fit for each CpG in the clock.}
#' }
#' @source <https://doi.org/10.1186/s13073-019-0667-1>
"Zhang2019_CpGs"
