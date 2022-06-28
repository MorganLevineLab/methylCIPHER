#' @title Lin Weights
#'
#' @description The CpGs and Weights to calculate the Lin Clock
#'
#' @format A data frame with 99 rows (CpGs) and 2 variables:
#' \describe{
#'   \item{id}{The CpG names from Illumina Array IDs used in the current clock}
#'   \item{coef}{The weighted linear regression beta value of fit for each CpG in the clock.}
#' }
#' @source <https://doi.org/10.18632/aging.100908>
"Lin_CpGs"
