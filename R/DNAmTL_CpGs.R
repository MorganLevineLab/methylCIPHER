#' @title DNAmTL Weights
#'
#' @description The CpGs and Weights to calculate the DNAmTL Clock
#'
#' @format A data frame with 140 rows (CpGs) and 2 variables:
#' \describe{
#'   \item{ID}{The CpG names from Illumina Array IDs used in the current clock}
#'   \item{Coef}{The weighted linear regression beta value of fit for each CpG in the clock.}
#' }
#' @source <https://www.aging-us.com/article/102173/text>
"DNAmTL_CpGs"
