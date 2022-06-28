#' @title Knight (Gestational Age) Predictor Weights
#'
#' @description The CpGs and Weights to calculate the Knight Gestational Age
#'
#' @format A data frame with 148 rows (CpGs) and 2 variables:
#' \describe{
#'   \item{CpG}{The CpG names from Illumina Array IDs used in the current clock}
#'   \item{coef}{The weighted linear regression beta value of fit for each CpG in the clock.}
#' }
#' @source <https://doi.org/10.1186/s13059-016-1068-z>
"Knight_CpGs"
