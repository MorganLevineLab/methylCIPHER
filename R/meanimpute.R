#' Mean Imputation for singly missing CpGs
#'
#' @param x A dataframe of CpG Betas with missing values
#'
#' @return Mean Imputed dataframe
#' @export
meanimpute <- function(x){
  apply(x,2,function(z)ifelse(is.na(z),mean(z,na.rm=T),z))
}
