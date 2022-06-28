#' calcUserClocks
#'
#' @description A function to calculate multiple epigenetic clocks at once, as defined by the user.
#'
#' @param clockList A character vector of clock calculation functions the user would like to compute at once. Options can be found by running the function clockOptions()
#' @param DNAm a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param pheno The sample phenotype data (also with samples as rows) that the clock will be appended to. If you don't have this, simply initialize a dataframe with a single column of sample IDs or rownumbers you can append clock values to.
#' @param CpGImputation An optional namesd vector for the mean value of each CpG that will be input from another dataset if such values are missing here (from sample cleaning)
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#'
#' @return A dataframe that has column names of the core clocks giving you multiple clocks at once in an easy to compute function. These will be appended onto the existing pheno dataframe as defined in the inputs.
#' @export
#'
#' @examples userClocks <- c("calcSmokingMcCartney","calcPhenoAge","calcEpiTOC2")
#' @examples calcUserClocks(userClocks, exampleBetas, examplePheno, imputation = F)
calcUserClocks <- function(clockList, DNAm, pheno , CpGImputation = NULL, imputation = F){

  for(i in 1:length(clockList)){

    x <- get(clockList[i])
    pheno <- x(DNAm, pheno, CpGImputation, imputation)

  }

  pheno

}
