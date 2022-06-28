#' calcBocklandt
#'
#' @description A function to calculate the Bocklandt clock.
#'
#' @param DNAm a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param pheno Optional: The sample phenotype data (also with samples as rows) that the clock will be appended to.
#' @param CpGImputation An optional namesd vector for the mean value of each CpG that will be input from another dataset if such values are missing here (from sample cleaning)
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#'
#' @return If you added the optional pheno input (preferred) the function appends a column with the clock calculation and returns the dataframe. Otherwise, it will return a vector of calculated clock values in order of the
#' @export
#'
#' @examples calcBocklandt(exampleBetas, examplePheno, imputation = F)
calcBocklandt <- function(DNAm, pheno = NULL, CpGImputation = NULL, imputation = F){

  #######################
  ### Read in the Data###
  #######################

  #data("Bocklandt_CpG")

  ###################################################
  ### Check if all necessary CpGs are in the data ###
  ###################################################

  CpGCheck <- length(Bocklandt_CpG) == sum(Bocklandt_CpG %in% colnames(DNAm))

  ###################################################################################
  ### The calculation will be performed or an error will be thrown as appropriate ###
  ###################################################################################

  if(CpGCheck == F && imputation == T){

    stop("Necessary CpG is missing!")

  } else if(CpGCheck == T){

    Bocklandt=DNAm[,Bocklandt_CpG]

    if(is.null(pheno)){
      Bocklandt
    } else{
      pheno$Bocklandt <- Bocklandt
      pheno
    }

  }

}
