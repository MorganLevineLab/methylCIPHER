#' calcKnight
#'
#' @description A function to calculate the Knight gestational age
#'
#' @param DNAm a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param pheno Optional: The sample phenotype data (also with samples as rows) that the clock will be appended to.
#' @param CpGImputation An optional namesd vector for the mean value of each CpG that will be input from another dataset if such values are missing here (from sample cleaning)
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#'
#' @return If you added the optional pheno input (preferred) the function appends a column with the clock calculation and returns the dataframe. Otherwise, it will return a vector of calculated clock values in order of the
#' @export
#'
#' @examples calcKnight(exampleBetas, examplePheno, imputation = F)
calcKnight <- function(DNAm, pheno = NULL, CpGImputation = NULL, imputation = F){

  #######################
  ### Read in the Data###
  #######################

  #data("Knight_CpGs")

  ###################################################
  ### Check if all necessary CpGs are in the data ###
  ###################################################

  CpGCheck <- length(Knight_CpGs$CpG) == sum(Knight_CpGs$CpG %in% colnames(DNAm))

  ###################################################################################
  ### The calculation will be performed or an error will be thrown as appropriate ###
  ###################################################################################

  if(CpGCheck == F && is.null(CpGImputation) && imputation == T){

    stop("Need to provide of named vector of CpG Imputations; Necessary CpGs are missing!")

  } else if(CpGCheck == T | imputation == F){

    present <- Knight_CpGs$CpG %in% colnames(DNAm)

    betas <- DNAm[,na.omit(match(Knight_CpGs$CpG,colnames(DNAm)))]
    tt <- rowSums(sweep(as.matrix(betas) ,MARGIN = 2, Knight_CpGs$coef[present], `*`), na.rm = T) + 41.7

    if(is.null(pheno)){
      tt
    } else{
      pheno$Knight <- tt
      pheno
    }

  } else {
    message("Imputation of mean CpG Values occured for Hannum")
    missingCpGs <- Knight_CpGs$CpG[!(Knight_CpGs$CpG %in% colnames(DNAm))]
    tempDNAm <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs))

    for(j in 1:length(missingCpGs)){
      meanVals <- CpGImputation[match(missingCpGs[j],names(CpGImputation))]
      tempDNAm[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    colnames(tempDNAm) <- missingCpGs
    DNAm <- cbind(DNAm,tempDNAm)

    betas <- DNAm[,match(Knight_CpGs$CpG,colnames(DNAm))]
    tt <- rowSums(sweep(betas ,MARGIN = 2, Knight_CpGs$coef, `*`)) + 41.7

    if(is.null(pheno)){
      tt
    } else{
      pheno$Knight <- tt
      pheno
    }

  }

}
