#' calcLeeRefinedRobust
#'
#' @description A function to calculate the Refined version of the Robust Predictor of Placental Age (Lee)
#'
#' @param DNAm a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param pheno Optional: The sample phenotype data (also with samples as rows) that the clock will be appended to.
#' @param CpGImputation An optional namesd vector for the mean value of each CpG that will be input from another dataset if such values are missing here (from sample cleaning)
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#'
#' @return If you added the optional pheno input (preferred) the function appends a column with the clock calculation and returns the dataframe. Otherwise, it will return a vector of calculated clock values in order of the
#' @export
#'
#' @examples calcLeeRefinedRobust(exampleBetas, examplePheno, imputation = F)
calcLeeRefinedRobust <- function(DNAm, pheno = NULL, CpGImputation = NULL, imputation = F){

  #######################
  ### Read in the Data###
  #######################

  #data("LeeRefinedRobust_CpGs")

  ###################################################
  ### Check if all necessary CpGs are in the data ###
  ###################################################
  CpGCheck <- length(LeeRefinedRobust_CpGs$CpG) == sum(LeeRefinedRobust_CpGs$CpG %in% colnames(DNAm))

  ###################################################################################
  ### The calculation will be performed or an error will be thrown as appropriate ###
  ###################################################################################

  if(CpGCheck == F && is.null(CpGImputation) && imputation == T){

    stop("Need to provide of named vector of CpG Imputations; Necessary CpGs are missing!")

  } else if(CpGCheck == T | imputation == F){

    present <- LeeRefinedRobust_CpGs$CpG %in% colnames(DNAm)

    betas <- DNAm[,na.omit(match(LeeRefinedRobust_CpGs$CpG,colnames(DNAm)))]
    tt <- rep(0, dim(DNAm)[1])

    tt <- rowSums(sweep(as.matrix(betas), MARGIN = 2, LeeRefinedRobust_CpGs$coef[present],`*`), na.rm = T) + 30.74966

    if(is.null(pheno)){
      tt
    } else{
      pheno$LeeRefinedRobust <- tt
      pheno
    }

  } else {
    message("Imputation of mean CpG Values occured for Lee Robust Refined")
    missingCpGs <- LeeRefinedRobust_CpGs$CpG[!(LeeRefinedRobust_CpGs$CpG %in% colnames(DNAm))]
    tempDNAm <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs))

    for(j in 1:length(missingCpGs)){
      meanVals <- CpGImputation[match(missingCpGs[j],names(CpGImputation))]
      tempDNAm[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    colnames(tempDNAm) <- missingCpGs
    DNAm <- cbind(DNAm,tempDNAm)

    betas <- DNAm[,match(LeeRefinedRobust_CpGs$CpG,colnames(DNAm))]
    tt <- rep(0, dim(DNAm)[1])
    tt <- rowSums(sweep(betas, MARGIN = 2, LeeRefinedRobust_CpGs$coef,`*`)) + 30.74966

    if(is.null(pheno)){
      tt
    } else{
      pheno$LeeRefinedRobust <- tt
      pheno
    }

  }

}
