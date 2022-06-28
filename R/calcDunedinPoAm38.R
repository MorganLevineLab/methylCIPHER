#' calcDunedinPoAm38
#'
#' @description A function to calculate the Dunedin Pace of Aging Measure (DunedinPoAm38). This accesses the library available at github.com/danbelsky/DunedinPoAm38
#'
#' @param DNAm a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param pheno Optional: The sample phenotype data (also with samples as rows) that the clock will be appended to.
#' @param CpGImputation An optional namesd vector for the mean value of each CpG that will be input from another dataset if such values are missing here (from sample cleaning)
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#'
#' @return If you added the optional pheno input (preferred) the function appends a column with the clock calculation and returns the dataframe. Otherwise, it will return a vector of calculated clock values in order of the
#' @export
#' @source <https://doi.org/10.7554/eLife.54870>
#'
#' @examples calcDunedinPoAm38(exampleBetas, examplePheno, imputation = F)
calcDunedinPoAm38 <- function(DNAm, pheno = NULL, CpGImputation = NULL, imputation = T){

  ###########################
  ### Read in the Library ###
  ###########################

  library(DunedinPoAm38)
  DunedinPoAmCpGs <- DunedinPoAm38::getRequiredProbes()[[1]]

  ###################################################
  ### Check if all necessary CpGs are in the data ###
  ###################################################

  newBeta <- t(DNAm)
  CpGCheck <- length(DunedinPoAmCpGs) == sum(DunedinPoAmCpGs %in% colnames(newBeta))

  ###################################################################################
  ### The calculation will be performed or an error will be thrown as appropriate ###
  ###################################################################################

  if(CpGCheck == F && is.null(CpGImputation) && imputation == T){

    stop("Need to provide of named vector of CpG Imputations; Necessary CpGs are missing!")

  } else if(CpGCheck == T | imputation == F){

    ### do DunedinPoAm38

    PoAm <- unlist(DunedinPoAm38::PoAmProjector(newBeta))

    if(is.null(pheno)){
      PoAm
    } else{
      pheno$DunedinPoAm38 <- PoAm
      pheno
    }

  } else {
    message("Imputation of mean CpG Values occured for epiTOC2")
    missingCpGs <- DunedinPoAmCpGs[!(DunedinPoAmCpGs %in% rownames(newBeta))]
    tempDNAm <- matrix(ncol = dim(newBeta)[2], nrow = length(missingCpGs))

    for(j in 1:length(missingCpGs)){
      meanVals <- CpGImputation[match(missingCpGs[j],names(CpGImputation))]
      tempDNAm[j,] <- rep(meanVals,dim(newBeta)[2])
    }
    rownames(tempDNAm) <- missingCpGs
    newBeta <- cbind(newBeta,tempDNAm)

    ### do epiTOC
    PoAm <- unlist(DunedinPoAm38::PoAmProjector(newBeta))

    if(is.null(pheno)){
      PoAm
    } else{
      pheno$DunedinPoAm38 <- PoAm
      pheno
    }

  }

}
