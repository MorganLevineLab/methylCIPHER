#' calcHorvath2
#'
#' @description A function to calculate the Horvath Skin & Blood epigenetic clock
#'
#' @param DNAm a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param pheno Optional: The sample phenotype data (also with samples as rows) that the clock will be appended to.
#' @param CpGImputation An optional namesd vector for the mean value of each CpG that will be input from another dataset if such values are missing here (from sample cleaning)
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#'
#' @return If you added the optional pheno input (preferred) the function appends a column with the clock calculation and returns the dataframe. Otherwise, it will return a vector of calculated clock values in order of the
#' @export
#'
#' @examples calcHorvath2(exampleBetas, examplePheno, imputation = F)
calcHorvath2 <- function(DNAm, pheno = NULL, CpGImputation = NULL, imputation = F){

  #######################
  ### Read in the Data###
  #######################

  if(!exists("anti.trafo") || !exists("trafo")){
    stop("Calculation of Horvath2 requires that the functions trafo and anti.trafo is loaded. \n Ensure that you loaded the entire calcAllCpGClocks library")
  }

  #data("Horvath2_CpGs")

  ###################################################
  ### Check if all necessary CpGs are in the data ###
  ###################################################

  CpGCheck <- length(Horvath2_CpGs$ID) == sum(Horvath2_CpGs$ID %in% colnames(DNAm))

  ###################################################################################
  ### The calculation will be performed or an error will be thrown as appropriate ###
  ###################################################################################

  if(CpGCheck == F && is.null(CpGImputation) && imputation == T){

    stop("Need to provide of named vector of CpG Imputations; Necessary CpGs are missing!")

  } else if(CpGCheck == T | imputation == F){

    present <- Horvath2_CpGs$ID %in% colnames(DNAm)

    betas <- DNAm[,na.omit(match(Horvath2_CpGs$ID,colnames(DNAm)))]
    tt <- sweep(as.matrix(betas), MARGIN = 2, Horvath2_CpGs$Coef[present], `*`)

    Horvath2 <- as.numeric(anti.trafo(rowSums(tt,na.rm=T)-0.447119319))
    if(is.null(pheno)){
      Horvath2
    } else{
      pheno$Horvath2 <- Horvath2
      pheno
    }

  } else {
    message("Imputation of mean CpG Values occured for Horvath2")
    missingCpGs <- Horvath2_CpGs$ID[-1][!(Horvath2_CpGs$ID %in% colnames(DNAm))]
    tempDNAm <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs))

    for(j in 1:length(missingCpGs)){
      meanVals <- CpGImputation[match(missingCpGs[j],names(CpGImputation))]
      tempDNAm[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    colnames(tempDNAm) <- missingCpGs
    DNAm <- cbind(DNAm,tempDNAm)

    betas <- DNAm[,match(Horvath2_CpGs$ID[-1],colnames(DNAm))]
    tt <- sweep(betas, MARGIN = 2, Horvath2_CpGs$Coef, `*`)

    Horvath2 <- as.numeric(anti.trafo(rowSums(tt,na.rm=T)-0.447119319))
    if(is.null(pheno)){
      Horvath2
    } else{
      pheno$Horvath2 <- Horvath2
      pheno
    }

  }

}
