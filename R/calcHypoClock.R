#' calcHypoClock
#'
#' @description A function to calculate the HypoClock
#'
#' @param DNAm a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param pheno Optional: The sample phenotype data (also with samples as rows) that the clock will be appended to.
#' @param CpGImputation An optional namesd vector for the mean value of each CpG that will be input from another dataset if such values are missing here (from sample cleaning)
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#'
#' @return If you added the optional pheno input (preferred) the function appends a column with the clock calculation and returns the dataframe. Otherwise, it will return a vector of calculated clock values in order of the
#' @export
#' @source <https://zenodo.org/record/2632938#.YfGA3S-B2Cg>
#'
#' @examples calcHypoClock(exampleBetas, examplePheno, imputation = F)
calcHypoClock <- function(DNAm, pheno = NULL, CpGImputation = NULL, imputation = F){

  #######################
  ### Read in the Data###
  #######################

  #data("hypoClock_CpGs")

  ###################################################
  ### Check if all necessary CpGs are in the data ###
  ###################################################

  CpGCheck <- length(hypoClock_CpGs) == sum(hypoClock_CpGs %in% colnames(DNAm))

  ###################################################################################
  ### The calculation will be performed or an error will be thrown as appropriate ###
  ###################################################################################

  if(CpGCheck == F && is.null(CpGImputation) && imputation == T){

    stop("Need to provide of named vector of CpG Imputations; Necessary CpGs are missing!")

  } else if(CpGCheck == T | imputation == F){

    ### do HypoClock
    common.v <- intersect(colnames(DNAm),hypoClock_CpGs);
    #print(paste("Number of represented solo-WCGWs (max=678)=",length(common.v),sep=""));
    map.idx <- match(common.v,colnames(DNAm));
    hypoSC.v <- rowMeans(DNAm[,map.idx],na.rm=TRUE);

    if(is.null(pheno)){
      hypoSC.v
    } else{
      pheno$hypoClock <- hypoSC.v
      pheno
    }

  } else {
    message("Imputation of mean CpG Values occured for epiTOC2")
    missingCpGs <- hypoClock_CpGs[!(hypoClock_CpGs %in% colnames(DNAm))]
    tempDNAm <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs))

    for(j in 1:length(missingCpGs)){
      meanVals <- CpGImputation[match(missingCpGs[j],names(CpGImputation))]
      tempDNAm[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    colnames(tempDNAm) <- missingCpGs
    DNAm <- cbind(DNAm,tempDNAm)

    ### do HypoClock
    common.v <- intersect(colnames(DNAm),hypoClock_CpGs);
    #print(paste("Number of represented solo-WCGWs (max=678)=",length(common.v),sep=""));
    map.idx <- match(common.v,colnames(DNAm));
    hypoSC.v <- rowMeans(DNAm[,map.idx],na.rm=TRUE);

    if(is.null(pheno)){
      hypoSC.v
    } else{
      pheno$hypoClock <- hypoSC.v
      pheno
    }

  }
}
