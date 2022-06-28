#' calcEpiTOC2
#'
#' @description A function to calculate the Epigenetic Time of Cancer 2 (EpiTOC2)
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
#' @examples calcEpiTOC2(exampleBetas, examplePheno, imputation = F)
calcEpiTOC2 <- function(DNAm, pheno = NULL, CpGImputation = NULL, imputation = T, approximated = F){

  #######################
  ### Read in the Data###
  #######################

  #data("EpiToc2_CpGs")

  ###################################################
  ### Check if all necessary CpGs are in the data ###
  ###################################################

  CpGCheck <- length(rownames(EpiToc2_CpGs)) == sum(rownames(EpiToc2_CpGs) %in% colnames(DNAm))

  ###################################################################################
  ### The calculation will be performed or an error will be thrown as appropriate ###
  ###################################################################################

  if(CpGCheck == F && is.null(CpGImputation) && imputation == T){

    stop("Need to provide of named vector of CpG Imputations; Necessary CpGs are missing!")

  } else if(CpGCheck == T | imputation == F){

    ### do epiTOC2
    map.idx <- match(rownames(EpiToc2_CpGs),colnames(DNAm));
    rep.idx <- which(is.na(map.idx)==FALSE);
    #print(paste("Number of represented epiTOC2 CpGs (max=163)=",length(rep.idx),sep=""))
    tmp.m <- as.matrix(DNAm[,map.idx[rep.idx]]);
    TNSC.v <- 2*colMeans(diag(1/(EpiToc2_CpGs[rep.idx,1]*(1-EpiToc2_CpGs[rep.idx,2])), nrow = length(rep.idx)) %*% (t(tmp.m) - EpiToc2_CpGs[rep.idx,2]),na.rm=TRUE);
    TNSC2.v <- 2*colMeans(diag(1/EpiToc2_CpGs[rep.idx,1], nrow = length(rep.idx)) %*% t(tmp.m),na.rm=TRUE);

    if(approximated == F){
      if(is.null(pheno)){
        TNSC.v
      } else{
        pheno$epiTOC2 <- TNSC.v
        pheno
      }
    } else if(approximated == T){
      if(is.null(pheno)){
        TNSC2.v
      } else{
        pheno$epiTOC2 <- TNSC2.v
        pheno
      }
    }


  } else {
    message("Imputation of mean CpG Values occured for epiTOC2")
    missingCpGs <- rownames(EpiToc2_CpGs)[!(rownames(EpiToc2_CpGs) %in% colnames(DNAm))]
    tempDNAm <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs))

    for(j in 1:length(missingCpGs)){
      meanVals <- CpGImputation[match(missingCpGs[j],names(CpGImputation))]
      tempDNAm[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    colnames(tempDNAm) <- missingCpGs
    DNAm <- cbind(DNAm,tempDNAm)

    ### do epiTOC2
    map.idx <- match(rownames(EpiToc2_CpGs),colnames(DNAm));
    rep.idx <- which(is.na(map.idx)==FALSE);
    #print(paste("Number of represented epiTOC2 CpGs (max=163)=",length(rep.idx),sep=""))
    tmp.m <- DNAm[,map.idx[rep.idx]];
    TNSC.v <- 2*colMeans(diag(1/(EpiToc2_CpGs[rep.idx,1]*(1-EpiToc2_CpGs[rep.idx,2]))) %*% (t(tmp.m) - EpiToc2_CpGs[rep.idx,2]),na.rm=TRUE);
    TNSC2.v <- 2*colMeans(diag(1/EpiToc2_CpGs[rep.idx,1]) %*% t(tmp.m),na.rm=TRUE);

    if(approximated == F){
      if(is.null(pheno)){
        TNSC.v
      } else{
        pheno$epiTOC2 <- TNSC.v
        pheno
      }
    } else if(approximated == T){
      if(is.null(pheno)){
        TNSC2.v
      } else{
        pheno$epiTOC2 <- TNSC2.v
        pheno
      }
    }
  }
}
