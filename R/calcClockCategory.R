#' calcClockCategory
#'
#' @description A function to calculate all clocks in a chose category
#'
#' @param DNAm a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param pheno The sample phenotype data (also with samples as rows) that the clock will be appended to. If you don't have this, simply initialize a dataframe with a single column of sample IDs or rownumbers you can append clock values to.
#' @param CpGImputation An optional namesd vector for the mean value of each CpG that will be input from another dataset if such values are missing here (from sample cleaning)
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#'
#' @section Warning! If you use this function, look for agreement between category functions, and if you use many or all, it is essential that you perform multiple testing correction.
#'
#' @return A dataframe that has column names of the core clocks giving you multiple clocks at once in an easy to compute function. These will be appended onto the existing pheno dataframe as defined in the inputs.
#' @export
#'
#' @examples calcClockCategory(exampleBetas, examplePheno, category = "chronological", imputation = F)
calcClockCategory <- function(DNAm, pheno , category = NULL,
                              CpGImputation = NULL, imputation = F){

  message("Please remember to cite all of the clocks you have used! Please refer to the README.md file for assistance.")

  if(suppressWarnings(as.logical("originals" %in% category))){

    pheno <- calcBocklandt(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcGaragnani(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcWeidner(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcVidalBralo(DNAm, pheno, CpGImputation, imputation)

  }
  if(suppressWarnings(as.logical("chronological" %in% category))){

    pheno <- calcHannum(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcHorvath1(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcLin(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcZhang2019(DNAm, pheno, CpGImputation, imputation)

  }
  if(suppressWarnings(as.logical("phenotypic" %in% category))){

    pheno <- calcPhenoAge(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcHRSInChPhenoAge(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcDunedinPoAm38(DNAm, pheno, CpGImputation, imputation)
    #pheno <- prcPhenoAge::calcPRCPhenoAge(DNAm, pheno, CpGImputation, imputation)
    #pheno <- prcPhenoAge::calcnonPRCPhenoAge(DNAm, pheno, CpGImputation, imputation)

  }
  if(suppressWarnings(as.logical("mitotic" %in% category))){

    pheno <- calcEpiTOC(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcMiAge(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcDNAmTL(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcEpiTOC2(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcHypoClock(DNAm, pheno, CpGImputation, imputation)

  }
  if(suppressWarnings(any(as.logical(category %in% c("pediatric", "gestational"))))){

    pheno <- calcBohlin(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcKnight(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcLeeControl(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcLeeRobust(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcLeeRefinedRobust(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcMayne(DNAm, pheno, CpGImputation, imputation)
    pheno <- calcPEDBE(DNAm, pheno, CpGImputation, imputation)

  }

  if(suppressWarnings(as.logical("adult" %in% category))){

    clockNames <- base::setdiff(clockOptions()[1:30],
                                c("calcBohlin","calcKnight","calcLeeControl",
                                  "calcLeeRobust","calcLeeRefinedRobust","calcMayne",
                                  "calcPEDBE"))
    for(i in 1:length(clockNames)){
      clockFun <- get(clockNames[i])
      pheno <- clockFun(DNAm, pheno, CpGImputation, imputation)
    }
  }

  if(suppressWarnings(as.logical("all" %in% category))){

    clockNames <- clockOptions()[1:30]
    for(i in 1:length(clockNames)){
      clockFun <- get(clockNames[i])
      pheno <- clockFun(DNAm, pheno, CpGImputation, imputation)
    }

  }

  pheno

}
