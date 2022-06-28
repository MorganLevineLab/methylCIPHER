#' subsetCG
#'
#' @description A function to quickly subset the CpGs you want to work with in a methylation dataframe
#'
#' @param dat The methylation Beta values you will need to subset, where columns are CpGs (and are named), and rows are samples.
#' @param cgSet The character vector of Illumina CpG probe IDs that you will be subsetting to.
#'
#' @return A new beta methylation matrix with columns of only the CpGs that you wanted to subset to.
#' @export
#'
#' @examples subsetCG(exampleBetas, HorvathOnlineRef$Name)
subsetCG <- function(dat,cgSet){
  match1 = match(cgSet,colnames(dat))
  if(any(is.na(match1))){warning(paste(sum(is.na(match1)), "CpGs in the requested subset are missing from your data."))}
  match1 <- match1[!is.na(match1)]
  datReduced = dat[,match1]
  datReduced
}
