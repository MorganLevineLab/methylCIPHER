#' Removal of empty CpG columns
#'
#' @param x A dataframe of CpG Betas with CpGs of all NA values
#'
#' @return Dataframe with NA columns removed
#' @export
removeNAcol <- function(x){

  removeCols <- apply(x,2,function(z)all(is.na(z)))
  return(x[,!removeCols])

}
