#' clockOptions
#'
#' @return A list of the available clocks to calculate, including those in the prcPhenoAge and DunedinPoAm packages
#' @export
#'
clockOptions <- function() {
   x <- unlist(strsplit(lsf.str("package:methylCIPHER", pattern = "^calc"),"[:]"))
   y <- setdiff(x, c("calcUserClocks","calcCoreClocks"))
   z <- union(y, c("prcPhenoAge::calcPRCPhenoAge","prcPhenoAge::calcnonPRCPhenoAge","DunedinPoAm38::PoAmProjector"))

   z
}
