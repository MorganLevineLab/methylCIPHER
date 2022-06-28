#' MiAge fr2
#'
#' #' @description A function used by the original MiAge code found <http://www.columbia.edu/~sw2206/softwares.htm>
#'
#' @return objective fuction summed over all CpG istes for patient j
#' @export
MiAge_fr2 <- function(x,b,c,d,betaj)  ## objective fuction summed over all CpG istes for patient j
{
  nj=x
  return(sum((c+b^(nj-1)*d-betaj)^2,na.rm=T))
}
