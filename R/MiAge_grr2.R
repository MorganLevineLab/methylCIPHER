#' MiAge grr2
#'
#' #' @description A function used by the original MiAge code found <http://www.columbia.edu/~sw2206/softwares.htm>
#'
#' @return The derivative of fr2 with respect to n_{j}
#' @export
MiAge_grr2 <- function(x,b,c,d,betaj)  ## derivative of fr.j with respect to n_{j}
{
  nj=x
  return(2*sum((c+b^(nj-1)*d-betaj)*b^(nj-1)*log(b)*d,na.rm=T))
}
