#' trafo
#'
#' @param x A vector of sample ages
#' @param adult.age Age boundary of adulthood--set to 20
#'
#' @return a vector of transformed ages
#' @export
trafo= function(x,adult.age= 20) {
  x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y
  }
