#' MiAge Mitotic Age
#'
#' @description A function used by the original MiAge code found <http://www.columbia.edu/~sw2206/softwares.htm>
#'
#' @return The MiAge estimation
#' @export
MiAge_mitotic.age <- function(beta,b,c,d) {
  library(methods)

  upperage=10000
  lowerage=10

  n=rep(500,ncol(beta)) #true parameter
  no.initial.n=5

  ### minimize the objective function "MiAge_fr2" for each patient in turn  #########
  for(j in 1:ncol(beta)) {
    current.value=MiAge_fr2(n[j],b,c,d,beta[,j])

    columnoptim=vector("list",no.initial.n)
    val=rep(NA,no.initial.n)
    for(jj in 1:(no.initial.n-1)) {
      temp=try(optim(par=lowerage+jj*(upperage-lowerage)/no.initial.n,fn=MiAge_fr2,gr=MiAge_grr2,b=b,c=c,d=d,betaj=beta[,j],method="L-BFGS-B",lower=lowerage,upper=upperage,control=list(factr=1)),silent=T)
      if(!is(temp,"try-error")){ columnoptim[[jj]]=temp; val[jj]=temp$value}
    }
    temp=try(optim(par=n[j],fn=MiAge_fr2,gr=MiAge_grr2,b=b,c=c,d=d,betaj=beta[,j],method="L-BFGS-B",lower=lowerage,upper=upperage,control=list(factr=1)),silent=T)
    if(!is(temp,"try-error")) { columnoptim[[no.initial.n]]=temp; val[no.initial.n]=temp$value}
    temp=columnoptim[[ which(val==min(val,na.rm=T))[1] ]]


    if(!is(temp,"try-error"))
    {
      n[j]=temp$par
    } else {print(2);print(temp)}

  }

  return(n)

}
