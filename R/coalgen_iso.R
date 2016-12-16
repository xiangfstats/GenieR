#' Simulate coalescent times for isochronous data.
#'
#' \code{coalgen_iso} simulates coalescent times for isochronous data.
#'
#' @param sample A two dimensional vector of number of individuals and the initial time.
#' @param trajectory A population growth function.
#' @param val_upper Upper end of time points to be simulated.
#'
#'
#' @return Coalescent intervals and lineages.
#'
#' @references \url{https://github.com/JuliaPalacios/coalsieve}.
#'
#' @export
#'
#' @examples
#' sample<-c(100,0)
#'
#' trajectory<-function(x)  exp(10*x)
#' example_iso<-coalgen_iso(sample, trajectory)
#'
#'@author Fei Xiang (\email{xf3087@@gmail.com})
#'
#'
coalgen_iso<-function(sample, trajectory,val_upper=10){
  # sample = is a matrix with 2 columns. The first column contains the number of samples collected at the time defined in the second column
  # trajectory = one over the effective population size function
  # this works for isochronous sampling
  s=sample[2]
  n<-sample[1]
  out<-rep(0,n-1)
  #   val_upper<-10*(1/choose(n+1,3))
  for (j in n:2){
    t=rexp(1,choose(j,2))
    # trajectory is the inverse of the effective population size function
    f <- function(x,t,s) integrate(trajectory, s, s+x)$value - t
    #--- I will probably need to catch an error here, for val_upper, it breaks if
    # val_upper is not large enough
    #    val_upper<-10
    y<-uniroot(f,t=t,s=s,lower=0,upper=val_upper)$root
    s<-s+y
    out[n-j+1]<-s
  }
  return(list(branches=c(out[1],diff(out)),lineages=seq(n,2,-1)))

}

