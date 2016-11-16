#' Simulate coalescent times for heterochronous data.
#'
#' \code{coalgen_hetero} simulates coalescent times for heterochronous data.
#'
#' @param sample A two columns matrix of number of individuals and the initial time.
#' @param trajectory A population growth function.
#' @param val_upper Upper end of time points to be simulated.
#'
#'
#' @return Coalescent intervals and lineages.
#'
#' @references \url{https://github.com/JuliaPalacios/coalsieve}.
#'
#'
#'
#' @examples
#' sample1<-cbind(c(9,1,2,1),c(0,.008,.03,.1))
#'
#' trajectory<-function(x)  exp(10*x)
#' example_hetero<-coalgen_hetero(sample1, trajectory)
#'@author Fei Xiang (\email{xf3087@@gmail.com})
#'
#'@export
#'
#'
#'
coalgen_hetero <-function(sample, trajectory,val_upper=10){
  # sample = is a matrix with 2 columns. The first column contains the number of samples collected at the time defined in the second column
  # trajectory = one over the effective population size function
  # this works for heterochronous sampling
  # assumes sample[1,1]>1
  s=sample[1,2]
  b<-sample[1,1]
  n<-sum(sample[,1])-1
  m<-n
  nsample<-nrow(sample)
  sample<-rbind(sample,c(0,10))
  out<-rep(0,n)
  branches<-rep(0,n)
  i<-1
  while (i<(nsample+1)){
    if (b==1) {break}
    if (b<2){
      b<-b+sample[i+1,1]
      s<-sample[i+1,2]
      i<-i+1
    }
    x<-rexp(1)
    f <- function(bran,u,x,s) .5*bran*(bran-1)*integrate(trajectory, s, s+u)$value - x
    y<-uniroot(f,bran=b,x=x,s=s,lower=0,upper=val_upper)$root
    while ( (s+y)>sample[i+1,2]) {
      #     f <- function(bran,u,x,s) .5*bran*(bran-1)*integrate(trajectory, s, s+u)$value - x
      #     y<-uniroot(f,bran=b,x=x,s=s,lower=0,upper=val_upper)$root
      x<-x-.5*b*(b-1)*integrate(trajectory,s,sample[i+1,2])$value
      b<-b+sample[i+1,1]
      s<-sample[i+1,2]
      i<-i+1
      f <- function(bran,u,x,s) .5*bran*(bran-1)*integrate(trajectory, s, s+u)$value - x
      y<-uniroot(f,bran=b,x=x,s=s,lower=0,upper=val_upper)$root
      if (i==nsample) {sample[nsample+1,2]<-10*(s+y)}
    }

    s<-s+y
    out[m-n+1]<-s
    branches[m-n+1]<-b
    n<-n-1
    b<-b-1
    if (i==nsample) {sample[nsample+1,2]<-10*(s+y)}

  }

  return(list(branches=c(out[1],diff(out)),lineages=branches))
}
