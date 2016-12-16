#'Sort out sampling times, coalescent times and sampling lineages from a phylogenetic tree
#'
##' \code{att} sorts out sampling times, coalescent times and sampling lineages from a phylogenetic tree.
##'
##' @param phy A phylogenetic tree.
##' @param eps Difference parameter to separate coalescent and sampling event.
##'
##' @return Sorted sampling times, coalescent times and sampling lineages.
##'
##' @references Palacios JA and Minin VN. Integrated nested Laplace approximation for Bayesian nonparametric phylodynamics, in Proceedings of the Twenty-Eighth Conference on Uncertainty in Artificial Intelligence, 2012.
##' @examples
##' library(ape)
##' t1=rcoal(20)
##' att(t1)
##'
##'
##' @author Simon Frost (\email{sdwfrost@@gmail.com})
##'
##'
##' @export


att <- function(phy,eps=1e-6){
  b.s.times = branching.sampling.times(phy)
  int.ind = which(as.numeric(names(b.s.times)) < 0)
  tip.ind = which(as.numeric(names(b.s.times)) > 0)
  num.tips = length(tip.ind)
  num.coal.events = length(int.ind)
  sampl.suf.stat = rep(NA, num.coal.events)
  coal.interval = rep(NA, num.coal.events)
  coal.lineages = rep(NA, num.coal.events)
  sorted.coal.times = sort(b.s.times[int.ind])
  names(sorted.coal.times) = NULL
  sampling.times = sort((b.s.times[tip.ind]))
  for (i in 2:length(sampling.times)){
    if ((sampling.times[i]-sampling.times[i-1])<eps){
      sampling.times[i]<-sampling.times[i-1]}
  }
  unique.sampling.times<-unique(sampling.times)
  sampled.lineages = NULL
  for (sample.time in unique.sampling.times){
    sampled.lineages = c(sampled.lineages,
                         sum(sampling.times == sample.time))
  }
  if(sum(sorted.coal.times %in% unique.sampling.times)>0){
    stop("Simultaneous sample and coalescence time")
  }
  all.times <- sort(unique(c(unique.sampling.times,sorted.coal.times)))
  # Check that first time is sampling
  if(!(all.times[1] %in% unique.sampling.times)){
    stop("Samples must be first (in reverse time)")
  }
  A <- rep(0,length(all.times))
  lastA <- 0
  for(i in 1:length(all.times)){
    is.sample <- match(all.times[i],unique.sampling.times)
    if(!is.na(is.sample)){
      ss <- sampled.lineages[is.sample]
      A[i] <- lastA + ss
    }else{
      A[i] <- lastA - 1
    }
    lastA <- A[i]
  }
  data.frame(t=all.times,A=A)
}
