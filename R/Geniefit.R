#' A function to fit coalescent models to a given phylogenetic tree.
#' @param phy A phylogenetic tree.
#' @param Model A Model choice from const (constant population size), expo (exponetial growth),expan (expansion growth), log (logistic growth), step (piecewise constant), pexpan (piecewise expansion growth) and plog (piecewise logistic growth).
#' @param start Initial values for the parameters to be optimized over.
#' @param lower, upper Bounds on the variables.
#' @return Parameters estimation of a given model, loglikelihood and AIC
#' @examples
#' library(ape)
#' t1=rcoal(20)
#' Geniefit(t1,Model="expo",start=c(100,.1,.1),upper=Inf,lower=0)
#'
#' @author Fei Xiang (\email{xf3087@@gmail.com})
#'
#' @export
#'
#'
#'
#'
#######one function to produce the fit######
Geniefit=function(phy,Model="user",start,upper,lower){
  #####wash the data from the tree file#########
  phy.times=heterochronous.gp.stat (phy)
  ##################times frame given the coalesent events#############
  n=length(phy.times$coal.times)
  coaltimes.pop=matrix(0,nrow=n,ncol=2)
  coaltimes.pop[,1]=phy.times$coal.times
  coaltimes.pop[,2]=-1
  #################times frame given the sampled events################
  xn=length(phy.times$sample.times)
  samptimes.pop=matrix(0,nrow=xn,ncol=2)
  samptimes.pop[,1]=phy.times$sample.times
  samptimes.pop[,2]=phy.times$sampled.lineages
  ######sorted time and index matrix#####
  times.pop=rbind(samptimes.pop,coaltimes.pop)
  sort.times.pop=times.pop
  sort.times.pop[,1]=times.pop[,1][order(times.pop[,1])]
  sort.times.pop[,2]=times.pop[,2][order(times.pop[,1])]
  #####population at diffrent times###
  pop.times=cumsum(sort.times.pop[,2])
  #####type of time, 0 denoting sampling event and 1 denoting coalesent event####
  type=c(rep(0,xn),rep(1,n))
  sort.type=type[order(times.pop[,1])]
  ntotal=length(sort.type)
  #####if statement to get rid of first event when it is sampling event##########
  if (pop.times[1]<2) {
    pop.times=pop.times[-1]
    sort.times.pop=sort.times.pop[-1,]
    sort.times.pop[,1]=sort.times.pop[,1]-sort.times.pop[1,1]
    ntotal=ntotal-1
    sort.type=sort.type[-1]
  }
  ######population trajectory function#########
  fnpar=function(parr){
    ####function of t for population trajectory#####
    fnt=function(t){
      if (Model=="const") {trajectory=parr[1]}
      if (Model=="expo")  {trajectory=parr[1]*exp(-parr[2]*t)}
      if (Model=="expan") {trajectory=parr[1]*(parr[3]+(1-parr[3])*exp(-parr[2]*t))}
      if (Model=="log")   {trajectory=parr[1]*((1+parr[3])/(1+parr[3]*exp(parr[2]*t)))}
      if (Model=="step")  {trajectory=ifelse(t<parr[3],parr[1],parr[1]*parr[2])}
      if (Model=="pexpan") {trajectory=ifelse(t<-log(parr[3])/parr[2],parr[1]*exp(-parr[2]*t),parr[1]*parr[3])}
      return(1/trajectory)
    }
    ######define the integral explicit function given fnt###########
    intfnt=function(lowerlim,upperlim){
      if (Model=="const") {intg=1/parr[1]*(upperlim-lowerlim)}
      if (Model=="expo")  {intg=1/parr[1]/parr[2]*(exp(parr[2]*upperlim)-exp(parr[2]*lowerlim))}
      if (Model=="expan") {intg=1/parr[1]/parr[2]/parr[3]*(log(parr[3]*exp(parr[2]*upperlim)+1-parr[3])-log(parr[3]*exp(parr[2]*lowerlim)+1-parr[3]))}
      if (Model=="log")   {intg=1/parr[1]/(1+parr[3])*(upperlim-lowerlim+parr[3]/parr[2]*(exp(parr[2]*upperlim)-exp(parr[2]*lowerlim)))}
      if (Model=="step")  {
        intg=ifelse(upperlim<parr[3],1/parr[1]*(upperlim-lowerlim),ifelse(lowerlim>parr[3],1/parr[1]/parr[2]*(upperlim-lowerlim),1/parr[1]*(parr[3]-lowerlim)+1/parr[1]/parr[2]*(upperlim-parr[3])) )
      }
      if (Model=="pexpan") {
        intg=(upperlim< -log(parr[3])/parr[2])*1/parr[1]/parr[2]*(exp(parr[2]*upperlim)-exp(parr[2]*lowerlim))+(lowerlim>-log(parr[3])/parr[2])*1/parr[1]/parr[3]*(upperlim-lowerlim)+(lowerlim< -log(parr[3])/parr[2] && -log(parr[3])/parr[2]<upperlim)*(1/parr[1]/parr[3]*(upperlim+log(parr[3])/parr[2])+1/parr[1]/parr[2]*(exp(parr[2]*
                                                                                                                                                                                                                                                                                                                                      -log(parr[3])/parr[2])-exp(parr[2]*lowerlim)))
      }
      return(intg)
    }
    logcoe=sort.type[-1]*(log(fnt(sort.times.pop[-1,1]))+log(pop.times[-ntotal]*(pop.times[-ntotal]-1)/2))
    logint=-pop.times[-ntotal]*(pop.times[-ntotal]-1)/2*intfnt(sort.times.pop[-ntotal,1],sort.times.pop[-1,1])
    return(-sum(logcoe)-sum(logint))

  }
  fn2 <- function(x){
    fnpar(exp(x))
  }
  require(minqa)
  fit2 <- bobyqa(log(start),fn2,lower=log(lower),upper=log(upper))
  return(list(parr=exp(fit2$par),loglikelihood=-fit2$fval,AIC=2*length(start)+2*fit2$fval))
}

