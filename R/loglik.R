require(ape)
Bayesfnt=function(t,par,Model){
  ##the function of a population trajectory given a model#####
  ##t is time units; par is parameters####
  if (Model=="const") {trajectory=par[1]}
  if (Model=="expo")  {trajectory=par[1]*exp(-par[2]*t)}
  if (Model=="expan") {trajectory=par[1]*(par[3]+(1-par[3])*exp(-par[2]*t))}
  if (Model=="log")   {trajectory=par[1]*((1+par[3])/(1+par[3]*exp(par[2]*t)))}
  if (Model=="step")  {trajectory=ifelse(t<=par[3],par[1],par[1]*par[2])}
  if (Model=="pexpan") {trajectory=ifelse(t< -log(par[3])/par[2],par[1]*exp(-par[2]*t),par[1]*par[3])}
  if (Model=="plog") {trajectory=ifelse(t<par[3],par[1],par[1]*exp(-par[2]*(t-par[3])))}
  return(1/trajectory)
}

Bayesintfnt=function(lowerlim,upperlim,parr,Model){
  if (Model=="const") {intg=1/parr[1]*(upperlim-lowerlim)}
  if (Model=="expo")  {intg=1/parr[1]/parr[2]*(exp(parr[2]*upperlim)-exp(parr[2]*lowerlim))}
  if (Model=="expan") {intg=1/parr[1]/parr[2]/parr[3]*(log(parr[3]*exp(parr[2]*upperlim)+1-parr[3])-log(parr[3]*exp(parr[2]*lowerlim)+1-parr[3]))}
  if (Model=="log")   {intg=1/parr[1]/(1+parr[3])*(upperlim-lowerlim+parr[3]/parr[2]*(exp(parr[2]*upperlim)-exp(parr[2]*lowerlim)))}
  if (Model=="step")  {
    intg=ifelse(upperlim<parr[3],1/parr[1]*(upperlim-lowerlim),ifelse(lowerlim>parr[3],1/parr[1]/parr[2]*(upperlim-lowerlim),1/parr[1]*(parr[3]-lowerlim)+1/parr[1]/parr[2]*(upperlim-
                                                                                                                                                                               parr[3])) )
  }
  if (Model=="pexpan") {
    intg=ifelse(upperlim< -log(parr[3])/parr[2], (exp(parr[2]*upperlim)-exp(parr[2]*lowerlim))/parr[1]/parr[2],ifelse(lowerlim>-log(parr[3])/parr[2], (upperlim-lowerlim)/parr[1]/parr[3],(upperlim+log(parr[3])/parr[2])/parr[1]/parr[3]+(1/parr[3]-exp(parr[2]*lowerlim))/parr[1]/parr[2]))
  }
  if (Model=="plog"){
    intg=(upperlim<parr[3])*(upperlim-lowerlim)/parr[1]+(lowerlim>parr[3])*(-exp(parr[2]*lowerlim-parr[2]*parr[3])+exp(parr[2]*upperlim-parr[2]*parr[3]))/parr[1]/parr[2]+(lowerlim<parr[3] && upperlim>parr[3])*((upperlim-parr[3])/parr[1]+(-exp(parr[2]*lowerlim-parr[2]*parr[3])+1)/parr[1]/parr[2])
  }
  return(intg)
}



#########the loglikelihoood function########
fntreeloglik=function(phy,para,Model){
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
  #if (pop.times[1]<2) {
  #pop.times=pop.times[-1]
  #sort.times.pop=sort.times.pop[-1,]
  #sort.times.pop[,1]=sort.times.pop[,1]-sort.times.pop[1,1]
  #ntotal=ntotal-1
  #sort.type=sort.type[-1]
  #}

  logcoe=sort.type[-1]*(log(Bayesfnt(sort.times.pop[-1,1],para,Model))+log(pop.times[-ntotal]*(pop.times[-ntotal]-1)/2))
  logcoe[which(sort.type[-1]==0)]=0
  logint=-pop.times[-ntotal]*(pop.times[-ntotal]-1)/2*Bayesintfnt(sort.times.pop[-ntotal,1],sort.times.pop[-1,1],para,Model)





  loglik=sum(logcoe+logint)
  return(-loglik)
}
#' A function to calculate the log likelihood of a phylogenetic tree given various coalescent models.
#' @param phy A phylogenetic tree.
#' @param Model A Model choice from const (constant population size), expo (exponetial growth),expan (expansion growth), log (logistic growth), step (piecewise constant), pexpan (piecewise expansion growth) and plog (piecewise logistic growth).
#' @param parr Parameters.
#' @param Rcpp Calculation is based on C++ code when it is True and on R code when it is False.
#' @return log likelihood of a phylogeny given a coalescent model and parameters.
#' @examples
#' data(village)
#' system.time(loglik(village,c(100,1),Model="expo",Rcpp=T))
#' system.time(loglik(village,c(100,1),Model="expo",Rcpp=F))
#' @author Fei Xiang (\email{xf3087@@gmail.com})
#'
#' @export
#'
#'
#'
#'

loglik=function(phy,parr,Model,Rcpp=F){
  if (Rcpp==F) return(-fntreeloglik(phy,para=parr,Model))
  if (Rcpp==T) {
    x=att(phy)
    if (Model=="const")  return(-negllc_const(parr,t=x$t,A=x$A))
    if (Model=="expo")   return(-negllc_expo(parr,t=x$t,A=x$A))
    if (Model=="pexpan") return(-negllc_pexpan(parr,t=x$t,A=x$A))
    if (Model=="log")    return(-negllc_log(parr,t=x$t,A=x$A))
    if (Model=="step")   return(-negllc_step(parr,t=x$t,A=x$A))
    if (Model=="expan")  return(-negllc_expan(parr,t=x$t,A=x$A))
    if (Model=="plog")   return(-negllc_plog(parr,t=x$t,A=x$A))
  }

}



