#' A function to fit coalescent models to a given phylogenetic tree.
#' @param phy A phylogenetic tree.
#' @param Model A Model choice from const (constant population size), expo (exponetial growth),expan (expansion growth), log (logistic growth), step (piecewise constant), pexpan (piecewise expansion growth) and plog (piecewise logistic growth).
#' @param start Initial values for the parameters to be optimized over.
#' @param lower, upper Bounds on the variables.
#' @param Rcpp Calculation is based on C++ code when it is True and on R code when it is False.
#' @param MCMC MCMC simulation is run when it is true, and not run when it is False. The default prior is uniform given the lower and upper.
#' @param sig MCMC simulation step size.
#' @param run Number of MCMC simulation.
#' @return Parameters estimation of a given model, loglikelihood and AIC
#' @examples
#' library(ape)
#' t1=rcoal(20)
#' Geniefit(t1,Model="expo",start=c(100,.1),upper=Inf,lower=0)
#' Geniefit(t1,Model="expo",start=c(100,.1),upper=Inf,lower=0,Rcpp=T)
#' ###a MCMC simulation included##
#' f=Geniefit(t1,Model="expo",start=c(100,.1),upper=Inf,lower=0,MCMC=T,sig=.1,run=10000)
#' acf(f$MCMC.simulation)
#' @author Fei Xiang (\email{xf3087@@gmail.com})
#'
#' @export
#'
#'
#'
#'



Geniefit=function(phy,Model,start,upper,lower,Rcpp=F,MCMC=F,sig=.1,run=10000){
  ###

  ###stops if number of parameters are not correct according to Model
  if (Model=="const") {
    if (length(start)!=1)  stop(paste("length of start must be 1 under a constant model", "\n", ""))
    }
  if (Model=="expo")  {
    if (length(start)!=2)  stop(paste("length of start must be 2 under an exponential model", "\n", ""))
    }
  if (Model=="expan") {
    if (length(start)!=3)  stop(paste("length of start must be 3 under an expansion model", "\n", ""))
    }
  if (Model=="log")   {
    if (length(start)!=3)  stop(paste("length of start must be 3 under a logistic model", "\n", ""))
    }
  if (Model=="step")  {
    if (length(start)!=3)  stop(paste("length of start must be 3 under a piecewise constant model", "\n", ""))
    }
  if (Model=="pexpan") {
    if (length(start)!=3)  stop(paste("length of start must be 3 under a piecewise expansion model", "\n", ""))
    }
  if (Model=="plog") {
    if (length(start)!=3)  stop(paste("length of start must be 3 under a piecewise logistic model", "\n", ""))
    }


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
      if (Model=="plog") {trajectory=ifelse(t<parr[3],parr[1],parr[1]*exp(-parr[2]*(t-parr[3])))}
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
      if (Model=="plog"){
        intg=(upperlim<parr[3])*(upperlim-lowerlim)/parr[1]+(lowerlim>parr[3])*(-exp(parr[2]*lowerlim-parr[2]*parr[3])+exp(parr[2]*upperlim-parr[2]*parr[3]))/parr[1]/parr[2]+(lowerlim<parr[3] && upperlim>parr[3])*((upperlim-parr[3])/parr[1]+(-exp(parr[2]*lowerlim-parr[2]*parr[3])+1)/parr[1]/parr[2])
      }
      return(intg)
    }
    logcoe=sort.type[-1]*(log(fnt(sort.times.pop[-1,1]))+log(pop.times[-ntotal]*(pop.times[-ntotal]-1)/2))
    logcoe[which(sort.type[-1]==0)]=0
    logint=-pop.times[-ntotal]*(pop.times[-ntotal]-1)/2*intfnt(sort.times.pop[-ntotal,1],sort.times.pop[-1,1])
    return(-sum(logcoe)-sum(logint))

  }

  if (Rcpp==F){
  fn2 <- function(x){
    fnpar(exp(x))
  }
  require(minqa)
  fit2 <- bobyqa(log(start),fn2,lower=log(lower),upper=log(upper))
  if (MCMC==T){
    mcmcsim=MCMCupdates(phy=phy,Model=Model,start=start,lower=lower,upper=upper,sig=sig,run=run)
    return(list(parr=exp(fit2$par),loglikelihood=-fit2$fval,AIC=2*length(start)+2*fit2$fval,MCMC.simulation=mcmcsim))
  } else{
    return(list(parr=exp(fit2$par),loglikelihood=-fit2$fval,AIC=2*length(start)+2*fit2$fval))
  }

  }


  if (Rcpp==T){
    require(dfoptim)
  x=att(phy)
  if (Model=="const")  fit.c=bobyqa(start,negllc_const,lower=lower,upper=upper,t=x$t,A=x$A)
  if (Model=="expo")   fit.c=nmkb(start,negllc_expo,t=x$t,A=x$A,upper = upper,lower = lower)
  if (Model=="pexpan") fit.c=nmkb(start,negllc_pexpan,t=x$t,A=x$A,upper = upper,lower = lower)
  if (Model=="log")    fit.c=nmkb(start,negllc_log,t=x$t,A=x$A,upper = upper,lower = lower)
  if (Model=="step")   fit.c=nmkb(start,negllc_step,t=x$t,A=x$A,upper = upper,lower = lower)
  if (Model=="expan")  fit.c=nmkb(start,negllc_expan,t=x$t,A=x$A,upper = upper,lower = lower)
  if (Model=="plog")   fit.c=nmkb(start,negllc_plog,t=x$t,A=x$A,upper = upper,lower = lower)

  if (MCMC==T){
    mcmcsim=MCMCupdates(phy=phy,Model=Model,start=start,lower=lower,upper=upper,sig=sig,run=run)
    return(list(parr=fit.c$par,loglikelihood=ifelse(Model=="const",-fit.c$fval,-fit.c$value),AIC=2*length(start)-2*ifelse(Model=="const",-fit.c$fval,-fit.c$value),MCMC.simulation=mcmcsim))
  } else{
    return(list(parr=fit.c$par,loglikelihood=ifelse(Model=="const",-fit.c$fval,-fit.c$value),AIC=2*length(start)-2*ifelse(Model=="const",-fit.c$fval,-fit.c$value)))
  }



  }



}



