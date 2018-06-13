#' A function to fit coalescent models to a given phylogenetic tree.
#' @description This function fits a user chosen coalescent model to fit a given phylogenetic tree by maximum likelihood estimation or Markov Chain Monte Carlo (MCMC) inference.
#'
#' Coalescent theory describes the mathematical properties
#' of intra-population phylogenies (Kingman (1982)). Griffiths and Tavare (1994) suggests a particular coalescent model,
#' which describes how the shape of the genealogy depends on the demograhic history of the sampled population. It provides
#' a probability distribution for the waiting times between coalescent events in a genealogy. This distribution depends on
#' the  demographic function N_e(t), which is the effective population size at time t before present. Pybus and Rambaut (2002)
#' represents the demographic function  using simple mathematical functions that characterize biologically plausible population
#' dynamic histories, such as constant, exponential, logistic, expension, piecewise constant, piecewise expansion,
#' and piecewise logistic growth.
#'
#' @param phy A phylogenetic tree, which is an object of class "phylo".
#' @param Model A model choice from const (constant population size), expo (exponetial growth),expan (expansion growth), log (logistic growth), step (piecewise constant), pexpan (piecewise expansion growth) and plog (piecewise logistic growth).
#' @param start A numeric vector of starting estimates of the parameters of the coalescent model.
#' @param lower A numeric vector of lower bounds on the parameters. If the length is 1 the single lower bound is applied to all parameters.
#' @param upper A numeric vector of upper bounds on the parameters. If the length is 1 the single upper bound is applied to all parameters.
#' @param Rcpp Calculation is based on C++ code when it is True and on R code when it is False.
#' @param MCMC MCMC simulation is run when it is true, and not run when it is False. The default prior is uniform given the lower and upper.
#' @param sig MCMC simulation step size.
#' @param run Number of MCMC simulation.
#' @return Parameters estimation of a given model, loglikelihood and AIC, when MCMC=F and additional MCMC simulations for parameters and loglikelihood will be presented when MCMC=T.
#' @examples
#' library(ape)
#' t1=rcoal(20)
#' geniefit(t1,Model="expo",start=c(100,.1),upper=Inf,lower=0)
#' geniefit(t1,Model="expo",start=c(100,.1),upper=Inf,lower=0,Rcpp=T)
#' ###a MCMC simulation included##
#' f=geniefit(t1,Model="expo",start=c(100,.1),upper=Inf,lower=0,MCMC=T,sig=.1,run=10000)
#' acf(f$MCMC.simulation)
#' @author Fei Xiang (\email{xf3087@@gmail.com})
#' @references
#' Kingman, J. F. C. (1982). On the Genealogy of Large Populations. Journal of Applied Probability 19, 27-43.
#'
#' Griffiths, R., and Tavare, S. (1994). Sampling Theory for Neutral Alleles in a Varying Environment. Philosophical Transactions: Biological Sciences, 344(1310), 403-410.
#'
#' Pybus, O. G., and Rambaut, A. (2002). GENIE: Estimating Demographic History from Molecular Phylogenies. Bioinformatics 18, 1404-1405.
#' @export
#'
#'
#'
#'



geniefit=function(phy,Model,start,upper,lower,Rcpp=F,MCMC=F,sig=.1,run=10000){
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
    mcmcsim=MCMC.updates(phy=phy,Model=Model,start=start,lower=lower,upper=upper,sig=sig,run=run)
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
    mcmcsim=MCMC.updates(phy=phy,Model=Model,start=start,lower=lower,upper=upper,sig=sig,run=run)
    return(list(parr=fit.c$par,loglikelihood=ifelse(Model=="const",-fit.c$fval,-fit.c$value),AIC=2*length(start)-2*ifelse(Model=="const",-fit.c$fval,-fit.c$value),MCMC.simulation=mcmcsim))
  } else{
    return(list(parr=fit.c$par,loglikelihood=ifelse(Model=="const",-fit.c$fval,-fit.c$value),AIC=2*length(start)-2*ifelse(Model=="const",-fit.c$fval,-fit.c$value)))
  }



  }



}



