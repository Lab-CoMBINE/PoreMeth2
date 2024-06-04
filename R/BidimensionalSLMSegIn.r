#' Truncated Gaussian Bidimensional SLM Segmentation
#' @NoRd


BidimensionalSLMSegIn <- function(DataMatrix,muk,mi,smu,sepsilon,omega,eta)
{
  ####  Calcolo il vettore delle covariate dipendenti dalla distanza ###

  PathSrc <- paste0(.libPaths()[1], "/PoreMeth2/libs/")
  dyn.load(paste0(PathSrc, "PoreMeth2.so"))

  T=ncol(DataMatrix)
  K0<-ncol(muk)
  etav<-log(rep(1,K0)*(1/K0))
  NExp<-nrow(DataMatrix)
  
  
  TruncCoef<-matrix(0,nrow=NExp,ncol=K0)
  for (jj in 1:NExp)
  {
    for (kk in 1:K0)
    {
      TruncCoef[jj,kk]<-sepsilon[jj]*(pnorm(1, mean = muk[jj,kk], sd = sepsilon[jj])-pnorm(-1, mean = muk[jj,kk], sd = sepsilon[jj]))
    }
  }
  
  
  ####  Transition and Emission Calculation ####
  P<-matrix(data=0,nrow=K0,ncol=K0)
  G<-matrix(data=0,nrow=K0,ncol=K0)
  emission<-matrix(data=0,nrow=K0,ncol=T)
  out<-.Fortran("transemisjslm",as.vector(muk),as.vector(mi),as.double(eta),as.matrix(DataMatrix),as.integer(K0),as.integer(NExp),as.vector(smu),as.vector(sepsilon),as.integer(T),as.matrix(G),as.matrix(P),as.matrix(emission),as.matrix(TruncCoef), PACKAGE = "PoreMeth2")
  
  P<-out[[11]]
  emission<-out[[12]]
  
  
  
  ####### Viterbi algorithm #########
  psi<-matrix(data=0,nrow=K0,ncol=T)
  path<-c(as.integer(rep(0,T)))
  out2 <- .Fortran("bioviterbi", as.vector(etav), as.matrix(P), as.matrix(emission), as.integer(T), as.integer(K0), as.vector(path), as.matrix(psi), PACKAGE = "PoreMeth2")
  s<-out2[[6]]
  
  
  sortResult <- SortState(s)
  TotalPredBreak<-sortResult[[3]]
  TotalPredBreak
}

