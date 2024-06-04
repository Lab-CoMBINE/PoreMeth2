#' Estimating Starting Parameters
#' @NoRd


ParamEstSeq <- function(DataMatrix,omega)
{
  T=ncol(DataMatrix)
  NExp<-nrow(DataMatrix)
  sigmax<-c()
  mi<-c()
  
  for (i in 1:NExp)
  {
    mi[i]<-0
    sigmax[i]<-mad(DataMatrix[i,])^2
  }
  
  smu<-sqrt(omega*sigmax)
  sepsilon<-sqrt((1-omega)*sigmax)
  Results<-list()
  Results$mi<-mi
  Results$smu<-smu
  Results$sepsilon<-sepsilon
  
  Results
}

