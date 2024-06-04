#' MUK Estimator
#' @NoRd


MukEst <- function(DataMatrix)
{
  NExp<-dim(DataMatrix)[1]
  if (NExp==1)
  {
    muk<-rbind(seq(-1,1,by=0.1))
  }
  if (NExp==2)
  {
    
    binsize=0.1
    binVec<-c(seq(-1,-0.1,by=binsize),0,seq(0.1,1,by=binsize))
    muk<-t(expand.grid( binVec, binVec))
    
  }
  muk
}
