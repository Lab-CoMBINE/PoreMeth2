#' Function to identify methylation state from Beta and Entropy Values
#'
#' description description description description description description
#'
#' @param TableIn 
#' @param PTStart numeric, 
#' @param ratiosd numeric, 
#' @param NormDist numeric, 
#' @param FW numeric, 
#'
#' @return what what what what what what what what what
#' 
#' @export

PoreMethState <- function(TableIn,PTStart=0.1,ratiosd=0.1,NormDist=1000000,FW=2)
{
  
  
  ChrIn<-TableIn$V1
  PosIn<-TableIn$V2
  EntropyValues<-TableIn$V3
  BetaValues<-TableIn$V5
  
  #########################################
  ###### Calculate Global Parameters  #####
  #########################################
  
  
  muk<-c(0,0.3,0.4,0.5,0.6,0.7,1)
  
  sdtot<-sd(BetaValues)
  sepsilon<-rep(sdtot*ratiosd,length(muk))
  sepsilon[2:6]<-sdtot*(1-ratiosd)
  
  MatrixMeth2Save<-c()
  
  ChrVec<-unique(ChrIn)
  
  ##############################################
  #### Starting Analysis across Chromosomes ####
  ##############################################
  
  for (jj in 1:length(ChrVec))
  {
    ChrSel<-ChrVec[jj]
    indC<-which(ChrIn==ChrSel)
    PosInC<-PosIn[indC]
    EntropyValuesC<-EntropyValues[indC]
    BetaValuesC<-BetaValues[indC]
    
    
    ############################
    ####### HMM  Analysis ######
    ############################
    
    
    Segmented<-PoreMeth2StateFinder(BetaValuesC,PosInC,muk,sepsilon,NormDist,PTStart)
    
    
    
    
    SegmentedS<-SortState(Segmented)
    callind<-SegmentedS[[3]]
    ResultSeg<-SegResults(rbind(BetaValuesC),callind)
    
    ##### Segmented Region Classification #####
    
    MatrixMethOut<-matrix(NA,ncol=8,nrow=length(callind)-1)
    for (pp in 1:(length(callind)-1)) 
    {
      PosDiffStart<-PosIn[callind[pp]+1]
      PosDiffEnd<-PosIn[callind[pp+1]]
      NumCPG<-callind[pp+1]-(callind[pp]+1)
      EntropyMean<-median(EntropyValuesC[(callind[pp]+1):callind[pp+1]])
      BetaMean<-ResultSeg[1,callind[pp]+1]
      
      if (BetaMean<0.5 & EntropyMean<=0.3)
      {
        SegmentedState<-"Unmethylated"
      }
      
      if (BetaMean>0.5 & EntropyMean<=0.3)
      {
        SegmentedState<-"Methylated"
      }
      
      if (EntropyMean>0.3)
      {
        SegmentedState<-"PMD"
      }
      MatrixMethOut[pp,]<-c(ChrVec[jj],PosDiffStart,PosDiffEnd,PosDiffEnd-PosDiffStart,ResultSeg[1,callind[pp]+1],EntropyMean,NumCPG,SegmentedState)
      
    }
    
    MatrixMethOutF<-MatrixMethOut[which(as.numeric(MatrixMethOut[,7])>=FW),]
    
    MatrixMeth2Save<-rbind(MatrixMeth2Save,MatrixMethOutF)
  }
  
  colnames(MatrixMeth2Save)<-c("chr","start","end","Size","Beta",
                               "Entropy","NumCpG","State")
  
  data.frame(MatrixMeth2Save)
}
