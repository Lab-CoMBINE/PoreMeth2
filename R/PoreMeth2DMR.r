#' Function for DMR, Beta and Entropy Detection
#'
#' Performs the double segmentation of Beta and S to identify DMRs and calculate Delta\Beta and Delta S between test and control sample. 
#'
#' @param TableTest Output from parse_nanopolish_entropy.pl, data.frame.
#' @param TableControl Output from parse_nanopolish_entropy.pl, data.frame.
#' @param omega Optional parameter that modulates the relative weight between the experimental and the biological variance. When omega is close to 1, the biological variance is much larger than the experimental one, while for values of omega close to 0 the experimental noise gives the leading contribution to the total variance. We suggest to use `omega` in the range 0.1-0.5, numeric.
#' @param eta Optional parameter that represents the baseline probability the mean process (m_i) changes its value for the HSLM algorithm. Suggested values are inside 10^{-7}-10^{−3} range, numeric.
#' @param FW The minimum number of datapoints for a DMR to be called (DMRs made of a number of CpGs smaller than `FW` are discarded), numeric.
#'
#' @return data.frame with identified DMR with corresponding Beta and Entropy values for each detected range.
#'
#' @export

PoreMeth2DMR <- function(TableTest,TableControl,omega=0.1,eta=1e-5,FW=3)
{
  
  PathSrc <- paste0(.libPaths()[1], "/PoreMeth2/libs/")
  dyn.load(paste0(PathSrc, "PoreMeth2.so"))

  ##### Extracting Info from Test and Control Tables ##### 
  ChrTest<-unlist(TableTest[,1])
  PosTest<-unlist(TableTest[,2])
  EntropyTest<-unlist(TableTest[,3])
  BetaTest<-unlist(TableTest[,5])
  
  ChrControl<-unlist(TableControl[,1])
  PosControl<-unlist(TableControl[,2])
  EntropyControl<-unlist(TableControl[,3])
  BetaControl<-unlist(TableControl[,5])
  
  
  ### Starting Analysis Chromosome by Chromosomes ####
  
  ChrVec<-intersect(unique(ChrTest),unique(ChrControl))
  
  MatrixDiffMethOut<-c()
  
  for (jj in 1:length(ChrVec))
  {
    ChrSel<-ChrVec[jj]
    indTestC<-which(ChrTest==ChrSel)
    PosTestC<-PosTest[indTestC]
    EntropyTestC<-EntropyTest[indTestC]
    BetaTestC<-BetaTest[indTestC]
    
    
    
    indControlC<-which(ChrControl==ChrSel)
    PosControlC<-PosControl[indControlC]
    EntropyControlC<-EntropyControl[indControlC]
    BetaControlC<-BetaControl[indControlC]
    
    
    ####### Finding common CpG between Test and Control Samples #####
    indMatch<-match(PosTestC,PosControlC)
    PosTestCMatch<-PosTestC[which(!is.na(indMatch))]
    PosControlCMatch<-PosControlC[indMatch[which(!is.na(indMatch))]]
    
    EntropyTestCMatch<-EntropyTestC[which(!is.na(indMatch))]
    EntropyControlCMatch<-EntropyControlC[indMatch[which(!is.na(indMatch))]]
    
    BetaTestCMatch<-BetaTestC[which(!is.na(indMatch))]
    BetaControlCMatch<-BetaControlC[indMatch[which(!is.na(indMatch))]]
    
    
    
    EntropyDiff<-EntropyTestCMatch-EntropyControlCMatch
    BetaDiff<-BetaTestCMatch-BetaControlCMatch
    
    #########################################
    ####### Bidimensional Segmentation ######
    #########################################
    
    MDSData<-rbind(BetaDiff,EntropyDiff)
    
        
    ### Calculating parameters of the HSLM ###
    ParamList <- ParamEstSeq(rbind(MDSData),omega)
    mi <- ParamList$mi
    smu <- ParamList$smu
    sepsilon <- ParamList$sepsilon
    muk <- MukEst(rbind(MDSData))
    
    ##### Performing bidimensional segmentation #####
    callind <- BidimensionalSLMSegIn(rbind(MDSData),muk,mi,smu,sepsilon,omega,eta)
    
    callind <- FilterSeg(callind,FW)
    DataSeg1 <- SegResults(rbind(MDSData),callind)
    
    
    MatrixDiffMeth<-matrix(NA,nrow=(length(callind)-1),ncol=11)
    for (pp in 1:(length(callind)-1)) 
    {
      PosDiffStart<-PosTestCMatch[callind[pp]+1]
      PosDiffEnd<-PosTestCMatch[callind[pp+1]]
      NumCPG<-callind[pp+1]-callind[pp]+1
      MatrixDiffMeth[pp,]<-c(ChrVec[jj],
                                PosDiffStart,
                                PosDiffEnd,
                                DataSeg1[1,callind[pp]+1],
                                DataSeg1[2,callind[pp]+1],
                                median(BetaTestCMatch[(callind[pp]+1):(callind[pp+1])]),
                                median(BetaControlCMatch[(callind[pp]+1):(callind[pp+1])]),
                                median(EntropyTestCMatch[(callind[pp]+1):(callind[pp+1])]),
                                median(EntropyControlCMatch[(callind[pp]+1):(callind[pp+1])]),
                                NumCPG,
                                wilcox.test(BetaTestCMatch[(callind[pp]+1):(callind[pp+1])],BetaControlCMatch[(callind[pp]+1):(callind[pp+1])])$p.value)
      
    }
    
    
    MatrixDiffMethOut<-rbind(MatrixDiffMethOut,MatrixDiffMeth)
  }
  
  colnames(MatrixDiffMethOut)<-c("chr","start","end","DeltaBeta","DeltaEntropy","BetaTest",
                                    "BetaControl","EntropyTest","EntropyControl","NumCpG","p")
  
  data.frame(MatrixDiffMethOut)
  
}
