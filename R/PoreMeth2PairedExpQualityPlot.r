#' Display `beta_cov` and `entropy_cov` stats for paired samples.
#'
#' With this function it is possible to display stats about `beta_cov` and `entropy_cov` for common CpG dineucleotides in the Test/Control samples pair.
#'
#' @param TableTestIn Output from `parse_nanopolish_entropy.pl`, data.frame.
#' @param TableControlIn Output from `parse_nanopolish_entropy.pl`, data.frame.
#'
#' @return Plot(s)
#'
#' @export



PoreMeth2PairedExpQualityPlot <- function(TableTestIn,TableControlIn)
{
  ##### Info Test #####
  ChrTest<-TableTestIn$V1
  PosTest<-TableTestIn$V2
  BetaValuesTest<-TableTestIn$V5
  EntropyValuesTest<-TableTestIn$V3
  ReadsBetaTest<-TableTestIn$V6
  ReadsEntropyTest<-TableTestIn$V4
  
  ##### Info Control #####
  ChrControl<-TableControlIn$V1
  PosControl<-TableControlIn$V2
  BetaValuesControl<-TableControlIn$V5
  EntropyValuesControl<-TableControlIn$V3
  ReadsBetaControl<-TableControlIn$V6
  ReadsEntropyControl<-TableControlIn$V4
  
  
  
  
  ####################################
  ###### Calculating Statistics ######
  ####################################
  
  MeanCoverageTest<-median(ReadsBetaTest)
  MeanCoverageControl<-median(ReadsBetaControl)
  
  MaxCoverage<-max(MeanCoverageTest*3,MeanCoverageControl*3)
  
  
  breaks2Count<-c(1:MaxCoverage)
  
  
  
  
  
  
  ################################################################
  ###### Searching Common Position, Beta and Entropy Values ######
  ################################################################
  
  
  indTestF<-which(ReadsBetaTest<MaxCoverage)
  indControlF<-which(ReadsBetaControl<MaxCoverage)
  
  ChrControlF<-ChrControl[indControlF]
  PosControlF<-PosControl[indControlF]
  ReadsBetaControlF<-ReadsBetaControl[indControlF]
  ReadsEntropyControlF<-ReadsEntropyControl[indControlF]
  
  ChrTestF<-ChrTest[indTestF]
  PosTestF<-PosTest[indTestF]
  ReadsBetaTestF<-ReadsBetaTest[indTestF]
  ReadsEntropyTestF<-ReadsEntropyTest[indTestF]
  
  
  CountReadsBetaMinCommon<-rep(0,length(breaks2Count)-1)
  CountReadsEntropyMinCommon<-rep(0,length(breaks2Count)-1)
  
  
  ChrCommonU<-intersect(unique(ChrTestF),unique(ChrControlF))
  for (xx in 1:length(ChrCommonU))
  {
    
    indTestFC<-which(ChrTestF==ChrCommonU[xx])
    indControlFC<-which(ChrControlF==ChrCommonU[xx])
    
    PosTestFC<-PosTestF[indTestFC]
    ReadsBetaTestFC<-ReadsBetaTestF[indTestFC]
    ReadsEntropyTestFC<-ReadsEntropyTestF[indTestFC]
    
    PosControlFC<-PosControlF[indControlFC]
    ReadsBetaControlFC<-ReadsBetaControlF[indControlFC]
    ReadsEntropyControlFC<-ReadsEntropyControlF[indControlFC]
    
    indMatch<-match(PosTestFC,PosControlFC)
    
    ReadsBetaControlFCMatch<-ReadsBetaControlFC[indMatch[which(!is.na(indMatch))]]
    ReadsEntropyControlFCMatch<-ReadsEntropyControlFC[indMatch[which(!is.na(indMatch))]]
    
    ReadsBetaTestFCMatch<-ReadsBetaTestFC[which(!is.na(indMatch))]
    ReadsEntropyTestFCMatch<-ReadsEntropyTestFC[which(!is.na(indMatch))]
    
    
    
    
    
    
    ReadsBetaMinFCMatch<-apply(cbind(ReadsBetaTestFCMatch,ReadsBetaControlFCMatch), 1, FUN = min)
    ReadsEntropyMinFCMatch<-apply(cbind(ReadsEntropyTestFCMatch,ReadsEntropyControlFCMatch), 1, FUN = min)
    
    
    CountReadsBetaMinCommon<-CountReadsBetaMinCommon+hist(ReadsBetaMinFCMatch,breaks = breaks2Count,plot=FALSE)$counts
    
    CountReadsEntropyMinCommon<-CountReadsEntropyMinCommon+hist(ReadsEntropyMinFCMatch,breaks = breaks2Count,plot=FALSE)$counts
  }
  
  
  
  
  
  
  CumulativeCoverageBetaCommon<-rev(cumsum(rev(CountReadsBetaMinCommon)))
  CumulativeCoverageEntropyCommon<-rev(cumsum(rev(CountReadsEntropyMinCommon)))
  
  
  MaxCpG<-22000000
  
  par(mfrow=c(1,2))
  
  
  
  ProportionEpigenomeBetaCommon10<-(CumulativeCoverageBetaCommon/MaxCpG)[9]
  plot(breaks2Count[-1],CumulativeCoverageBetaCommon/MaxCpG,pch=19,cex=1,xlab="Depth",ylab="Fraction of Covered CpGs",ylim=c(0,1),xlim=c(0,MaxCoverage),main=expression(beta))
  lines(breaks2Count[-1],CumulativeCoverageBetaCommon/MaxCpG,lwd=2,col="black")
  abline(v=10,lty=2,col="red")
  abline(h=ProportionEpigenomeBetaCommon10,lty=2,col="blue")
  text(18,round(ProportionEpigenomeBetaCommon10,digit=3)+0.05,paste(round(ProportionEpigenomeBetaCommon10,digit=3)*100,"%"),col="blue")
  
  ProportionEpigenomeBetaCommon10<-(CumulativeCoverageEntropyCommon/MaxCpG)[9]
  plot(breaks2Count[-1],CumulativeCoverageEntropyCommon/MaxCpG,pch=19,cex=1,xlab="Depth",ylab="Fraction of Covered CpGs",ylim=c(0,1),xlim=c(0,MaxCoverage),main=expression(Entropy))
  lines(breaks2Count[-1],CumulativeCoverageEntropyCommon/MaxCpG,lwd=2,col="black")
  abline(v=10,lty=2,col="red")
  abline(h=ProportionEpigenomeBetaCommon10,lty=2,col="blue")
  text(18,round(ProportionEpigenomeBetaCommon10,digit=3)+0.05,paste(round(ProportionEpigenomeBetaCommon10,digit=3)*100,"%"),col="blue")
  
}

