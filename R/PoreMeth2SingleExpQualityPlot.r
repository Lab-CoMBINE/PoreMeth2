#' Functions for printing plots for the evaluation of input data quality
#' 
#' This function prints plots for the evaluation of input data quality by displaying stats about `beta_cov` and `entropy_cov` (see Data Preparation) and beta and S density functions.
#'
#' @param TableIn Output from `parse_nanopolish_entropy.pl`
#'
#' @return Plot(s)
#'
#' @export


PoreMeth2SingleExpQualityPlot <- function(TableIn)
{

  
  Chr<-TableIn$V1
  Pos<-TableIn$V2
  BetaValues<-TableIn$V5
  EntropyValues<-TableIn$V3
  ReadsBeta<-TableIn$V6
  ReadsEntropy<-TableIn$V4
  
  
  ####################################
  ###### Calculating Statistics ######
  ####################################
  
  MeanCoverage<-median(ReadsBeta)
  MaxCoverage<-MeanCoverage*3
  
  breaks2Count<-c(1:MaxCoverage)
  ReadsBetaF<-ReadsBeta[which(ReadsBeta<MaxCoverage)]
  ReadsEntropyF<-ReadsEntropy[which(ReadsEntropy<MaxCoverage)]
  
  CountReadsBeta<-hist(ReadsBetaF,breaks = breaks2Count,plot=FALSE)$counts
  CountReadsEntropy<-hist(ReadsEntropyF,breaks = breaks2Count,plot=FALSE)$counts
  
  CumulativeCoverageBeta<-rev(cumsum(rev(CountReadsBeta)))
  CumulativeCoverageEntropy<-rev(cumsum(rev(CountReadsEntropy)))
  
  MaxCpG<-22000000
  
  par(mfrow=c(2,2))
  #########################################
  ####### Plotting Coverage vs Beta #######
  #########################################
  
  ProportionEpigenomeBeta10<-(CumulativeCoverageBeta/MaxCpG)[9]
  plot(breaks2Count[-1],CumulativeCoverageBeta/MaxCpG,pch=19,cex=1,xlab="Depth",ylab="Fraction of Covered Genome",ylim=c(0,1),xlim=c(0,MaxCoverage),main=expression(beta))
  lines(breaks2Count[-1],CumulativeCoverageBeta/MaxCpG,lwd=2,col="black")
  abline(v=10,lty=2,col="red")
  abline(h=ProportionEpigenomeBeta10,lty=2,col="blue")
  text(18,round(ProportionEpigenomeBeta10,digit=3)+0.05,paste(round(ProportionEpigenomeBeta10,digit=3)*100,"%"),col="blue")
  
  ############################################
  ####### Plotting Coverage vs Entropy #######
  ############################################
  
  ProportionEpigenomeEntropy10<-(CumulativeCoverageEntropy/MaxCpG)[9]
  plot(breaks2Count[-1],CumulativeCoverageEntropy/MaxCpG,pch=19,cex=1,xlab="Depth",ylab="Fraction of Covered Genome",ylim=c(0,1),xlim=c(0,MaxCoverage),main="ME")
  lines(breaks2Count[-1],CumulativeCoverageEntropy/MaxCpG,lwd=2,col="black")
  abline(v=10,lty=2,col="red")
  abline(h=ProportionEpigenomeEntropy10,lty=2,col="blue")
  text(18,round(ProportionEpigenomeEntropy10,digit=3)+0.05,paste(round(ProportionEpigenomeEntropy10,digit=3)*100,"%"),col="blue")
  
  
  ####################################################
  ######## Plotting  Beta and Entropy Density ########
  ####################################################
  
  plot(density(BetaValues,na.rm=TRUE),main="",xlab=expression(beta),lwd=2)
  
  plot(density(EntropyValues,na.rm=TRUE),main="",xlab="ME",lwd=2)
  
  
  
}
