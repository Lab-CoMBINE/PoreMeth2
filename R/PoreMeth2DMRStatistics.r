#' Summarize Results 
#'
#' Returns a summary of hyper/hypo-methylated regions with different values of Delta S across different genic features, CpG Islands and enhancers.
#'
#' @param TableInDMR Is the output from `PoreMeth2DMR`, data.frame.
#' @param Assembly Optional parameter that specifies the reference version to use for statistics (`"hg19"` or `"hg38"`), character.
#' @param BetaThr Delta beta threshold applied for DMRs' classification, numeric.
#' @param EntropyThr Delta S threshold applied for DMRs' classification, numeric.
#' @param PValueThr The p.value threshold to consider a DMR, numeric.
#' @param AnalysisClass define the summary to return ("All", "Beta", "Entropy"). Default = "All", character.
#' 
#' @return data.frame with Summarized Results.
#'
#' @export

PoreMeth2DMRStatistics <- function(TableInDMR,Assembly="hg19",
                                  BetaThr=0.2,
                                  EntropyThr=0.1,
                                  PValueThr=0.05,
                                  AnalysisClass="All")
{
  
  PathSrc <- paste0(.libPaths()[1], "/PoreMeth2/libs/")
  dyn.load(paste0(PathSrc, "PoreMeth2.so"))


  PathDBIn <- paste0(.libPaths()[1], "/PoreMeth2/data/")
  
  indF<-which(abs(as.numeric(TableInDMR$DeltaBeta))>BetaThr & as.numeric(TableInDMR$p)<PValueThr)
  TableInDMR<-TableInDMR[indF,]
  
  
  
  
  
  ChrDMRString<-TableInDMR$chr[1]
  


  if (nchar(ChrDMRString)>2)
  {
    ChrDMR<-as.character(TableInDMR$chr)
  }
  if (nchar(ChrDMRString)<2)
  {
    
    ChrDMR<-paste("chr",as.character(TableInDMR$chr),sep="")
  }
  
  StartDMR<-as.numeric(TableInDMR$start)
  EndDMR<-as.numeric(TableInDMR$end)
  DeltaBetaDMR<-as.numeric(TableInDMR$DeltaBeta)
  DeltaEntropyDMR<-as.numeric(TableInDMR$DeltaEntropy)
  
  #################################################
  #################################################
  ###### Loading Genomic Features Databases #######
  #################################################
  #################################################
  
  ###########################
  ##### Flat Intergenic #####
  ###########################
  
  
  FileIntergenicIn<-file.path(PathDBIn,paste("InterGenic_",Assembly,".rds",sep=""))
  TableIntergenicIn<-readRDS(FileIntergenicIn)
  
  ChrIntergenic<-TableIntergenicIn$Chr
  StartIntergenic<-as.numeric(TableIntergenicIn$Start)
  EndIntergenic<-as.numeric(TableIntergenicIn$End)
  
  #######################
  ##### Flat GenCode #####
  #######################
  
  
  FileGenCodeIn<-file.path(PathDBIn,paste("GencodeTable_",Assembly,".rds",sep=""))
  TableGenCodeIn<-readRDS(FileGenCodeIn)
  
  ChrGenCode<-TableGenCodeIn$Chr
  StartGenCode<-as.numeric(TableGenCodeIn$Start)
  EndGenCode<-as.numeric(TableGenCodeIn$End)
  FeatureGenCode<-TableGenCodeIn$`Genic Element`
  StrandGenCode<-TableGenCodeIn$Strand
  SymbolGenCode<-TableGenCodeIn$`Gene Symbol`
  TypeGenCode<-TableGenCodeIn$`Transcript Type`
  
  FeatureGenCodeSub<-substr(FeatureGenCode,start=1,stop=4)
  
  
  ############################
  ##### File CpG Islands #####
  ############################
  
  FileCGIIn<-file.path(PathDBIn,paste("CGIAnno_",Assembly,".rds",sep=""))
  TableCGIIn<-readRDS(FileCGIIn)
  
  ChrCGI<-as.character(TableCGIIn$Chr)
  StartCGI<-as.numeric(TableCGIIn$Start)
  EndCGI<-as.numeric(TableCGIIn$End)
  NameCGI<-as.character(TableCGIIn$FeatureName)
  
  
  #################################
  ######## File Enhancers  ########
  #################################
  
  FileEnhancerIn<-file.path(PathDBIn,paste("EnhancerAnno_",Assembly,".rds",sep=""))
  TableEnhancerIn<-readRDS(FileEnhancerIn)
  
  ChrEnhancer<-as.character(TableEnhancerIn$Chr)
  StartEnhancer<-as.numeric(TableEnhancerIn$Start)
  EndEnhancer<-as.numeric(TableEnhancerIn$End)
  NameEnhancer<-as.character(TableEnhancerIn$FeatureName)
  
  
  
  
  
  ###################################
  ####### Generating OutPut #########
  ###################################
  
  #################################
  ###### Analysis Chr by Chr ######
  #################################
  FeatureNameOut<-c("Intergenic","Promoter","First Exon","Intron","Exon","Last Exon","CGI","Enhancer")
  BaseNameFeature<-c("Prom","Firs","Intr","Exon","Last")
  
  
  if (AnalysisClass=="Entropy")
  {
    CoordDMRCHighSumGencode<-rep(0,length(BaseNameFeature))
    CoordDMRCMidSumGencode<-rep(0,length(BaseNameFeature))
    CoordDMRCLowSumGencode<-rep(0,length(BaseNameFeature))
    
    CoordDMRCHighSumTotal<-0
    CoordDMRCMidSumTotal<-0
    CoordDMRCLowSumTotal<-0
    
    CoordDMRCHighSumCGI<-0
    CoordDMRCMidSumCGI<-0
    CoordDMRCLowSumCGI<-0
    
    CoordDMRCHighSumIntergenic<-0
    CoordDMRCMidSumIntergenic<-0
    CoordDMRCLowSumIntergenic<-0
    
    CoordDMRCHighSumEnhancer<-0
    CoordDMRCMidSumEnhancer<-0
    CoordDMRCLowSumEnhancer<-0
  }
  
  
  
  if (AnalysisClass=="All")
  {
    CoordDMRCHyperHighSumGencode<-rep(0,length(BaseNameFeature))
    CoordDMRCHypoHighSumGencode<-rep(0,length(BaseNameFeature))
    CoordDMRCHyperMidSumGencode<-rep(0,length(BaseNameFeature))
    CoordDMRCHypoMidSumGencode<-rep(0,length(BaseNameFeature))
    CoordDMRCHyperLowSumGencode<-rep(0,length(BaseNameFeature))
    CoordDMRCHypoLowSumGencode<-rep(0,length(BaseNameFeature))
    CoordDMRCHyperHighSumTotal<- 0
    CoordDMRCHypoHighSumTotal<- 0
    CoordDMRCHyperMidSumTotal<- 0
    CoordDMRCHypoMidSumTotal<- 0
    CoordDMRCHyperLowSumTotal<- 0
    CoordDMRCHypoLowSumTotal<- 0
    
    CoordDMRCHyperHighSumCGI<- 0
    CoordDMRCHypoHighSumCGI<- 0
    CoordDMRCHyperMidSumCGI<- 0
    CoordDMRCHypoMidSumCGI<- 0
    CoordDMRCHyperLowSumCGI<- 0
    CoordDMRCHypoLowSumCGI<- 0
    
    CoordDMRCHyperHighSumIntergenic<- 0
    CoordDMRCHypoHighSumIntergenic<- 0
    CoordDMRCHyperMidSumIntergenic<- 0
    CoordDMRCHypoMidSumIntergenic<- 0
    CoordDMRCHyperLowSumIntergenic<- 0
    CoordDMRCHypoLowSumIntergenic<- 0
    
    CoordDMRCHyperHighSumEnhancer<- 0
    CoordDMRCHypoHighSumEnhancer<- 0
    CoordDMRCHyperMidSumEnhancer<- 0
    CoordDMRCHypoMidSumEnhancer<- 0
    CoordDMRCHyperLowSumEnhancer<- 0
    CoordDMRCHypoLowSumEnhancer<- 0
  }
  
  if (AnalysisClass=="Beta")
  {
    
    CoordDMRCHyperSumGencode<-rep(0,length(BaseNameFeature))
    CoordDMRCHypoSumGencode<-rep(0,length(BaseNameFeature))
    CoordDMRCHyperSumTotal<-0
    CoordDMRCHypoSumTotal<-0
    CoordDMRCHyperSumCGI<-0
    CoordDMRCHypoSumCGI<-0
    
    CoordDMRCHyperSumIntergenic<-0
    CoordDMRCHypoSumIntergenic<-0
    
    CoordDMRCHyperSumEnhancer<-0
    CoordDMRCHypoSumEnhancer<-0
    
  }
  
  #dyn.load("/home/alberto/NanoporeMethylation/AnnotationTool/OverlapLibrary.so")
  
  
  ChrDMRInU<-unique(ChrDMR)
  for (i in 1:length(ChrDMRInU))
  {
    
    
    #### Selecting Chr Coordinates for DMRs ####
    indDMRC<-which(ChrDMR==ChrDMRInU[i])
    StartDMRC<-StartDMR[indDMRC]
    EndDMRC<-EndDMR[indDMRC]
    DeltaBetaDMRC<-DeltaBetaDMR[indDMRC]
    DeltaEntropyDMRC<-DeltaEntropyDMR[indDMRC]
    CoordDMRC<-cbind(as.numeric(StartDMRC),as.numeric(EndDMRC))
    
    ##################################################
    #### Selecting Chr Coordinates for Intergenic ####
    ##################################################
    indIntergenicC<-which(ChrIntergenic==ChrDMRInU[i])
    StartIntergenicC<-StartIntergenic[indIntergenicC]
    EndIntergenicC<-EndIntergenic[indIntergenicC]
    
    CoordIntergenicC<-cbind(as.numeric(StartIntergenicC),as.numeric(EndIntergenicC))
    
    ###############################################
    #### Selecting Chr Coordinates for GenCode ####
    ###############################################
    
    
    
    indGenCodeC<-which(ChrGenCode==ChrDMRInU[i])
    StartGenCodeC<-StartGenCode[indGenCodeC]
    EndGenCodeC<-EndGenCode[indGenCodeC]
    FeatureGenCodeC<-FeatureGenCode[indGenCodeC]
    StrandGenCodeC<-StrandGenCode[indGenCodeC]
    SymbolGenCodeC<-SymbolGenCode[indGenCodeC]
    TypeGenCodeC<-TypeGenCode[indGenCodeC]
    FeatureGenCodeSubC<-FeatureGenCodeSub[indGenCodeC]
    
    CoordGenCodeC<-cbind(as.numeric(StartGenCodeC),as.numeric(EndGenCodeC))
    
    ###########################################
    #### Selecting Chr Coordinates for CGI ####
    ###########################################
    
    indCGIC<-which(ChrCGI==ChrDMRInU[i])
    StartCGIC<-StartCGI[indCGIC]
    EndCGIC<-EndCGI[indCGIC]
    NameCGIC<-NameCGI[indCGIC]
    
    CoordCGIC<-cbind(as.numeric(StartCGIC),as.numeric(EndCGIC))
    
    ###########################################
    #### Selecting Chr Coordinates for Enhancer ####
    ###########################################
    
    indEnhancerC<-which(ChrEnhancer==ChrDMRInU[i])
    StartEnhancerC<-StartEnhancer[indEnhancerC]
    EndEnhancerC<-EndEnhancer[indEnhancerC]
    NameEnhancerC<-NameEnhancer[indEnhancerC]
    
    CoordEnhancerC<-cbind(as.numeric(StartEnhancerC),as.numeric(EndEnhancerC))
    
    if (AnalysisClass=="All")
    {
      CoordDMRCHyperHigh<-rbind(CoordDMRC[which(DeltaBetaDMRC>0 & DeltaEntropyDMRC>=EntropyThr),])
      CoordDMRCHypoHigh<-rbind(CoordDMRC[which(DeltaBetaDMRC<0 & DeltaEntropyDMRC>=EntropyThr),])
      CoordDMRCHyperMid<-rbind(CoordDMRC[which(DeltaBetaDMRC>0 & DeltaEntropyDMRC> -EntropyThr & DeltaEntropyDMRC< EntropyThr),])
      CoordDMRCHypoMid<-rbind(CoordDMRC[which(DeltaBetaDMRC<0 & DeltaEntropyDMRC> -EntropyThr & DeltaEntropyDMRC< EntropyThr),])
      CoordDMRCHyperLow<-rbind(CoordDMRC[which(DeltaBetaDMRC>0 & DeltaEntropyDMRC<= -EntropyThr),])
      CoordDMRCHypoLow<-rbind(CoordDMRC[which(DeltaBetaDMRC<0 & DeltaEntropyDMRC<= -EntropyThr),])
      
      
      if (length(CoordDMRCHyperHigh)>0)
      {
        CoordDMRCHyperHighSumTotal<- CoordDMRCHyperHighSumTotal + sum(CoordDMRCHyperHigh[,2]-CoordDMRCHyperHigh[,1])
      }
      if (length(CoordDMRCHypoHigh)>0)
      {
        CoordDMRCHypoHighSumTotal<- CoordDMRCHypoHighSumTotal + sum(CoordDMRCHypoHigh[,2]-CoordDMRCHypoHigh[,1])
      }
      if (length(CoordDMRCHyperMid)>0)
      {
        CoordDMRCHyperMidSumTotal<- CoordDMRCHyperMidSumTotal + sum(CoordDMRCHyperMid[,2]-CoordDMRCHyperMid[,1])
      }
      if (length(CoordDMRCHypoMid)>0)
      {
        CoordDMRCHypoMidSumTotal<- CoordDMRCHypoMidSumTotal + sum(CoordDMRCHypoMid[,2]-CoordDMRCHypoMid[,1])
      }
      if (length(CoordDMRCHyperLow)>0)
      {
        CoordDMRCHyperLowSumTotal<- CoordDMRCHyperLowSumTotal + sum(CoordDMRCHyperLow[,2]-CoordDMRCHyperLow[,1])
      }
      if (length(CoordDMRCHypoLow)>0)
      {
        CoordDMRCHypoLowSumTotal<- CoordDMRCHypoLowSumTotal + sum(CoordDMRCHypoLow[,2]-CoordDMRCHypoLow[,1])
      }
      
      ###########################################
      #### Searching CGI and DMRs Overlap #######
      ###########################################
      
      if (length(CoordDMRCHyperHigh)>0)
      {
      N<-nrow(CoordDMRCHyperHigh)
      M<-nrow(CoordCGIC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyperHigh), as.matrix(CoordCGIC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHyperHighSumCGI = CoordDMRCHyperHighSumCGI + out[[5]]
      }
      
      if (length(CoordDMRCHypoHigh)>0)
      {
        
      N<-nrow(CoordDMRCHypoHigh)
      M<-nrow(CoordCGIC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypoHigh), as.matrix(CoordCGIC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHypoHighSumCGI = CoordDMRCHypoHighSumCGI + out[[5]]
      }
      
      if (length(CoordDMRCHyperMid)>0)
      {
        
      N<-nrow(CoordDMRCHyperMid)
      M<-nrow(CoordCGIC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyperMid), as.matrix(CoordCGIC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHyperMidSumCGI = CoordDMRCHyperMidSumCGI + out[[5]]
      }
      
      if (length(CoordDMRCHypoMid)>0)
      {
        
      N<-nrow(CoordDMRCHypoMid)
      M<-nrow(CoordCGIC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypoMid), as.matrix(CoordCGIC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHypoMidSumCGI = CoordDMRCHypoMidSumCGI + out[[5]]
      }
      
      if (length(CoordDMRCHyperLow)>0)
      {
        
      N<-nrow(CoordDMRCHyperLow)
      M<-nrow(CoordCGIC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyperLow), as.matrix(CoordCGIC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHyperLowSumCGI = CoordDMRCHyperLowSumCGI + out[[5]]
      }
      
      if (length(CoordDMRCHypoLow)>0)
      {
        
      N<-nrow(CoordDMRCHypoLow)
      M<-nrow(CoordCGIC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypoLow), as.matrix(CoordCGIC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHypoLowSumCGI = CoordDMRCHypoLowSumCGI + out[[5]]
      }
      
      ###########################################
      #### Searching Intergenic and DMRs Overlap #######
      ###########################################
      
      if (length(CoordDMRCHyperHigh)>0)
      {
        N<-nrow(CoordDMRCHyperHigh)
        M<-nrow(CoordIntergenicC)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyperHigh), as.matrix(CoordIntergenicC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHyperHighSumIntergenic = CoordDMRCHyperHighSumIntergenic + out[[5]]
      }
      
      if (length(CoordDMRCHypoHigh)>0)
      {
        
        N<-nrow(CoordDMRCHypoHigh)
        M<-nrow(CoordIntergenicC)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypoHigh), as.matrix(CoordIntergenicC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHypoHighSumIntergenic = CoordDMRCHypoHighSumIntergenic + out[[5]]
      }
      
      if (length(CoordDMRCHyperMid)>0)
      {
        
        N<-nrow(CoordDMRCHyperMid)
        M<-nrow(CoordIntergenicC)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyperMid), as.matrix(CoordIntergenicC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHyperMidSumIntergenic = CoordDMRCHyperMidSumIntergenic + out[[5]]
      }
      
      if (length(CoordDMRCHypoMid)>0)
      {
        
        N<-nrow(CoordDMRCHypoMid)
        M<-nrow(CoordIntergenicC)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypoMid), as.matrix(CoordIntergenicC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHypoMidSumIntergenic = CoordDMRCHypoMidSumIntergenic + out[[5]]
      }
      
      if (length(CoordDMRCHyperLow)>0)
      {
        
        N<-nrow(CoordDMRCHyperLow)
        M<-nrow(CoordIntergenicC)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyperLow), as.matrix(CoordIntergenicC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHyperLowSumIntergenic = CoordDMRCHyperLowSumIntergenic + out[[5]]
      }
      
      if (length(CoordDMRCHypoLow)>0)
      {
        
        N<-nrow(CoordDMRCHypoLow)
        M<-nrow(CoordIntergenicC)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypoLow), as.matrix(CoordIntergenicC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHypoLowSumIntergenic = CoordDMRCHypoLowSumIntergenic + out[[5]]
      }
      ###########################################
      #### Searching Enhancer and DMRs Overlap #######
      ###########################################
      if (length(CoordDMRCHyperHigh)>0)
      {
        
      N<-nrow(CoordDMRCHyperHigh)
      M<-nrow(CoordEnhancerC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyperHigh), as.matrix(CoordEnhancerC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHyperHighSumEnhancer = CoordDMRCHyperHighSumEnhancer + out[[5]]
      }
      
      if (length(CoordDMRCHypoHigh)>0)
      {
        
      N<-nrow(CoordDMRCHypoHigh)
      M<-nrow(CoordEnhancerC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypoHigh), as.matrix(CoordEnhancerC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHypoHighSumEnhancer = CoordDMRCHypoHighSumEnhancer + out[[5]]
      }
      if (length(CoordDMRCHyperMid)>0)
      {
        
      
      N<-nrow(CoordDMRCHyperMid)
      M<-nrow(CoordEnhancerC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyperMid), as.matrix(CoordEnhancerC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHyperMidSumEnhancer = CoordDMRCHyperMidSumEnhancer + out[[5]]
      }
      
      if (length(CoordDMRCHypoMid)>0)
      {
        
      N<-nrow(CoordDMRCHypoMid)
      M<-nrow(CoordEnhancerC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypoMid), as.matrix(CoordEnhancerC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHypoMidSumEnhancer = CoordDMRCHypoMidSumEnhancer + out[[5]]
      }
      
      if (length(CoordDMRCHyperLow)>0)
      {
        
      N<-nrow(CoordDMRCHyperLow)
      M<-nrow(CoordEnhancerC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyperLow), as.matrix(CoordEnhancerC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHyperLowSumEnhancer = CoordDMRCHyperLowSumEnhancer + out[[5]]
      }
      
      if (length(CoordDMRCHypoLow)>0)
      {
        
      N<-nrow(CoordDMRCHypoLow)
      M<-nrow(CoordEnhancerC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypoLow), as.matrix(CoordEnhancerC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHypoLowSumEnhancer = CoordDMRCHypoLowSumEnhancer + out[[5]]
      }
      #########################################################
      ##### Calculating Total DMR in each GenCode Feature #####
      #########################################################
      
      
      
      for (kk in 1:length(BaseNameFeature))
      {
        ind2Feat<-which(FeatureGenCodeSubC==BaseNameFeature[kk])
        
        CoordGenCodeC2Feat<-CoordGenCodeC[ind2Feat,]
        
        #### Searching GenCode and DMRs Overlap #######
        if (length(CoordDMRCHyperHigh)>0)
        {
        N<-nrow(CoordDMRCHyperHigh)
        M<-nrow(CoordGenCodeC2Feat)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyperHigh), as.matrix(CoordGenCodeC2Feat), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHyperHighSumGencode[kk] = CoordDMRCHyperHighSumGencode[kk]+ out[[5]]
        }
        
        #### Searching GenCode and DMRs Overlap #######
        if (length(CoordDMRCHypoHigh)>0)
        {
          
        N<-nrow(CoordDMRCHypoHigh)
        M<-nrow(CoordGenCodeC2Feat)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypoHigh), as.matrix(CoordGenCodeC2Feat), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHypoHighSumGencode[kk] = CoordDMRCHypoHighSumGencode[kk] + out[[5]]
        }
        #### Searching GenCode and DMRs Overlap #######
        if (length(CoordDMRCHyperMid)>0)
        {
          
        N<-nrow(CoordDMRCHyperMid)
        M<-nrow(CoordGenCodeC2Feat)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyperMid), as.matrix(CoordGenCodeC2Feat), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHyperMidSumGencode[kk] = CoordDMRCHyperMidSumGencode[kk] + out[[5]]
        }
        
        #### Searching GenCode and DMRs Overlap #######
        if (length(CoordDMRCHypoMid)>0)
        {
          
        N<-nrow(CoordDMRCHypoMid)
        M<-nrow(CoordGenCodeC2Feat)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypoMid), as.matrix(CoordGenCodeC2Feat), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHypoMidSumGencode[kk] = CoordDMRCHypoMidSumGencode[kk] + out[[5]]
        }
        #### Searching GenCode and DMRs Overlap #######
        if (length(CoordDMRCHyperLow)>0)
        {
          
        N<-nrow(CoordDMRCHyperLow)
        M<-nrow(CoordGenCodeC2Feat)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyperLow), as.matrix(CoordGenCodeC2Feat), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHyperLowSumGencode[kk] = CoordDMRCHyperLowSumGencode[kk] + out[[5]]
        }
        
        #### Searching GenCode and DMRs Overlap #######
        if (length(CoordDMRCHypoLow)>0)
        {
          
        N<-nrow(CoordDMRCHypoLow)
        M<-nrow(CoordGenCodeC2Feat)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypoLow), as.matrix(CoordGenCodeC2Feat), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHypoLowSumGencode[kk] = CoordDMRCHypoLowSumGencode[kk] + out[[5]]
        }
      }
      
      
      
    }
    
    
    
    if (AnalysisClass=="Beta")
    {
      CoordDMRCHyper<-rbind(CoordDMRC[which(DeltaBetaDMRC>0),])
      CoordDMRCHypo<-rbind(CoordDMRC[which(DeltaBetaDMRC<0),])
      
      
      
      if (length(CoordDMRCHyper)>0)
      {
      CoordDMRCHyperSumTotal<- CoordDMRCHyperSumTotal + sum(CoordDMRCHyper[,2]-CoordDMRCHyper[,1])
      }
      if (length(CoordDMRCHypo)>0)
      {
      CoordDMRCHypoSumTotal<- CoordDMRCHypoSumTotal + sum(CoordDMRCHypo[,2]-CoordDMRCHypo[,1])
      }
      
      ###########################################
      #### Searching CGI and DMRs Overlap #######
      ###########################################
      if (length(CoordDMRCHyper)>0)
      {
        
      N<-nrow(CoordDMRCHyper)
      M<-nrow(CoordCGIC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyper), as.matrix(CoordCGIC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHyperSumCGI = CoordDMRCHyperSumCGI + out[[5]]
      }
      if (length(CoordDMRCHypo)>0)
      {
        
      N<-nrow(CoordDMRCHypo)
      M<-nrow(CoordCGIC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypo), as.matrix(CoordCGIC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHypoSumCGI = CoordDMRCHypoSumCGI + out[[5]]
      }
      
      ###########################################
      #### Searching Intergenic and DMRs Overlap #######
      ###########################################
      if (length(CoordDMRCHyper)>0)
      {
        
        N<-nrow(CoordDMRCHyper)
        M<-nrow(CoordIntergenicC)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyper), as.matrix(CoordIntergenicC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHyperSumIntergenic = CoordDMRCHyperSumIntergenic + out[[5]]
      }
      if (length(CoordDMRCHypo)>0)
      {
        
        N<-nrow(CoordDMRCHypo)
        M<-nrow(CoordIntergenicC)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypo), as.matrix(CoordIntergenicC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHypoSumIntergenic = CoordDMRCHypoSumIntergenic + out[[5]]
      }
      ###########################################
      #### Searching Enhancer and DMRs Overlap #######
      ###########################################
      if (length(CoordDMRCHyper)>0)
      {
        
      N<-nrow(CoordDMRCHyper)
      M<-nrow(CoordEnhancerC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyper), as.matrix(CoordEnhancerC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHyperSumEnhancer = CoordDMRCHyperSumEnhancer + out[[5]]
      }
      
      if (length(CoordDMRCHypo)>0)
      {
        
      N<-nrow(CoordDMRCHypo)
      M<-nrow(CoordEnhancerC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypo), as.matrix(CoordEnhancerC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHypoSumEnhancer = CoordDMRCHypoSumEnhancer + out[[5]]
      }
      
      #########################################################
      ##### Calculating Total DMR in each GenCode Feature #####
      #########################################################
      
      
      
      for (kk in 1:length(BaseNameFeature))
      {
        ind2Feat<-which(FeatureGenCodeSubC==BaseNameFeature[kk])
        
        CoordGenCodeC2Feat<-CoordGenCodeC[ind2Feat,]
        
        #### Searching GenCode and DMRs Overlap #######
        if (length(CoordDMRCHyper)>0)
        {
          
        N<-nrow(CoordDMRCHyper)
        M<-nrow(CoordGenCodeC2Feat)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyper), as.matrix(CoordGenCodeC2Feat), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHyperSumGencode[kk] = CoordDMRCHyperSumGencode[kk]+ out[[5]]
        }
        
        #### Searching GenCode and DMRs Overlap #######
        if (length(CoordDMRCHypo)>0)
        {
          
        N<-nrow(CoordDMRCHypo)
        M<-nrow(CoordGenCodeC2Feat)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypo), as.matrix(CoordGenCodeC2Feat), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHypoSumGencode[kk] = CoordDMRCHypoSumGencode[kk] + out[[5]]
        }
      }
      
      
    }
    
    
    if (AnalysisClass=="Entropy")
    {
      CoordDMRCHigh<-rbind(CoordDMRC[which(DeltaBetaDMRC>0 & DeltaEntropyDMRC>=EntropyThr),])
      CoordDMRCMid<-rbind(CoordDMRC[which(DeltaBetaDMRC>0 & DeltaEntropyDMRC> -EntropyThr & DeltaEntropyDMRC< EntropyThr),])
      CoordDMRCLow<-rbind(CoordDMRC[which(DeltaBetaDMRC>0 & DeltaEntropyDMRC<= -EntropyThr),])
      
      if (length(CoordDMRCHigh)>0)
      {
        CoordDMRCHighSumTotal<- CoordDMRCHighSumTotal + sum(CoordDMRCHigh[,2]-CoordDMRCHigh[,1])
      }
      if (length(CoordDMRCMid)>0)
      {
        CoordDMRCMidSumTotal<- CoordDMRCMidSumTotal + sum(CoordDMRCMid[,2]-CoordDMRCMid[,1])
      }
      if (length(CoordDMRCLow)>0)
      {
        CoordDMRCLowSumTotal<- CoordDMRCLowSumTotal + sum(CoordDMRCLow[,2]-CoordDMRCLow[,1])
      }
      
      ###########################################
      #### Searching CGI and DMRs Overlap #######
      ###########################################
      if (length(CoordDMRCHigh)>0)
      {
        
      N<-nrow(CoordDMRCHigh)
      M<-nrow(CoordCGIC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHigh), as.matrix(CoordCGIC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHighSumCGI = CoordDMRCHighSumCGI + out[[5]]
      }
      
      if (length(CoordDMRCMid)>0)
      {
        
      N<-nrow(CoordDMRCMid)
      M<-nrow(CoordCGIC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCMid), as.matrix(CoordCGIC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCMidSumCGI = CoordDMRCMidSumCGI + out[[5]]
      }
      
      if (length(CoordDMRCLow)>0)
      {
        
      N<-nrow(CoordDMRCLow)
      M<-nrow(CoordCGIC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCLow), as.matrix(CoordCGIC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCLowSumCGI = CoordDMRCLowSumCGI + out[[5]]
      }
      
      
      ###########################################
      #### Searching Intergenic and DMRs Overlap #######
      ###########################################
      if (length(CoordDMRCHigh)>0)
      {
        
        N<-nrow(CoordDMRCHigh)
        M<-nrow(CoordIntergenicC)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHigh), as.matrix(CoordIntergenicC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHighSumIntergenic = CoordDMRCHighSumIntergenic + out[[5]]
      }
      
      if (length(CoordDMRCMid)>0)
      {
        
        N<-nrow(CoordDMRCMid)
        M<-nrow(CoordIntergenicC)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCMid), as.matrix(CoordIntergenicC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCMidSumIntergenic = CoordDMRCMidSumIntergenic + out[[5]]
      }
      
      if (length(CoordDMRCLow)>0)
      {
        
        N<-nrow(CoordDMRCLow)
        M<-nrow(CoordIntergenicC)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCLow), as.matrix(CoordIntergenicC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCLowSumIntergenic = CoordDMRCLowSumIntergenic + out[[5]]
      }
      
      
      ###########################################
      #### Searching Enhancer and DMRs Overlap #######
      ###########################################
      if (length(CoordDMRCHigh)>0)
      {
        
      N<-nrow(CoordDMRCHigh)
      M<-nrow(CoordEnhancerC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHigh), as.matrix(CoordEnhancerC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHighSumEnhancer = CoordDMRCHighSumEnhancer + out[[5]]
      }
      
      if (length(CoordDMRCMid)>0)
      {
        
      N<-nrow(CoordDMRCMid)
      M<-nrow(CoordEnhancerC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCMid), as.matrix(CoordEnhancerC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCMidSumEnhancer = CoordDMRCMidSumEnhancer + out[[5]]
      }
      
      if (length(CoordDMRCLow)>0)
      {
        
      N<-nrow(CoordDMRCLow)
      M<-nrow(CoordEnhancerC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCLow), as.matrix(CoordEnhancerC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCLowSumEnhancer = CoordDMRCLowSumEnhancer + out[[5]]
      }
      
      #########################################################
      ##### Calculating Total DMR in each GenCode Feature #####
      #########################################################
      
      
      
      for (kk in 1:length(BaseNameFeature))
      {
        ind2Feat<-which(FeatureGenCodeSubC==BaseNameFeature[kk])
        
        CoordGenCodeC2Feat<-CoordGenCodeC[ind2Feat,]
        
        #### Searching GenCode and DMRs Overlap #######
        if (length(CoordDMRCHigh)>0)
        {
          
        N<-nrow(CoordDMRCHigh)
        M<-nrow(CoordGenCodeC2Feat)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHigh), as.matrix(CoordGenCodeC2Feat), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHighSumGencode[kk] = CoordDMRCHighSumGencode[kk]+ out[[5]]
        }
        
        #### Searching GenCode and DMRs Overlap #######
        if (length(CoordDMRCMid)>0)
        {
          
        N<-nrow(CoordDMRCMid)
        M<-nrow(CoordGenCodeC2Feat)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCMid), as.matrix(CoordGenCodeC2Feat), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCMidSumGencode[kk] = CoordDMRCMidSumGencode[kk] + out[[5]]
        
        }
        
        #### Searching GenCode and DMRs Overlap #######
        if (length(CoordDMRCLow)>0)
        {
          
        N<-nrow(CoordDMRCLow)
        M<-nrow(CoordGenCodeC2Feat)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCLow), as.matrix(CoordGenCodeC2Feat), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCLowSumGencode[kk] = CoordDMRCLowSumGencode[kk] + out[[5]]
        }
      }
    }
    
  }
  
  
  
  
  
  
  if (AnalysisClass=="Beta")
  {
    
    TotalSizeHyper<-c(CoordDMRCHyperSumIntergenic,CoordDMRCHyperSumGencode,CoordDMRCHyperSumCGI,CoordDMRCHyperSumEnhancer)
    TotalSizeHypo<-c(CoordDMRCHypoSumIntergenic,CoordDMRCHypoSumGencode,CoordDMRCHypoSumCGI,CoordDMRCHypoSumEnhancer)
    MatrixSizeOut<-rbind(TotalSizeHyper,
                         TotalSizeHypo)
    rownames(MatrixSizeOut)<-c("Hyper","Hypo")
    colnames(MatrixSizeOut)<-FeatureNameOut
  }
  
  
  if (AnalysisClass=="Entropy")
  {
    
    TotalSizeHigh<-c(CoordDMRCHighSumIntergenic,CoordDMRCHighSumGencode,CoordDMRCHighSumCGI,CoordDMRCHighSumEnhancer)
    TotalSizeMid<-c(CoordDMRCMidSumIntergenic,CoordDMRCMidSumGencode,CoordDMRCMidSumCGI,CoordDMRCMidSumEnhancer)
    TotalSizeLow<-c(CoordDMRCLowSumIntergenic,CoordDMRCLowSumGencode,CoordDMRCLowSumCGI,CoordDMRCLowSumEnhancer)
    MatrixSizeOut<-rbind(TotalSizeHigh,
                         TotalSizeMid,
                         TotalSizeLow)
    rownames(MatrixSizeOut)<-c("High","Mid","Low")
    colnames(MatrixSizeOut)<-FeatureNameOut
  }
  
  if (AnalysisClass=="All")
  {
    
    TotalSizeHyperHigh<-c(CoordDMRCHyperHighSumIntergenic,CoordDMRCHyperHighSumGencode,CoordDMRCHyperHighSumCGI,CoordDMRCHyperHighSumEnhancer)
    TotalSizeHyperMid<-c(CoordDMRCHyperMidSumIntergenic,CoordDMRCHyperMidSumGencode,CoordDMRCHyperMidSumCGI,CoordDMRCHyperMidSumEnhancer)
    TotalSizeHyperLow<-c(CoordDMRCHyperLowSumIntergenic,CoordDMRCHyperLowSumGencode,CoordDMRCHyperLowSumCGI,CoordDMRCHyperLowSumEnhancer)
    TotalSizeHypoHigh<-c(CoordDMRCHypoHighSumIntergenic,CoordDMRCHypoHighSumGencode,CoordDMRCHypoHighSumCGI,CoordDMRCHypoHighSumEnhancer)
    TotalSizeHypoMid<-c(CoordDMRCHypoMidSumIntergenic,CoordDMRCHypoMidSumGencode,CoordDMRCHypoMidSumCGI,CoordDMRCHypoMidSumEnhancer)
    TotalSizeHypoLow<-c(CoordDMRCHypoLowSumIntergenic,CoordDMRCHypoLowSumGencode,CoordDMRCHypoLowSumCGI,CoordDMRCHypoLowSumEnhancer)
    
    
    MatrixSizeOut<-rbind(TotalSizeHyperHigh,
                         TotalSizeHyperMid,
                         TotalSizeHyperLow,
                         TotalSizeHypoHigh,
                         TotalSizeHypoMid,
                         TotalSizeHypoLow)
    rownames(MatrixSizeOut)<-c("HyperHigh","HyperMid","HyperLow","HypoHigh","HypoMid","HypoLow")
    colnames(MatrixSizeOut)<-FeatureNameOut
  }
  
  
  MatrixSizeOut
  
}
