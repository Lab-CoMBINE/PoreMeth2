#' Functions for Calculating the size of DMR in each genic feature
#'
#' description description description description description description
#'
#' @param TableInDMR
#' @param Assembly ,
#' @param BetaThr numeric,
#' @param PValueThr numeric,
#'
#' @return what what what what what what what what what
#'
#' @export


PoreMeth1DMRStatistics <- function(TableInDMR,Assembly="hg19",
                                   BetaThr=0.2,
                                   PValueThr=0.05)
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
    
  
  #dyn.load("/home/alberto/NanoporeMethylation/AnnotationTool/OverlapLibrary.so")
  
  
  ChrDMRInU<-unique(ChrDMR)
  for (i in 1:length(ChrDMRInU))
  {
    
    
    #### Selecting Chr Coordinates for DMRs ####
    indDMRC<-which(ChrDMR==ChrDMRInU[i])
    StartDMRC<-StartDMR[indDMRC]
    EndDMRC<-EndDMR[indDMRC]
    DeltaBetaDMRC<-DeltaBetaDMR[indDMRC]
    
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
    

    
    
    
      CoordDMRCHyper<-CoordDMRC[which(DeltaBetaDMRC>0),]
      CoordDMRCHypo<-CoordDMRC[which(DeltaBetaDMRC<0),]
      
      
      
      CoordDMRCHyperSumTotal<- CoordDMRCHyperSumTotal + sum(CoordDMRCHyper[,2]-CoordDMRCHyper[,1])
      CoordDMRCHypoSumTotal<- CoordDMRCHypoSumTotal + sum(CoordDMRCHypo[,2]-CoordDMRCHypo[,1])
      
      
      ###########################################
      #### Searching Intergenic and DMRs Overlap #######
      ###########################################
      
      N<-nrow(CoordDMRCHyper)
      M<-nrow(CoordIntergenicC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyper), as.matrix(CoordIntergenicC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHyperSumIntergenic = CoordDMRCHyperSumIntergenic + out[[5]]
      
      
      N<-nrow(CoordDMRCHypo)
      M<-nrow(CoordIntergenicC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypo), as.matrix(CoordIntergenicC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHypoSumIntergenic = CoordDMRCHypoSumIntergenic + out[[5]]
      
      
      
      
      
      ###########################################
      #### Searching CGI and DMRs Overlap #######
      ###########################################
      
      N<-nrow(CoordDMRCHyper)
      M<-nrow(CoordCGIC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyper), as.matrix(CoordCGIC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHyperSumCGI = CoordDMRCHyperSumCGI + out[[5]]
      
      N<-nrow(CoordDMRCHypo)
      M<-nrow(CoordCGIC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypo), as.matrix(CoordCGIC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHypoSumCGI = CoordDMRCHypoSumCGI + out[[5]]
      
      
      
      ###########################################
      #### Searching Enhancer and DMRs Overlap #######
      ###########################################
      
      N<-nrow(CoordDMRCHyper)
      M<-nrow(CoordEnhancerC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyper), as.matrix(CoordEnhancerC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHyperSumEnhancer = CoordDMRCHyperSumEnhancer + out[[5]]
      
      N<-nrow(CoordDMRCHypo)
      M<-nrow(CoordEnhancerC)
      OS<-0
      
      out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypo), as.matrix(CoordEnhancerC), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
      
      
      CoordDMRCHypoSumEnhancer = CoordDMRCHypoSumEnhancer + out[[5]]
      
      
      #########################################################
      ##### Calculating Total DMR in each GenCode Feature #####
      #########################################################
      
      
      
      for (kk in 1:length(BaseNameFeature))
      {
        ind2Feat<-which(FeatureGenCodeSubC==BaseNameFeature[kk])
        
        CoordGenCodeC2Feat<-CoordGenCodeC[ind2Feat,]
        
        #### Searching GenCode and DMRs Overlap #######
        N<-nrow(CoordDMRCHyper)
        M<-nrow(CoordGenCodeC2Feat)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHyper), as.matrix(CoordGenCodeC2Feat), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHyperSumGencode[kk] = CoordDMRCHyperSumGencode[kk]+ out[[5]]
        
        
        #### Searching GenCode and DMRs Overlap #######
        N<-nrow(CoordDMRCHypo)
        M<-nrow(CoordGenCodeC2Feat)
        OS<-0
        
        out = .Fortran("matrixoverlapsum", as.matrix(CoordDMRCHypo), as.matrix(CoordGenCodeC2Feat), as.integer(N), as.integer(M), as.numeric(OS), PACKAGE = "PoreMeth2")
        
        
        CoordDMRCHypoSumGencode[kk] = CoordDMRCHypoSumGencode[kk] + out[[5]]
        
      }
      
      
    

  
  
  
  
  

    TotalSizeHyper<-c(CoordDMRCHyperSumIntergenic,CoordDMRCHyperSumGencode,CoordDMRCHyperSumCGI,CoordDMRCHyperSumEnhancer)
    TotalSizeHypo<-c(CoordDMRCHypoSumIntergenic,CoordDMRCHypoSumGencode,CoordDMRCHypoSumCGI,CoordDMRCHypoSumEnhancer)
    MatrixSizeOut<-rbind(TotalSizeHyper,
                         TotalSizeHypo)
    rownames(MatrixSizeOut)<-c("Hyper","Hypo")
    colnames(MatrixSizeOut)<-FeatureNameOut
  }
  
  
  
  
  MatrixSizeOut
  
}

