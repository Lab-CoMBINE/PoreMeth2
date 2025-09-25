#' Main Function for DMRs Annotation
#'
#' Function for genomic and regulatory annotation of results from PoreMeth2DMR or PoreMeth1DMR
#'
#' @param TableMethIn Output from PoreMeth2DMR (or PoreMeth1DMR), data.frame.
#' @param NumProc Optional argument for the number of cores to use in parallel, numeric.
#' @param AnnotationType Optional argument that specifies whether to annotate DMRs on genic elements only (`"Genes"`) or genic elements and regulatory features (`"GenesReg"`), character.
#' @param Assembly Optional parameter that specifies the reference version to use for annotation (`"hg19"` or `"hg38"`), character.
#'
#' @return data.frame with DMRs annotated to genomic and regulatory elements.
#'
#' @export

PoreMethAnnotate <- function(TableMethIn,FileOut,NumProc=5,AnnotationType="Genes",Assembly="hg19")
{
  Out <- list()
  PathDBIn <- paste0(.libPaths()[1], "/PoreMeth2/data/")
  
  
  cl <- parallel::makePSOCKcluster(NumProc)
  parallel::setDefaultCluster(cl)
  parallel::clusterExport(NULL, c("RegionAnnotateCompact"))
  
  
  ############################
  ##### Loading Meth File #####
  ############################
  
  
  StringChr<-TableMethIn$chr[1]
  
  
  if (length(grep("chr",StringChr))!=0)
  {
    ChrMethIn<-TableMethIn$chr
  }
  
  if (length(grep("chr",StringChr))==0)
  {
    ChrMethIn<-paste("chr",TableMethIn$chr,sep="")
  }
  
  
  StartMethIn<-as.numeric(TableMethIn$start)
  EndMethIn<-as.numeric(TableMethIn$end)
  
  NFieldMethIn<-ncol(TableMethIn)
  FieldNamesMethIn<-colnames(TableMethIn)
  
  MatrixMethResults<-TableMethIn[,4:NFieldMethIn]
  
  ChrMethInU<-unique(ChrMethIn)
  
  
  
  #################################################
  #################################################
  ###### Loading Genomic Features Databases #######
  #################################################
  #################################################
  
  ########################
  ##### Flat GenCode #####
  ########################
  
  FileGenCodeIn<-file.path(PathDBIn,paste("GencodeTable_",Assembly,".rds",sep=""))
  TableGenCodeIn<-readRDS(FileGenCodeIn)
  
  ChrGenCode<-TableGenCodeIn$Chr
  StartGenCode<-as.numeric(TableGenCodeIn$Start)
  EndGenCode<-as.numeric(TableGenCodeIn$End)
  FeatureGenCode<-TableGenCodeIn$`Genic Element`
  StrandGenCode<-TableGenCodeIn$Strand
  SymbolGenCode<-TableGenCodeIn$`Gene Symbol`
  TypeGenCode<-TableGenCodeIn$`Transcript Type`
  
  ######################################
  ######################################
  #### Loading Regulatory Databases ####
  ######################################
  ######################################
  
  if (AnnotationType=="GenesReg")
  {
    ############################
    ##### File CpG Islands #####
    ############################
    
    FileCGIIn<-file.path(PathDBIn,paste("CGIAnno_",Assembly,".rds",sep=""))
    TableCGIIn<-readRDS(FileCGIIn)
    
    ChrCGI<-as.character(TableCGIIn$Chr)
    StartCGI<-as.numeric(TableCGIIn$Start)
    EndCGI<-as.numeric(TableCGIIn$End)
    NameCGI<-as.character(TableCGIIn$FeatureName)
    
    
    ############################
    ######## File DNASE ########
    ############################
    
    FileDNASEIn<-file.path(PathDBIn,paste("DNASEAnno_",Assembly,".rds",sep=""))
    TableDNASEIn<-readRDS(FileDNASEIn)
    
    ChrDNASE<-as.character(TableDNASEIn$Chr)
    StartDNASE<-as.numeric(TableDNASEIn$Start)
    EndDNASE<-as.numeric(TableDNASEIn$End)
    NameDNASE<-as.character(TableDNASEIn$FeatureName)
    
    ############################
    ######## File TFBS  ########
    ############################
    
    FileTFBSIn<-file.path(PathDBIn,paste("TFBSAnno_",Assembly,".rds",sep=""))
    TableTFBSIn<-readRDS(FileTFBSIn)
    
    ChrTFBS<-as.character(TableTFBSIn$Chr)
    StartTFBS<-as.numeric(TableTFBSIn$Start)
    EndTFBS<-as.numeric(TableTFBSIn$End)
    NameTFBS<-as.character(TableTFBSIn$FeatureName)
    
    
    #################################
    ######## File Enhancers  ########
    #################################
    
    FileEnhancerIn<-file.path(PathDBIn,paste("EnhancerAnno_",Assembly,".rds",sep=""))
    TableEnhancerIn<-readRDS(FileEnhancerIn)
    
    ChrEnhancer<-as.character(TableEnhancerIn$Chr)
    StartEnhancer<-as.numeric(TableEnhancerIn$Start)
    EndEnhancer<-as.numeric(TableEnhancerIn$End)
    NameEnhancer<-as.character(TableEnhancerIn$FeatureName)
    
  }
  
    
  
  
  
  
  ###################################
  ####### Generating OutPut #########
  ###################################
  
  
  GenCodeAnnoName<-c("chr.GenCode",
                     "start.GenCode",
                     "end.GenCode",
                     "feature.GenCode",
                     "strand.GenCode",
                     "symbol.GenCode",
                     "type.GenCode",
                     "chr.GenCode.overlap",
                     "start.GenCode.overlap",
                     "end.GenCode.overlap",
                     "ratio1.GenCode.overlap",
                     "ratio2.GenCode.overlap")
  
  
  DNASEAnnoName<-c("name.DNASE","chr.DNASE",
                   "start.DNASE","end.DNASE",
                   "chr.DNASE.overlap",
                   "start.DNASE.overlap",
                   "end.DNASE.overlap")
  
  CGIAnnoName<-c("name.CGI","chr.CGI",
                 "start.CGI","end.CGI",
                 "chr.CGI.overlap",
                 "start.CGI.overlap",
                 "end.CGI.overlap")
  
  TFBSAnnoName<-c("name.TFBS","chr.TFBS",
                  "start.TFBS","end.TFBS",
                  "chr.TFBS.overlap",
                  "start.TFBS.overlap",
                  "end.TFBS.overlap")
  
  EnhancerAnnoName<-c("name.Enhancer","chr.Enhancer",
                      "start.Enhancer","end.Enhancer",
                      "chr.Enhancer.overlap",
                      "start.Enhancer.overlap",
                      "end.Enhancer.overlap")
  
  
  
  if (AnnotationType=="Genes")
  {
    
    
    
    TotalNameOut<-c(FieldNamesMethIn,
                    GenCodeAnnoName)
    
  }
  
  
  if (AnnotationType=="GenesReg")
  {
    
    TotalNameOut<-c(FieldNamesMethIn,
                    GenCodeAnnoName,
                    CGIAnnoName,
                    EnhancerAnnoName,
                    DNASEAnnoName,
                    TFBSAnnoName)
  }
  
    Out <- append(Out, list(rbind(TotalNameOut)))

  ###############################################
  ###############################################
  ###### Analysis Chromosome by Chromosome ######
  ###############################################
  ###############################################
  
  for (i in 1:length(ChrMethInU))
  {
    ############################################
    #### Selecting Chr Coordinates for Meths ####
    ############################################
    
    
    indMethC<-which(ChrMethIn==ChrMethInU[i])
    StartMethInC<-StartMethIn[indMethC]
    EndMethInC<-EndMethIn[indMethC]
    MatrixMethResultsC<-MatrixMethResults[indMethC,]
    
    ###############################################
    #### Selecting Chr Coordinates for GenCode ####
    ###############################################
    
    indGenCodeC<-which(ChrGenCode==ChrMethInU[i])
    StartGenCodeC<-StartGenCode[indGenCodeC]
    EndGenCodeC<-EndGenCode[indGenCodeC]
    FeatureGenCodeC<-FeatureGenCode[indGenCodeC]
    StrandGenCodeC<-StrandGenCode[indGenCodeC]
    SymbolGenCodeC<-SymbolGenCode[indGenCodeC]
    TypeGenCodeC<-TypeGenCode[indGenCodeC]
    
    CoordGenCodeC<-cbind(StartGenCodeC,EndGenCodeC)
    
    if (AnnotationType=="GenesReg")
    {
      ###########################################
      #### Selecting Chr Coordinates for CGI ####
      ###########################################
      
      indCGIC<-which(ChrCGI==ChrMethInU[i])
      StartCGIC<-StartCGI[indCGIC]
      EndCGIC<-EndCGI[indCGIC]
      NameCGIC<-NameCGI[indCGIC]
      
      CoordCGIC<-cbind(StartCGIC,EndCGIC)
      
      
      ################################################
      #### Selecting Chr Coordinates for Enhancer ####
      ################################################
      
      indEnhancerC<-which(ChrEnhancer==ChrMethInU[i])
      StartEnhancerC<-StartEnhancer[indEnhancerC]
      EndEnhancerC<-EndEnhancer[indEnhancerC]
      NameEnhancerC<-NameEnhancer[indEnhancerC]
      
      CoordEnhancerC<-cbind(StartEnhancerC,EndEnhancerC)
      
      
      #############################################
      #### Selecting Chr Coordinates for DNASE ####
      #############################################
      
      indDNASEC<-which(ChrDNASE==ChrMethInU[i])
      StartDNASEC<-StartDNASE[indDNASEC]
      EndDNASEC<-EndDNASE[indDNASEC]
      NameDNASEC<-NameDNASE[indDNASEC]
      
      
      CoordDNASEC<-cbind(StartDNASEC,EndDNASEC)
      
      
      ###########################################
      #### Selecting Chr Coordinates for TFBS ####
      ###########################################
      
      indTFBSC<-which(ChrTFBS==ChrMethInU[i])
      StartTFBSC<-StartTFBS[indTFBSC]
      EndTFBSC<-EndTFBS[indTFBSC]
      NameTFBSC<-NameTFBS[indTFBSC]
      
      
      CoordTFBSC<-cbind(StartTFBSC,EndTFBSC)
    }
    
    if (AnnotationType=="Genes")
    {
      CoordCGIC<-c()
      NameCGIC<-c()
      CoordEnhancerC<-c()
      NameEnhancerC<-c()
      CoordDNASEC<-c()
      NameDNASEC<-c()
      CoordTFBSC<-c()
      NameTFBSC<-c()
    }
    #################################
    #################################
    ###### Parallel Annotation ######
    #################################
    #################################
    
    
    
    RegionAnnotateCompactOut<-parallel::parLapply(cl,c(1:length(StartMethInC)),RegionAnnotateCompact,ChrMethInU[i],StartMethInC,EndMethInC,
                                        MatrixMethResultsC,
                                        CoordGenCodeC,
                                        CoordCGIC,
                                        NameCGIC,
                                        CoordEnhancerC,
                                        NameEnhancerC,
                                        CoordDNASEC,
                                        NameDNASEC,
                                        CoordTFBSC,
                                        NameTFBSC,
                                        FeatureGenCodeC,
                                        StrandGenCodeC,
                                        SymbolGenCodeC,
                                        TypeGenCodeC,
                                        AnnotationType)
    
    
    NColOut<-ncol(RegionAnnotateCompactOut[[1]])
    

    for (kk in 1:length(StartMethInC))
    {
      Out <- append(Out, list(RegionAnnotateCompactOut[[kk]]))
    }

}
  OutDF <- data.frame(do.call("rbind", Out[-1]), row.names = NULL)
  colnames(OutDF) <- Out[[1]]
  return(OutDF)
  parallel::stopCluster(cl)

}
