#' Parallel Function for Compact Annotation
#' @NoRd
#' @export


RegionAnnotateCompact <- function(j,ChrMethInUSel,StartMethInC,EndMethInC,
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
{
  
  
  CoordMeth2Test<-c(StartMethInC[j],EndMethInC[j])
  MatrixMeth2Test<-MatrixMethResultsC[j,]
  
  
  
  ################################################
  ################################################
  ##### Searching for GenCode Features in Meth #####
  ################################################
  ################################################
  
  
  ##### Loading Fortran Library ######

 PathSrc <- paste0(.libPaths()[1], "/PoreMeth2/libs/")
 dyn.load(paste0(PathSrc, "PoreMeth2.so"))
  
  
  M<-nrow(CoordGenCodeC)
  vecOverlapA = rep(0,M)
  vecOverlapB = rep(0,M)
  
  CoordOverlap<-cbind(rep(0,M),rep(0,M))
  
  out = .Fortran("matrixcompare", as.vector(CoordMeth2Test), as.matrix(CoordGenCodeC), as.integer(M), as.matrix(CoordOverlap), as.numeric(vecOverlapA), as.numeric(vecOverlapB), PACKAGE = "PoreMeth2")
  
  
  CoordOverlap<-out[[4]]
  vecOverlapA = out[[5]]
  vecOverlapB = out[[6]]
  
  indMethOverlap<-which(vecOverlapA!=0)
  
  if (AnnotationType=="Genes")
  {
    if (length(indMethOverlap)>1)
    {
      
      ##################################
      ##### Creating Matrix GenCode #####
      ##################################
      
      CoordMeth2Test2Save<-matrix(rep((c(ChrMethInUSel,CoordMeth2Test)),each=length(indMethOverlap)),nrow=length(indMethOverlap))
      
      
      
      MatrixMeth2Test2Save<-matrix(rep(as.character(MatrixMeth2Test),each=length(indMethOverlap)),nrow=length(indMethOverlap))
      
      MatrixGenCodeOverlap<-cbind(rep(ChrMethInUSel,length(indMethOverlap)),
                                  CoordGenCodeC[indMethOverlap,1],
                                  CoordGenCodeC[indMethOverlap,2],
                                  FeatureGenCodeC[indMethOverlap],
                                  StrandGenCodeC[indMethOverlap],
                                  SymbolGenCodeC[indMethOverlap],
                                  TypeGenCodeC[indMethOverlap],
                                  rep(ChrMethInUSel,length(indMethOverlap)),
                                  CoordOverlap[indMethOverlap,],
                                  vecOverlapA[indMethOverlap],
                                  vecOverlapB[indMethOverlap])
      
      
      
      MatrixOut<-cbind(CoordMeth2Test2Save,
                       MatrixMeth2Test2Save,
                       MatrixGenCodeOverlap)
    }
    
    if (length(indMethOverlap)==1)
    {
      
      ##################################
      ##### Creating Matrix GenCode #####
      ##################################
      
      CoordMeth2Test2Save<-matrix(rep((c(ChrMethInUSel,CoordMeth2Test)),each=length(indMethOverlap)),nrow=length(indMethOverlap))
      
      
      
      MatrixMeth2Test2Save<-matrix(rep(as.character(MatrixMeth2Test),each=length(indMethOverlap)),nrow=length(indMethOverlap))
      
      MatrixGenCodeOverlap<-c(rep(ChrMethInUSel,length(indMethOverlap)),
                              CoordGenCodeC[indMethOverlap,1],
                              CoordGenCodeC[indMethOverlap,2],
                              FeatureGenCodeC[indMethOverlap],
                              StrandGenCodeC[indMethOverlap],
                              SymbolGenCodeC[indMethOverlap],
                              TypeGenCodeC[indMethOverlap],
                              rep(ChrMethInUSel,length(indMethOverlap)),
                              CoordOverlap[indMethOverlap,],
                              vecOverlapA[indMethOverlap],
                              vecOverlapB[indMethOverlap])
      
      
      
      MatrixOut<-rbind(c(CoordMeth2Test2Save,
                         MatrixMeth2Test2Save,
                         MatrixGenCodeOverlap))
    }
    if (length(indMethOverlap)==0)
    {
      ##################################
      ##### Creating Matrix GenCode #####
      ##################################
      
      CoordMeth2Test2Save<-matrix(rep((c(ChrMethInUSel,CoordMeth2Test)),each=1),nrow=1)
      
      
      
      MatrixMeth2Test2Save<-matrix(rep(as.character(MatrixMeth2Test),each=1),nrow=1)
      
      MatrixGenCodeOverlap<-c("*",
                              "*",
                              "*",
                              "InterGenic",
                              "*",
                              "*",
                              "*",
                              "*",
                              "*",
                              "*",
                              "*",
                              "*")
      
      
      
      
      MatrixOut<-rbind(c(CoordMeth2Test2Save,
                         MatrixMeth2Test2Save,
                         MatrixGenCodeOverlap))
      
    }
    
  }
  
  
  if (AnnotationType=="GenesReg")
  {
    if (length(indMethOverlap)>1)
    {
      ###################################################################
      #### Annotating Genic Feature Overlap with Regulatory Elements ####
      ###################################################################
      
      
      CoordOverlapMeth<-CoordOverlap[indMethOverlap,]
      
      
      #####################################################
      ##### Search for Meths and GenCode in CpG Islands #####
      #####################################################
      
      M<-nrow(CoordOverlapMeth)
      vecOverlapCGI = rep(0,M)
      
      
      CoordOverlapCGI<-cbind(rep(0,M),rep(0,M))
      
      out = .Fortran("matrixcompare3", as.matrix(CoordOverlapMeth), as.matrix(CoordCGIC), as.integer(M), as.integer(nrow(CoordCGIC)), as.numeric(vecOverlapCGI), as.matrix(CoordOverlapCGI), PACKAGE = "PoreMeth2")
      
      
      CoordOverlapCGI<-out[[6]]
      vecOverlapCGI = out[[5]]
      
      
      indCGIOverlap<-which(vecOverlapCGI!=0)
      MatrixOverlapCGI<-cbind(rep("*",length(vecOverlapCGI)),
                              rep("*",length(vecOverlapCGI)),
                              rep("*",length(vecOverlapCGI)),
                              rep("*",length(vecOverlapCGI)),
                              rep("*",length(vecOverlapCGI)),
                              CoordOverlapCGI)
      if (length(indCGIOverlap)!=0)
      {
        MatrixOverlapCGI[indCGIOverlap,1]<-NameCGIC[vecOverlapCGI[indCGIOverlap]]
        MatrixOverlapCGI[indCGIOverlap,2]<-rep(ChrMethInUSel,length(indCGIOverlap))
        MatrixOverlapCGI[indCGIOverlap,3]<-CoordCGIC[vecOverlapCGI[indCGIOverlap],1]
        MatrixOverlapCGI[indCGIOverlap,4]<-CoordCGIC[vecOverlapCGI[indCGIOverlap],2]
        MatrixOverlapCGI[indCGIOverlap,5]<-rep(ChrMethInUSel,length(indCGIOverlap))
      }
      
      
      #####################################################
      #####  Search for Meths and GenCode in Enhancers  #####
      #####################################################
      
      M<-nrow(CoordOverlapMeth)
      vecOverlapEnhancer = rep(0,M)
      
      
      CoordOverlapEnhancer<-cbind(rep(0,M),rep(0,M))
      
      out = .Fortran("matrixcompare3", as.matrix(CoordOverlapMeth), as.matrix(CoordEnhancerC), as.integer(M), as.integer(nrow(CoordEnhancerC)), as.numeric(vecOverlapEnhancer), as.matrix(CoordOverlapEnhancer), PACKAGE = "PoreMeth2")
      
      
      CoordOverlapEnhancer<-out[[6]]
      vecOverlapEnhancer = out[[5]]
      
      
      indEnhancerOverlap<-which(vecOverlapEnhancer!=0)
      MatrixOverlapEnhancer<-cbind(rep("*",length(vecOverlapEnhancer)),
                                   rep("*",length(vecOverlapEnhancer)),
                                   rep("*",length(vecOverlapEnhancer)),
                                   rep("*",length(vecOverlapEnhancer)),
                                   rep("*",length(vecOverlapEnhancer)),
                                   CoordOverlapEnhancer)
      if (length(indEnhancerOverlap)!=0)
      {
        MatrixOverlapEnhancer[indEnhancerOverlap,1]<-NameEnhancerC[vecOverlapEnhancer[indEnhancerOverlap]]
        MatrixOverlapEnhancer[indEnhancerOverlap,2]<-rep(ChrMethInUSel,length(indEnhancerOverlap))
        MatrixOverlapEnhancer[indEnhancerOverlap,3]<-CoordEnhancerC[vecOverlapEnhancer[indEnhancerOverlap],1]
        MatrixOverlapEnhancer[indEnhancerOverlap,4]<-CoordEnhancerC[vecOverlapEnhancer[indEnhancerOverlap],2]
        MatrixOverlapEnhancer[indEnhancerOverlap,5]<-rep(ChrMethInUSel,length(indEnhancerOverlap))
      }
      
      ###############################################
      ##### Search for Meths and GenCode in DNASE #####
      ###############################################
      
      M<-nrow(CoordOverlapMeth)
      vecOverlapDNASE = rep(0,M)
      
      
      CoordOverlapDNASE<-cbind(rep(0,M),rep(0,M))
      
      out = .Fortran("matrixcompare3", as.matrix(CoordOverlapMeth), as.matrix(CoordDNASEC), as.integer(M), as.integer(nrow(CoordDNASEC)), as.numeric(vecOverlapDNASE), as.matrix(CoordOverlapDNASE), PACKAGE = "PoreMeth2")
      
      
      CoordOverlapDNASE<-out[[6]]
      vecOverlapDNASE = out[[5]]
      
      
      indDNASEOverlap<-which(vecOverlapDNASE!=0)
      MatrixOverlapDNASE<-cbind(rep("*",length(vecOverlapDNASE)),
                                rep("*",length(vecOverlapDNASE)),
                                rep("*",length(vecOverlapDNASE)),
                                rep("*",length(vecOverlapDNASE)),
                                rep("*",length(vecOverlapDNASE)),
                                CoordOverlapDNASE)
      if (length(indDNASEOverlap)!=0)
      {
        MatrixOverlapDNASE[indDNASEOverlap,1]<-NameDNASEC[vecOverlapDNASE[indDNASEOverlap]] 
        MatrixOverlapDNASE[indDNASEOverlap,2]<-rep(ChrMethInUSel,length(indDNASEOverlap))
        MatrixOverlapDNASE[indDNASEOverlap,3]<-CoordDNASEC[vecOverlapDNASE[indDNASEOverlap],1]
        MatrixOverlapDNASE[indDNASEOverlap,4]<-CoordDNASEC[vecOverlapDNASE[indDNASEOverlap],2]
        MatrixOverlapDNASE[indDNASEOverlap,5]<-rep(ChrMethInUSel,length(indDNASEOverlap))
      }
      
      ###############################################
      ##### Search for Meths and GenCode in TFBS #####
      ###############################################
      
      M<-nrow(CoordOverlapMeth)
      vecOverlapTFBS = rep(0,M)
      
      
      CoordOverlapTFBS<-cbind(rep(0,M),rep(0,M))
      
      out = .Fortran("matrixcompare3", as.matrix(CoordOverlapMeth), as.matrix(CoordTFBSC), as.integer(M), as.integer(nrow(CoordTFBSC)), as.numeric(vecOverlapTFBS), as.matrix(CoordOverlapTFBS), PACKAGE = "PoreMeth2")
      
      
      CoordOverlapTFBS<-out[[6]]
      vecOverlapTFBS = out[[5]]
      
      
      indTFBSOverlap<-which(vecOverlapTFBS!=0)
      MatrixOverlapTFBS<-cbind(rep("*",length(vecOverlapTFBS)),
                               rep("*",length(vecOverlapTFBS)),
                               rep("*",length(vecOverlapTFBS)),
                               rep("*",length(vecOverlapTFBS)),
                               rep("*",length(vecOverlapTFBS)),
                               CoordOverlapTFBS)
      if (length(indTFBSOverlap)!=0)
      {
        MatrixOverlapTFBS[indTFBSOverlap,1]<-NameTFBSC[vecOverlapTFBS[indTFBSOverlap]] 
        MatrixOverlapTFBS[indTFBSOverlap,2]<-rep(ChrMethInUSel,length(indTFBSOverlap))
        MatrixOverlapTFBS[indTFBSOverlap,3]<-CoordTFBSC[vecOverlapTFBS[indTFBSOverlap],1]
        MatrixOverlapTFBS[indTFBSOverlap,4]<-CoordTFBSC[vecOverlapTFBS[indTFBSOverlap],2]
        MatrixOverlapTFBS[indTFBSOverlap,5]<-rep(ChrMethInUSel,length(indTFBSOverlap))
      }
      
      
      ##################################
      ##### Creating Matrix GenCode #####
      ##################################
      
      CoordMeth2Test2Save<-matrix(rep((c(ChrMethInUSel,CoordMeth2Test)),each=length(indMethOverlap)),nrow=length(indMethOverlap))
      
      
      
      MatrixMeth2Test2Save<-matrix(rep(as.character(MatrixMeth2Test),each=length(indMethOverlap)),nrow=length(indMethOverlap))
      
      MatrixGenCodeOverlap<-cbind(rep(ChrMethInUSel,length(indMethOverlap)),
                                  CoordGenCodeC[indMethOverlap,1],
                                  CoordGenCodeC[indMethOverlap,2],
                                  FeatureGenCodeC[indMethOverlap],
                                  StrandGenCodeC[indMethOverlap],
                                  SymbolGenCodeC[indMethOverlap],
                                  TypeGenCodeC[indMethOverlap],
                                  rep(ChrMethInUSel,length(indMethOverlap)),
                                  CoordOverlap[indMethOverlap,],
                                  vecOverlapA[indMethOverlap],
                                  vecOverlapB[indMethOverlap])
      
      
      
      MatrixOut<-cbind(CoordMeth2Test2Save,
                       MatrixMeth2Test2Save,
                       MatrixGenCodeOverlap,
                       MatrixOverlapCGI,
                       MatrixOverlapEnhancer,
                       MatrixOverlapDNASE,
                       MatrixOverlapTFBS)
      
    }
    
    if (length(indMethOverlap)==1)
    {
      CoordOverlapMeth<-rbind(CoordOverlap[indMethOverlap,])
      
      
      #####################################################
      ##### Search for Meths and GenCode in CpG Islands #####
      #####################################################
      
      M<-nrow(CoordOverlapMeth)
      vecOverlapCGI = rep(0,M)
      
      
      CoordOverlapCGI<-cbind(rep(0,M),rep(0,M))
      
      out = .Fortran("matrixcompare3", as.matrix(CoordOverlapMeth), as.matrix(CoordCGIC), as.integer(M), as.integer(nrow(CoordCGIC)), as.numeric(vecOverlapCGI), as.matrix(CoordOverlapCGI), PACKAGE = "PoreMeth2")
      
      
      CoordOverlapCGI<-out[[6]]
      vecOverlapCGI = out[[5]]
      
      
      indCGIOverlap<-which(vecOverlapCGI!=0)
      MatrixOverlapCGI<-c(rep("*",length(vecOverlapCGI)),
                          rep("*",length(vecOverlapCGI)),
                          rep("*",length(vecOverlapCGI)),
                          rep("*",length(vecOverlapCGI)),
                          rep("*",length(vecOverlapCGI)),
                          CoordOverlapCGI)
      if (length(indCGIOverlap)!=0)
      {
        MatrixOverlapCGI[1]<-NameCGIC[vecOverlapCGI[indCGIOverlap]]
        MatrixOverlapCGI[2]<-rep(ChrMethInUSel,length(indCGIOverlap))
        MatrixOverlapCGI[3]<-CoordCGIC[vecOverlapCGI[indCGIOverlap],1]
        MatrixOverlapCGI[4]<-CoordCGIC[vecOverlapCGI[indCGIOverlap],2]
        MatrixOverlapCGI[5]<-rep(ChrMethInUSel,length(indCGIOverlap))
      }
      
      #####################################################
      #####  Search for Meths and GenCode in Enhancers  #####
      #####################################################
      
      M<-nrow(CoordOverlapMeth)
      vecOverlapEnhancer = rep(0,M)
      
      
      CoordOverlapEnhancer<-cbind(rep(0,M),rep(0,M))
      
      out = .Fortran("matrixcompare3", as.matrix(CoordOverlapMeth), as.matrix(CoordEnhancerC), as.integer(M), as.integer(nrow(CoordEnhancerC)), as.numeric(vecOverlapEnhancer), as.matrix(CoordOverlapEnhancer), PACKAGE = "PoreMeth2")
      
      
      CoordOverlapEnhancer<-out[[6]]
      vecOverlapEnhancer = out[[5]]
      
      
      indEnhancerOverlap<-which(vecOverlapEnhancer!=0)
      MatrixOverlapEnhancer<-c(rep("*",length(vecOverlapEnhancer)),
                               rep("*",length(vecOverlapEnhancer)),
                               rep("*",length(vecOverlapEnhancer)),
                               rep("*",length(vecOverlapEnhancer)),
                               rep("*",length(vecOverlapEnhancer)),
                               CoordOverlapEnhancer)
      if (length(indEnhancerOverlap)!=0)
      {
        MatrixOverlapEnhancer[1]<-NameEnhancerC[vecOverlapEnhancer[indEnhancerOverlap]]
        MatrixOverlapEnhancer[2]<-rep(ChrMethInUSel,length(indEnhancerOverlap))
        MatrixOverlapEnhancer[3]<-CoordEnhancerC[vecOverlapEnhancer[indEnhancerOverlap],1]
        MatrixOverlapEnhancer[4]<-CoordEnhancerC[vecOverlapEnhancer[indEnhancerOverlap],2]
        MatrixOverlapEnhancer[5]<-rep(ChrMethInUSel,length(indEnhancerOverlap))
      }
      
      ###############################################
      ##### Search for Meths and GenCode in DNASE #####
      ###############################################
      
      M<-nrow(CoordOverlapMeth)
      vecOverlapDNASE = rep(0,M)
      
      
      CoordOverlapDNASE<-cbind(rep(0,M),rep(0,M))
      
      out = .Fortran("matrixcompare3", as.matrix(CoordOverlapMeth), as.matrix(CoordDNASEC), as.integer(M), as.integer(nrow(CoordDNASEC)), as.numeric(vecOverlapDNASE), as.matrix(CoordOverlapDNASE), PACKAGE = "PoreMeth2")
      
      
      CoordOverlapDNASE<-out[[6]]
      vecOverlapDNASE = out[[5]]
      
      
      indDNASEOverlap<-which(vecOverlapDNASE!=0)
      MatrixOverlapDNASE<-c(rep("*",length(vecOverlapDNASE)),
                            rep("*",length(vecOverlapDNASE)),
                            rep("*",length(vecOverlapDNASE)),
                            rep("*",length(vecOverlapDNASE)),
                            rep("*",length(vecOverlapDNASE)),
                            CoordOverlapDNASE)
      if (length(indDNASEOverlap)!=0)
      {
        MatrixOverlapDNASE[1]<-NameDNASEC[vecOverlapDNASE[indDNASEOverlap]]
        MatrixOverlapDNASE[2]<-rep(ChrMethInUSel,length(indDNASEOverlap))
        MatrixOverlapDNASE[3]<-CoordDNASEC[vecOverlapDNASE[indDNASEOverlap],1]
        MatrixOverlapDNASE[4]<-CoordDNASEC[vecOverlapDNASE[indDNASEOverlap],2]
        MatrixOverlapDNASE[5]<-rep(ChrMethInUSel,length(indDNASEOverlap))
      }
      
      ###############################################
      ##### Search for Meths and GenCode in TFBS #####
      ###############################################
      
      M<-nrow(CoordOverlapMeth)
      vecOverlapTFBS = rep(0,M)
      
      
      CoordOverlapTFBS<-cbind(rep(0,M),rep(0,M))
      
      out = .Fortran("matrixcompare3", as.matrix(CoordOverlapMeth), as.matrix(CoordTFBSC), as.integer(M), as.integer(nrow(CoordTFBSC)), as.numeric(vecOverlapTFBS), as.matrix(CoordOverlapTFBS), PACKAGE = "PoreMeth2")
      
      
      CoordOverlapTFBS<-out[[6]]
      vecOverlapTFBS = out[[5]]
      
      
      indTFBSOverlap<-which(vecOverlapTFBS!=0)
      MatrixOverlapTFBS<-c(rep("*",length(vecOverlapTFBS))
                           ,rep("*",length(vecOverlapTFBS)),
                           rep("*",length(vecOverlapTFBS)),
                           rep("*",length(vecOverlapTFBS)),
                           rep("*",length(vecOverlapTFBS)),
                           CoordOverlapTFBS)
      if (length(indTFBSOverlap)!=0)
      {
        MatrixOverlapTFBS[1]<-NameTFBSC[vecOverlapTFBS[indTFBSOverlap]]
        MatrixOverlapTFBS[2]<-rep(ChrMethInUSel,length(indTFBSOverlap))
        MatrixOverlapTFBS[3]<-CoordTFBSC[vecOverlapTFBS[indTFBSOverlap],1]
        MatrixOverlapTFBS[4]<-CoordTFBSC[vecOverlapTFBS[indTFBSOverlap],2]
        MatrixOverlapTFBS[5]<-rep(ChrMethInUSel,length(indTFBSOverlap))
      }
      
      
      ##################################
      ##### Creating Matrix GenCode #####
      ##################################
      
      CoordMeth2Test2Save<-matrix(rep((c(ChrMethInUSel,CoordMeth2Test)),each=length(indMethOverlap)),nrow=length(indMethOverlap))
      
      
      
      MatrixMeth2Test2Save<-matrix(rep(as.character(MatrixMeth2Test),each=length(indMethOverlap)),nrow=length(indMethOverlap))
      
      MatrixGenCodeOverlap<-c(rep(ChrMethInUSel,length(indMethOverlap)),
                              CoordGenCodeC[indMethOverlap,1],
                              CoordGenCodeC[indMethOverlap,2],
                              FeatureGenCodeC[indMethOverlap],
                              StrandGenCodeC[indMethOverlap],
                              SymbolGenCodeC[indMethOverlap],
                              TypeGenCodeC[indMethOverlap],
                              rep(ChrMethInUSel,length(indMethOverlap)),
                              CoordOverlap[indMethOverlap,],
                              vecOverlapA[indMethOverlap],
                              vecOverlapB[indMethOverlap])
      
      
      
      MatrixOut<-rbind(c(CoordMeth2Test2Save,
                         MatrixMeth2Test2Save,
                         MatrixGenCodeOverlap,
                         MatrixOverlapCGI,
                         MatrixOverlapEnhancer,
                         MatrixOverlapDNASE,
                         MatrixOverlapTFBS))
      
    }
    
    if (length(indMethOverlap)==0)
    {
      
      CoordOverlapMeth<-rbind(CoordMeth2Test)
      
      
      #####################################################
      ##### Search for Meths and GenCode in CpG Islands #####
      #####################################################
      
      M<-nrow(CoordOverlapMeth)
      vecOverlapCGI = rep(0,M)
      
      
      CoordOverlapCGI<-cbind(rep(0,M),rep(0,M))
      
      out = .Fortran("matrixcompare3", as.matrix(CoordOverlapMeth), as.matrix(CoordCGIC), as.integer(M), as.integer(nrow(CoordCGIC)), as.numeric(vecOverlapCGI), as.matrix(CoordOverlapCGI), PACKAGE = "PoreMeth2")
      
      
      CoordOverlapCGI<-out[[6]]
      vecOverlapCGI = out[[5]]
      
      
      indCGIOverlap<-which(vecOverlapCGI!=0)
      MatrixOverlapCGI<-c(rep("*",length(vecOverlapCGI)),
                          rep("*",length(vecOverlapCGI)),
                          rep("*",length(vecOverlapCGI)),
                          rep("*",length(vecOverlapCGI)),
                          rep("*",length(vecOverlapCGI)),
                          CoordOverlapCGI)
      if (length(indCGIOverlap)!=0)
      {
        MatrixOverlapCGI[1]<-NameCGIC[vecOverlapCGI[indCGIOverlap]]
        MatrixOverlapCGI[2]<-rep(ChrMethInUSel,length(indCGIOverlap))
        MatrixOverlapCGI[3]<-CoordCGIC[vecOverlapCGI[indCGIOverlap],1]
        MatrixOverlapCGI[4]<-CoordCGIC[vecOverlapCGI[indCGIOverlap],2]
        MatrixOverlapCGI[5]<-rep(ChrMethInUSel,length(indCGIOverlap))
      }
      
      
      #####################################################
      ##### Search for Meths and GenCode in Enhancers #####
      #####################################################
      
      M<-nrow(CoordOverlapMeth)
      vecOverlapEnhancer = rep(0,M)
      
      
      CoordOverlapEnhancer<-cbind(rep(0,M),rep(0,M))
      
      out = .Fortran("matrixcompare3", as.matrix(CoordOverlapMeth), as.matrix(CoordEnhancerC), as.integer(M), as.integer(nrow(CoordEnhancerC)), as.numeric(vecOverlapEnhancer), as.matrix(CoordOverlapEnhancer), PACKAGE = "PoreMeth2")
      
      
      CoordOverlapEnhancer<-out[[6]]
      vecOverlapEnhancer = out[[5]]
      
      
      indEnhancerOverlap<-which(vecOverlapEnhancer!=0)
      MatrixOverlapEnhancer<-c(rep("*",length(vecOverlapEnhancer)),
                               rep("*",length(vecOverlapEnhancer)),
                               rep("*",length(vecOverlapEnhancer)),
                               rep("*",length(vecOverlapEnhancer)),
                               rep("*",length(vecOverlapEnhancer)),
                               CoordOverlapEnhancer)
      if (length(indEnhancerOverlap)!=0)
      {
        MatrixOverlapEnhancer[1]<-NameEnhancerC[vecOverlapEnhancer[indEnhancerOverlap]]
        MatrixOverlapEnhancer[2]<-rep(ChrMethInUSel,length(indEnhancerOverlap))
        MatrixOverlapEnhancer[3]<-CoordEnhancerC[vecOverlapEnhancer[indEnhancerOverlap],1]
        MatrixOverlapEnhancer[4]<-CoordEnhancerC[vecOverlapEnhancer[indEnhancerOverlap],2]
        MatrixOverlapEnhancer[5]<-rep(ChrMethInUSel,length(indEnhancerOverlap))
      }
      
      
      
      ###############################################
      ##### Search for Meths and GenCode in DNASE #####
      ###############################################
      
      M<-nrow(CoordOverlapMeth)
      vecOverlapDNASE = rep(0,M)
      
      
      CoordOverlapDNASE<-cbind(rep(0,M),rep(0,M))
      
      out = .Fortran("matrixcompare3", as.matrix(CoordOverlapMeth), as.matrix(CoordDNASEC), as.integer(M), as.integer(nrow(CoordDNASEC)), as.numeric(vecOverlapDNASE), as.matrix(CoordOverlapDNASE), PACKAGE = "PoreMeth2")
      
      
      CoordOverlapDNASE<-out[[6]]
      vecOverlapDNASE = out[[5]]
      
      
      indDNASEOverlap<-which(vecOverlapDNASE!=0)
      MatrixOverlapDNASE<-c(rep("*",length(vecOverlapDNASE)),
                            rep("*",length(vecOverlapDNASE)),
                            rep("*",length(vecOverlapDNASE)),
                            rep("*",length(vecOverlapDNASE)),
                            rep("*",length(vecOverlapDNASE)),
                            CoordOverlapDNASE)
      if (length(indDNASEOverlap)!=0)
      {
        MatrixOverlapDNASE[1]<-NameDNASEC[vecOverlapDNASE[indDNASEOverlap]]
        MatrixOverlapDNASE[2]<-rep(ChrMethInUSel,length(indDNASEOverlap))
        MatrixOverlapDNASE[3]<-CoordDNASEC[vecOverlapDNASE[indDNASEOverlap],1]
        MatrixOverlapDNASE[4]<-CoordDNASEC[vecOverlapDNASE[indDNASEOverlap],2]
        MatrixOverlapDNASE[5]<-rep(ChrMethInUSel,length(indDNASEOverlap))
      }
      
      ###############################################
      ##### Search for Meths and GenCode in TFBS #####
      ###############################################
      
      M<-nrow(CoordOverlapMeth)
      vecOverlapTFBS = rep(0,M)
      
      
      CoordOverlapTFBS<-cbind(rep(0,M),rep(0,M))
      
      out = .Fortran("matrixcompare3", as.matrix(CoordOverlapMeth), as.matrix(CoordTFBSC), as.integer(M), as.integer(nrow(CoordTFBSC)), as.numeric(vecOverlapTFBS), as.matrix(CoordOverlapTFBS), PACKAGE = "PoreMeth2")
      
      
      CoordOverlapTFBS<-out[[6]]
      vecOverlapTFBS = out[[5]]
      
      
      indTFBSOverlap<-which(vecOverlapTFBS!=0)
      MatrixOverlapTFBS<-c(rep("*",length(vecOverlapTFBS)),
                           rep("*",length(vecOverlapTFBS)),
                           rep("*",length(vecOverlapTFBS)),
                           rep("*",length(vecOverlapTFBS)),
                           rep("*",length(vecOverlapTFBS)),
                           CoordOverlapTFBS)
      if (length(indTFBSOverlap)!=0)
      {
        MatrixOverlapTFBS[1]<-NameTFBSC[vecOverlapTFBS[indTFBSOverlap]]
        MatrixOverlapTFBS[2]<-rep(ChrMethInUSel,length(indTFBSOverlap))
        MatrixOverlapTFBS[3]<-CoordTFBSC[vecOverlapTFBS[indTFBSOverlap],1]
        MatrixOverlapTFBS[4]<-CoordTFBSC[vecOverlapTFBS[indTFBSOverlap],2]
        MatrixOverlapTFBS[5]<-rep(ChrMethInUSel,length(indTFBSOverlap))
      }
      
      
      ##################################
      ##### Creating Matrix GenCode #####
      ##################################
      
      CoordMeth2Test2Save<-matrix(rep((c(ChrMethInUSel,CoordMeth2Test)),each=1),nrow=1)
      
      
      
      MatrixMeth2Test2Save<-matrix(rep(as.character(MatrixMeth2Test),each=1),nrow=1)
      
      MatrixGenCodeOverlap<-c("*",
                              "*",
                              "*",
                              "InterGenic",
                              "*",
                              "*",
                              "*",
                              "*",
                              "*",
                              "*",
                              "*",
                              "*")
      
      
      
      MatrixOut<-rbind(c(CoordMeth2Test2Save,
                         MatrixMeth2Test2Save,
                         MatrixGenCodeOverlap,
                         MatrixOverlapCGI,
                         MatrixOverlapEnhancer,
                         MatrixOverlapDNASE,
                         MatrixOverlapTFBS))
      
    }
  }
  MatrixOut
}
