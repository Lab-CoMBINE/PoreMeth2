#' Annotation Priority 
#' @NoRd

RegionAnnotatePriority <- function(Input,TableDMRAnnoIn,
                                   ChrStartEndString,
                                   ChrStartEndStringU,
                                   indFinal)
{
  MatrixGeneOut<-c()
  indSel<-which(ChrStartEndString==ChrStartEndStringU[Input])
  TableDMRAnnoInSel<-TableDMRAnnoIn[indSel,]
  indGene<-which(TableDMRAnnoInSel$type.GenCode!="*")
  if (length(indGene)!=0)
  {
    TableDMRAnnoInSelGenes<-rbind(TableDMRAnnoInSel[indGene,])
    GeneSel<-TableDMRAnnoInSelGenes$symbol.GenCode
    GeneSelU<-unique(GeneSel)
    for (kk in 1:length(GeneSelU))
    {
      indG<-which(GeneSel==GeneSelU[kk])
      TableDMRAnnoInSelG<-rbind(TableDMRAnnoInSelGenes[indG,])
      FeatureG<-TableDMRAnnoInSelG$feature.GenCode
      CGIG<-TableDMRAnnoInSelG$name.CGI
      TFBSG<-TableDMRAnnoInSelG$name.TFBS
      DNASEG<-TableDMRAnnoInSelG$name.DNASE
      EnhancerG<-TableDMRAnnoInSelG$name.Enhancer
     
      indReg5<-which(FeatureG=="Promoter" | FeatureG=="FirstExon")
      if (length(indReg5)!=0)
      {
       
        TableDMRAnnoInSelGReg5<-rbind(TableDMRAnnoInSelG[indReg5,])
       
        CGIGReg5<-CGIG[indReg5]
        indCGI<-which(CGIGReg5!="*")
        if (length(indCGI)!=0)
        {
          MatrixGeneOut<-rbind(MatrixGeneOut,as.character(c(TableDMRAnnoInSelGReg5[indCGI[1],1:indFinal],
                                                            "5'Reg",
                                                            TableDMRAnnoInSelGReg5$symbol.GenCode[indCGI[1]],
                                                            TableDMRAnnoInSelGReg5$type.GenCode[indCGI[1]],
                                                            CGIGReg5[indCGI[1]],
                                                            TableDMRAnnoInSelGReg5$name.Enhancer[indCGI[1]],
                                                            TableDMRAnnoInSelGReg5$name.DNASE[indCGI[1]],
                                                            TableDMRAnnoInSelGReg5$name.TFBS[indCGI[1]])))
         
        }
        if (length(indCGI)==0)
        {
         
         
          MatrixGeneOut<-rbind(MatrixGeneOut,as.character(c(TableDMRAnnoInSelGReg5[1,1:indFinal],
                                                            "5'Reg",
                                                            TableDMRAnnoInSelGReg5$symbol.GenCode[1],
                                                            TableDMRAnnoInSelGReg5$type.GenCode[1],
                                                            CGIGReg5[1],
                                                            TableDMRAnnoInSelGReg5$name.Enhancer[1],
                                                            TableDMRAnnoInSelGReg5$name.DNASE[1],
                                                            TableDMRAnnoInSelGReg5$name.TFBS[1])))
         
        }
      }
      if (length(indReg5)==0)
      {
        indLE<-which(FeatureG=="LastExon")
        if (length(indLE)!=0)
        {
          TableDMRAnnoInSelGLE<-rbind(TableDMRAnnoInSelG[indLE,])
         
          CGIGLE<-CGIG[indLE]
          indCGI<-which(CGIGLE!="*")
          if (length(indCGI)!=0)
          {
            MatrixGeneOut<-rbind(MatrixGeneOut,as.character(c(TableDMRAnnoInSelGLE[indCGI[1],1:indFinal],
                                                              "3'UTR",
                                                              TableDMRAnnoInSelGLE$symbol.GenCode[indCGI[1]],
                                                              TableDMRAnnoInSelGLE$type.GenCode[indCGI[1]],
                                                              CGIGLE[indCGI[1]],
                                                              TableDMRAnnoInSelGLE$name.Enhancer[indCGI[1]],
                                                              TableDMRAnnoInSelGLE$name.DNASE[indCGI[1]],
                                                              TableDMRAnnoInSelGLE$name.TFBS[indCGI[1]])))
           
          }
          if (length(indCGI)==0)
          {
            MatrixGeneOut<-rbind(MatrixGeneOut,as.character(c(TableDMRAnnoInSelGLE[1,1:indFinal],
                                                              "3'UTR",
                                                              TableDMRAnnoInSelGLE$symbol.GenCode[1],
                                                              TableDMRAnnoInSelGLE$type.GenCode[1],
                                                              CGIGLE[1],
                                                              TableDMRAnnoInSelGLE$name.Enhancer[1],
                                                              TableDMRAnnoInSelGLE$name.DNASE[1],
                                                              TableDMRAnnoInSelGLE$name.TFBS[1])))
          }
        }
        if (length(indLE)==0)
        {
          indGB<-c(grep('^Exon',FeatureG),grep('^Intron',FeatureG))
          if (length(indGB)!=0)
          {
            TableDMRAnnoInSelGGB<-rbind(TableDMRAnnoInSelG[indGB,])
           
            CGIGGB<-CGIG[indGB]
            indCGI<-which(CGIGGB!="*")
            if (length(indCGI)!=0)
            {
              MatrixGeneOut<-rbind(MatrixGeneOut,as.character(c(TableDMRAnnoInSelGGB[indCGI[1],1:indFinal],
                                                                "GB",
                                                                TableDMRAnnoInSelGGB$symbol.GenCode[indCGI[1]],
                                                                TableDMRAnnoInSelGGB$type.GenCode[indCGI[1]],
                                                                CGIGGB[indCGI[1]],
                                                                TableDMRAnnoInSelGGB$name.Enhancer[indCGI[1]],
                                                                TableDMRAnnoInSelGGB$name.DNASE[indCGI[1]],
                                                                TableDMRAnnoInSelGGB$name.TFBS[indCGI[1]])))
             
            }
            if (length(indCGI)==0)
            {
              MatrixGeneOut<-rbind(MatrixGeneOut,as.character(c(TableDMRAnnoInSelGGB[1,1:indFinal],
                                                                "GB",
                                                                TableDMRAnnoInSelGGB$symbol.GenCode[1],
                                                                TableDMRAnnoInSelGGB$type.GenCode[1],
                                                                CGIGGB[1],
                                                                TableDMRAnnoInSelGGB$name.Enhancer[1],
                                                                TableDMRAnnoInSelGGB$name.DNASE[1],
                                                                TableDMRAnnoInSelGGB$name.TFBS[1])))
            }
          }
        }
      }
    }
   
   
  }
  MatrixGeneOut
}