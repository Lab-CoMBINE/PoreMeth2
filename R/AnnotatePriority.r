#' Annotation priorization based on genomic feature
#'
#' Returns a annotation priorizatized by genomic feature
#'
#' @param TableDMRAnnoIn Is the output from `PoreMethAnnotate`, data.frame.
#' @param NumProc Number of threads. Default = 1, numeric.
#' 
#' @return data.frame with Prioritized Annotation.
#'
#' @export


AnnotatePriority <- function(TableDMRAnnoIn, NumProc=1)
{


    FieldNames<-colnames(TableDMRAnnoIn)
    indFinal<-which(FieldNames=="PValue")


    FieldNamesOut<-c(FieldNames[1:indFinal],"Genic.Element","symbol.GenCode","type.GenCode","name.CGI","name.Enhancer","name.DNASE","name.TFBS")
    ChrStartEndString<-paste(TableDMRAnnoIn$chr,TableDMRAnnoIn$start,TableDMRAnnoIn$end,sep="")
    ChrStartEndStringU<-unique(ChrStartEndString)

    cl <- parallel::makePSOCKcluster(NumProc)
    parallel::setDefaultCluster(cl)
    parallel::clusterExport(NULL, c('RegionAnnotatePriority'))



    RegionAnnotatePriorityOut<-parLapply(cl,c(1:length(ChrStartEndStringU)),RegionAnnotatePriority,
                                        TableDMRAnnoIn,
                                        ChrStartEndString,
                                        ChrStartEndStringU,
                                        indFinal)

    parallel::stopCluster(cl)

    return(do.call("rbind", RegionAnnotatePriorityOut))

}