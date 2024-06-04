#' Function for plotting DMRs, Entropy levels and Gene Annotations
#' 
#' This function plot the results from PoreMeth2. It presents various parameters for plot customization. 
#'
#' @param Input Either a string character of coordinates to plot (chr:start-end) or gene symbol
#' @param AnnotatedRes Annotated results from PoreMeth2. Data.frame
#' @param PoreMeth2DMRResults Poremeth results
#' @param Meth1 Output from perl parser. Data.frame, chr, pos, entropy, number of reads for entropy estimation, beta level, number of reads for beta levels estimation
#' @param Meth2 Output from perl parser. Data.frame, chr, pos, entropy, number of reads for entropy estimation, beta level, number of reads for beta levels estimation
#' @param CGI Shows CpG Islands. TRUE or FALSE. Default TRUE
#' @param DNAse Shows DNAse sites. TRUE or FALSE. Default TRUE
#' @param Enhancers Shows Enhancers sites. TRUE or FALSE. Default TRUE
#' @param TFBS Shows TFBS sites. TRUE or FALSE. Default TRUE
#' @param CpGPositions Shows CpGPositions sites. TRUE or FALSE.  Default FALSE
#' @param AddRange Zoom out from input coordinates/gene. Can be "auto" or an integer. Default "auto"
#' @param CoefLM Coeficient for polynomial regression for beta and entropy levels. Numeric value. Default 0.3
#' @param Assembly Reference annotation to use. Can be "hg19" or "hg38". Default "hg19"
#' @param Palette Palette colors to use. Can be "Pal1" to "Pal12", otherwhise a custom color vector can be used (length > 5). Default = "Pal12".
#' @param TxtGeneColor Color for gene symbols. May improve readibility. Default "black"
#' @param DMRBeta Shows delta beta levels. TRUE or FALSE. Default TRUE
#' @param DMREntropy Shows delta entropy levels. TRUE or FALSE. Default TRUE
#' @param DMRBetaPoly Shows polynomial regression line of beta levels for test and control. TRUE or FALSE. Default TRUE
#' @param DMRPoints Shows the delta betas for each CpG on which DMRs are calculated. TRUE or FALSE. Default FALSE
#' @param DMRLines Shows raw lines linking the points on which DMRs are calculated. TRUE or FALSE. Default FALSE
#' @param FeatureLines As DMRLines but in a different plot. TRUE or FALSE, Default = FALSE 
#' @param FeaturePoints As DMRPoints but in a different plot. TRUE or FALSE, Default = FALSE
#' @param DMREntropyPoly Shows polynomial regression line for entropy levels for test and control. TRUE or FALSE. Default TRUE
#' @param EntropyPoints Shows the delta entropy for each CpG on which DMRs are calculated. TRUE or FALSE. Default FALSE
#' @param Grid Shows horizontal lines for better value estimation. TRUE or FALSE, Default FALSE
#' @param CoordLines Shows vertical lines for DMRs. It helps to compare plots. Default = TRUE 
#' @param ShowCpGNames Shows GpG IDs. TRUE or FALSE, Default FALSE
#' @param DMRLinesW Line width for DMREntropyPoly. Default 1
#' @param DMRLinesT Line type for DMREntropyPoly. Default 1
#' @param PointsCex Cex for test Points. Default 0.5
#' @param FeatureLinesW Line width for FeatureLines. Default 1
#' @param FeatureLinesT Line type for FeatureLines. Default 2
#' @param GeneNameOffset Offset on X axis of the gene names. Use in case the gene name ovelaps gene models. Default = 0. 
#' @return Plot(s)
#'
#' @export

PoreMeth2Plot <- function(
    Input = NA, AnnotatedRes, PoreMeth2DMRResults, Meth1, Meth2, 
    CGI = TRUE, DNAse = TRUE, Enhancers = TRUE, TFBS = TRUE, CpGPositions = FALSE,
    AddRange = 'auto', CoefLM = 0.3, Assembly = "hg19", Palette = "Pal12", TxtGeneColor = "black",
    DMRBeta = TRUE, DMREntropy = FALSE, DMRBetaPoly = FALSE, DMRPoints = TRUE, DMRLines = FALSE, FeatureLines = FALSE, FeaturePoints = FALSE, DMREntropyPoly = FALSE, EntropyPoints = TRUE,
    Grid = TRUE, CoordLines = TRUE, ShowCpGNames = FALSE,
    DMRLinesW = 1, DMRLinesT = 1, PointsCex = 0.5,
    FeatureLinesW = 1, FeatureLinesT = 2, FeatureLinesC1 = NULL, FeatureLinesC2 = NULL, FeaturePointsC = NULL, GeneNameOffset = 0
){
    # Generating Palettes
    if(class(Palette) == "character"){
        Pals <- list(
            Pal1 = c("#153B2C","#1C5C42","#97A397","#D4D7C4","#E3BD6A","#DCC695"),
            Pal2 = c("#264653","#2A9D8F","#E9C46A","#F4A261","#EE8959","#E76F51"),
            Pal3 = c("#A36361","#D3A29D","#E8B298","#EECC8C","#BDD1C5","#9EABA2"),
            Pal4 = c("#003049","#6B2C39","#D62828","#F77F00","#FCBF49","#EAE2B7"),
            Pal5 = c("#361C0E","#570211","#7E3110","#004540","#032C4D","#360825"),
            Pal6 = c("#264653","#2A9D8F","#E9C46A","#F4A261","#E76F51","#EC8C74"),
            Pal7 = rev(c("#153B2C","#1C5C42","#97A397","#D4D7C4","#E3BD6A","#DCC695")),
            Pal8 = rev(c("#264653","#2A9D8F","#E9C46A","#F4A261","#EE8959","#E76F51")),
            Pal9 = rev(c("#A36361","#D3A29D","#E8B298","#EECC8C","#BDD1C5","#9EABA2")),
            Pal10 = rev(c("#003049","#6B2C39","#D62828","#F77F00","#FCBF49","#EAE2B7")),
            Pal11 = rev(c("#361C0E","#570211","#7E3110","#004540","#032C4D","#360825")),
            Pal12 = rev(c("#264653","#2A9D8F","#E9C46A","#F4A261","#E76F51","#EC8C74"))
        )
        Pal <- Pals[[Palette]]
    } else {
        Pal <- Palette
    }

    # Loading Annotation Data
    PathDBIn <- paste0(.libPaths()[1], "/PoreMeth2/data/")

    if(all(DMRBetaPoly, DMRLines)){
        stop("DMRBetaPoly and DMRLines can not be both true")
    }

    if(!exists("FileGenCodeIn")){
        FileGenCodeIn <- file.path(PathDBIn, paste("GencodeTable_", Assembly, ".rds", sep = ""))
        TableGenCodeIn <<- data.frame(readRDS(FileGenCodeIn))
    }

    if(!exists("FileCGIIn") & CGI){
        FileCGIIn <- file.path(PathDBIn, paste("CGIAnno_", Assembly, ".rds", sep = ""))
        TableCGIIn <<- data.frame(readRDS(FileCGIIn))
    }

    if(!exists("FileDNAseIn") & DNAse){
        FileDNAseIn <- file.path(PathDBIn, paste("DNASEAnno_", Assembly, ".rds", sep = ""))
        TableDNAseIn <<- data.frame(readRDS(FileDNAseIn))
    }

    if(!exists("FileEnhancerIn") & Enhancers){
        FileEnhancerIn <- file.path(PathDBIn, paste("EnhancerAnno_", Assembly, ".rds", sep = ""))
        TableEnhancerIn <<- data.frame(readRDS(FileEnhancerIn))
    }

    if(!exists("FileTFBSIn") & TFBS){
        FileTFBSIn <- file.path(PathDBIn, paste("TFBSAnno_", Assembly, ".rds", sep = ""))
        TableTFBSIn <<-  data.frame(readRDS(FileTFBSIn))
    }

    if(class(AnnotatedRes)[1]!="data.frame"){
        AnnotatedRes <<- data.frame(AnnotatedRes)
    }

    # Check the input: is it a coordinate or a gene name?
    if(grepl("\\d{5,}|[-:]", Input)){
        Coord <- unlist(strsplit(Input, "[:\\-]"))
        Coord[2] <- as.numeric(Coord[2])
        Coord[3] <- as.numeric(Coord[3])
        if(!grepl("chr", Coord[1])) Coord[1] <- paste0("chr", Coord[1])
        IndDMRTarget <- suppressWarnings(which(AnnotatedRes[,1] == Coord[1] & as.numeric(AnnotatedRes[,2]) >= as.numeric(Coord[2]) & as.numeric(AnnotatedRes[,3]) <= as.numeric(Coord[3])))
        TargetFeatureCoord <- AnnotatedRes[IndDMRTarget, 12:14]
        if(nrow(TargetFeatureCoord)<1) stop("No genomic features to plot, check the input coordinates")
    } else {
        CoordDF <- subset(TableGenCodeIn, TableGenCodeIn$Gene.Symbol %in% Input)
        Coord <- c(CoordDF[1,1], range(c(CoordDF[,2], CoordDF[,3])))
        if(!grepl("chr", Coord[1])) Coord[1] <- paste0("chr", Coord[1])
        IndDMRTarget <- suppressWarnings(which(AnnotatedRes[,1] == Coord[1] & as.numeric(AnnotatedRes[,13]) >= as.numeric(Coord[2]) & as.numeric(AnnotatedRes[,14]) <= as.numeric(Coord[3])))
        TargetFeatureCoord <- data.frame(rbind(Coord))
        if(nrow(TargetFeatureCoord)<1) stop("No genomic features to plot, check the input coordinates")
    }

    # Extracting Info From Annotation
    TargetDMRCoord <- AnnotatedRes[IndDMRTarget, c(1:3)]
    TargetGene <- AnnotatedRes[IndDMRTarget, "symbol.GenCode"]
    TargetGeneUnique <- unique(TargetGene)
    TargetStrand <- AnnotatedRes[IndDMRTarget, "strand.GenCode"]
    TargetGeneAnnotation <- subset(TableGenCodeIn, TableGenCodeIn[, 4] %in% TargetGeneUnique)
    TargetGeneAnnotation <- split(TargetGeneAnnotation, TargetGeneAnnotation$"Gene.Symbol")[TargetGeneUnique]
    TargetBeta <- AnnotatedRes[IndDMRTarget, c(4, 6, 7)]
    TargetEntropy <- AnnotatedRes[IndDMRTarget, c(5, 8, 9)]
    TargetChrMeth1 <- subset(Meth1, Meth1[, 1] %in% gsub("chr", "", Coord[1]))
    TargetChrMeth2 <- subset(Meth2, Meth2[, 1] %in% gsub("chr", "", Coord[1]))

    # Other Info from Full Annotations
    tryCatch({
        TargetCGICoord <- AnnotatedRes[IndDMRTarget, c("chr.CGI", "start.CGI", "end.CGI", "name.CGI")]
        TargetEnhCoord <- AnnotatedRes[IndDMRTarget, c("chr.Enhancer", "start.Enhancer", "end.Enhancer", "name.Enhancer")]
        TargetDNAseCoord <- AnnotatedRes[IndDMRTarget, c("chr.DNASE", "start.DNASE", "end.DNASE")]
        TargetTFBSCoord <- AnnotatedRes[IndDMRTarget, c("chr.TFBS", "start.TFBS", "end.TFBS")]
    },
        error = function(e) {
    })

    # Plot Ranges and Rows
    if(AddRange == 'auto') AddRange <- suppressWarnings(round((max(as.numeric(TargetFeatureCoord[, 3]), na.rm = T) - min(as.numeric(TargetFeatureCoord[, 2]), na.rm = T)) * 0.1,0))
    PlotRange <- suppressWarnings(as.numeric(c(min(as.numeric(TargetFeatureCoord[, 2]), na.rm = T), max(as.numeric(TargetFeatureCoord[, 3]), na.rm = T))) + (rep(AddRange,2) * c(-1,1)))
    PlotRows <- length(TargetGeneUnique)

    MatrixLayout <- matrix(
        c(
            rep(1, any(DMRBeta, DMRBetaPoly, DMRPoints, DMRLines) * 2),
            rep(any(DMRBeta, DMRBetaPoly, DMRPoints, DMRLines) + 1, any(FeatureLines, FeaturePoints) * 2),
            rep(any(DMRBeta, DMRBetaPoly, DMRPoints, DMRLines) +  any(FeatureLines, FeaturePoints) + 1, any(DMREntropy, DMREntropyPoly) * 2),
            rep(any(DMRBeta, DMRBetaPoly, DMRPoints, DMRLines) + any(FeatureLines, FeaturePoints) + any(DMREntropy, DMREntropyPoly) + 1, 1 * 1),
            rep(any(DMRBeta, DMRBetaPoly, DMRPoints, DMRLines) + any(FeatureLines, FeaturePoints) + any(DMREntropy, DMREntropyPoly) + 2, any(CpGPositions, CGI, DNAse, Enhancers, TFBS) * 1)

        ),
        ncol = 1
    )

    # Window Methylation / Entropy levels
    if(DMRBeta){
        CurPoreMeth2DMRResults <- subset(PoreMeth2DMRResults, PoreMeth2DMRResults[, 1] == gsub("chr", "", Coord[1]))
        CurPoreMeth2DMRResults <- subset(CurPoreMeth2DMRResults, as.numeric(CurPoreMeth2DMRResults$p) <= 0.05)
        WindowBLevels <- rbind(
            subset(CurPoreMeth2DMRResults, as.numeric(CurPoreMeth2DMRResults[, 2]) >= PlotRange[1] & as.numeric(CurPoreMeth2DMRResults[, 3]) <= PlotRange[2]),
            subset(CurPoreMeth2DMRResults, as.numeric(CurPoreMeth2DMRResults[, 3]) >= PlotRange[1] & as.numeric(CurPoreMeth2DMRResults[, 3]) <= PlotRange[2]),
            subset(CurPoreMeth2DMRResults, as.numeric(CurPoreMeth2DMRResults[, 2]) >= PlotRange[1] & as.numeric(CurPoreMeth2DMRResults[, 2]) <= PlotRange[2])
        )
        WindowBLevels <- unique(WindowBLevels)
    }

    if(any(DMRPoints, DMRLines, DMRBetaPoly, CpGPositions, FeaturePoints, FeatureLines)){
        # Retrieving Info about common mCs within the DMR
        WindowCpGLevels1 <- subset(TargetChrMeth1, TargetChrMeth1[, 2] >= PlotRange[1])
        WindowCpGLevels1 <- subset(WindowCpGLevels1, WindowCpGLevels1[, 2] <= PlotRange[2])

        WindowCpGLevels2 <- subset(TargetChrMeth2, TargetChrMeth2[, 2] >= PlotRange[1])
        WindowCpGLevels2 <- subset(WindowCpGLevels2, WindowCpGLevels2[, 2] <= PlotRange[2])

        WindowCpGLevels1 <- subset(WindowCpGLevels1, WindowCpGLevels1[, 2] %in% WindowCpGLevels2[, 2])
        WindowCpGLevels2 <- subset(WindowCpGLevels2, WindowCpGLevels2[, 2] %in% WindowCpGLevels1[, 2])
    }

    # Extra annotations
    if(CGI) {
        WindowCGI <- subset(TableCGIIn, TableCGIIn[, 1] == Coord[1])
        WindowCGI <- subset(WindowCGI, (WindowCGI[, 2] >= PlotRange[1] & WindowCGI[, 2] <= PlotRange[2]) | (WindowCGI[, 3] >= PlotRange[1] & WindowCGI[, 3] <= PlotRange[2]))
    }

    if (DNAse) {
        WindowDNAse <- subset(TableDNAseIn, TableDNAseIn[, 1] == Coord[1])
        WindowDNAse <- subset(WindowDNAse, (WindowDNAse[, 2] >= PlotRange[1] & WindowDNAse[, 2] <= PlotRange[2]) | (WindowDNAse[, 3] >= PlotRange[1] & WindowDNAse[, 3] <= PlotRange[2]))
    }

    if (Enhancers){
        WindowEnhancers <- subset(TableEnhancerIn, TableEnhancerIn[, 1] == Coord[1])
        WindowEnhancers <- subset(WindowEnhancers, (WindowEnhancers[, 2] >= PlotRange[1] & WindowEnhancers[, 2] <= PlotRange[2]) | (WindowEnhancers[, 3] >= PlotRange[1] & WindowEnhancers[, 3] <= PlotRange[2]))
    }

    if (TFBS){
        WindowTFBS <- subset(TableTFBSIn, gsub("chr", "", TableTFBSIn[, 1]) == Coord[1])
        WindowTFBS <- subset(WindowTFBS, (WindowTFBS[, 2] >= PlotRange[1] & WindowTFBS[, 2] <= PlotRange[2]) | (WindowTFBS[, 3] >= PlotRange[1] & WindowTFBS[, 3] <= PlotRange[2]))
    }

    # Generating Plot
    layout(MatrixLayout)

    # Colors
    ColInd <- sum(
        DMRBeta,
        (any(DMRBetaPoly, DMRPoints, DMRLines) * 2),
        DMREntropy,
        (any(DMREntropyPoly, FeatureLines, FeaturePoints) * 2)
    )

    # Beta plot
    if(any(DMRBeta, DMRBetaPoly, DMRPoints, DMRLines)) {
        par(mar = c(0.7, 4.1, 5.1, 1.1))
        plot(x = 1, y = 1, col = "white", xlim = c(PlotRange[1], PlotRange[2]), ylim = c(-1, 1), main = Input, ylab = expression(paste(Delta, Beta)), xlab = "", bty = "n", xaxt = "n")
        if (Grid) {
            abline(h = seq(-1, 1, 0.1), col = "grey85", lwd = 0.1)
        }

        if (DMRPoints) {
            DMRPointsC1 <- Pal[ColInd]
            DMRPointsC2 <- Pal[ColInd - 1]
            ys <- WindowCpGLevels1[, 5] - WindowCpGLevels2[, 5]
            points(x = WindowCpGLevels2[, 2], y = ys, xlim = c(PlotRange[1], PlotRange[2]), ylim = c(-.1, 1), col = adjustcolor(ifelse(ys < 0, DMRPointsC1, DMRPointsC2), alpha.f = .15), pch = 21, cex = PointsCex)
        }

        if (DMRBetaPoly) {
            DMRPolyC2 <- Pal[ColInd]
            DMRPolyC1 <- Pal[ColInd - 1]

            CoefLMFract <- round(min(length(unique(WindowCpGLevels1[, 5])) - 1, length(unique(WindowCpGLevels2[, 5])) - 1) * CoefLM, 0)

            while(TRUE) {
                # try 
                tryCatch({
                    Model1 <- lm(WindowCpGLevels1[, 5] ~ poly(WindowCpGLevels1[, 2], CoefLMFract))
                    Model2 <- lm(WindowCpGLevels2[, 5] ~ poly(WindowCpGLevels2[, 2], CoefLMFract))
                    break
                },
                # if error, decrease
                    error = function(e){
                        CoefLMFract <<- CoefLMFract - 1
                    }
                )
            }

            Model1YFit <- predict(Model1)
            lines(WindowCpGLevels1[, 2], Model1YFit, col = DMRPolyC1, lwd = DMRLinesW, lty = DMRLinesT)

            Model2YFit <- predict(Model2)
            lines(WindowCpGLevels2[, 2], Model2YFit, col = DMRPolyC2, lwd = DMRLinesW, lty = DMRLinesT)
        }

        if (DMRLines) {
            DMRLinesC1 <- Pal[ColInd]
            DMRLinesC2 <- Pal[ColInd - 1]
            lines(x = WindowCpGLevels1[, 2], y = WindowCpGLevels1[, 5], xlim = c(PlotRange[1], PlotRange[2]), ylim = c(-.1, 1), col = DMRLinesC1, lwd = DMRLinesW, lty = DMRLinesT)
            lines(x = WindowCpGLevels2[, 2], y = WindowCpGLevels2[, 5], xlim = c(PlotRange[1], PlotRange[2]), ylim = c(-.1, 1), col = DMRLinesC2, lwd = DMRLinesW, lty = DMRLinesT)
        }

        if (DMRBeta) {
            CurDMRCols1 <- Pal[ColInd]
            CurDMRCols2 <- Pal[ColInd - 1]
            rect(xright = WindowBLevels[, 3], xleft = WindowBLevels[, 2], ybottom = 0, ytop = WindowBLevels[, 4], col = adjustcolor(ifelse(WindowBLevels[, 4] > 0 , CurDMRCols1, CurDMRCols2), alpha.f = 0.7), border = ifelse(WindowBLevels[, 4] > 0 , CurDMRCols1, CurDMRCols2), lwd = 1)
            abline(h = 0)
            ColInd <- ColInd - 2
        }

        if(CoordLines) abline(v = c(WindowBLevels[, 3], WindowBLevels[, 2]), col = 'grey85', lty = 2, lwd = 0.2)
    }

    if (any(DMREntropy, DMREntropyPoly, EntropyPoints)) {
        par(mar = c(0.7, 4.1, 1.5, 1.1))
        plot(x = 1, y = 1, col = "white", xlim = c(PlotRange[1], PlotRange[2]), ylim = c(-1, 1), ylab = "Entropy Levels", xlab = "", xaxt = "n", bty = "n", main = Input, col.main = ifelse(any(DMRBeta, DMRBetaPoly, DMRPoints, DMRLines), "white", "black"))
        if (Grid) {
            abline(h = seq(-1, 1, 0.1), col = "grey95", lwd = 0.1)
        }

        if (EntropyPoints) {
            EntropyPointsC1 <- Pal[ColInd]
            EntropyPointsC2 <- Pal[ColInd - 1]
            ys <- WindowCpGLevels1[, 3] - WindowCpGLevels2[, 3]
            points(x = WindowCpGLevels2[, 2], y = ys, xlim = c(PlotRange[1], PlotRange[2]), ylim = c(-.1, 1), col = adjustcolor(ifelse(ys < 0, EntropyPointsC1, EntropyPointsC2), alpha.f = .2), pch = 21, cex = PointsCex)
        }

        if(DMREntropyPoly){
            DMREntropyPolyC2 <- Pal[ColInd]
            DMREntropyPolyC1 <- Pal[ColInd - 1]
            CoefLMFract <- round(min(length(unique(WindowCpGLevels1[, 3])) - 1, length(unique(WindowCpGLevels2[, 3])) - 1) * CoefLM, 0)
            while(TRUE){
                tryCatch({
                    Model1 <- lm(WindowCpGLevels1[, 3] ~ poly(WindowCpGLevels1[, 2], CoefLMFract))
                    Model2 <- lm(WindowCpGLevels2[, 3] ~ poly(WindowCpGLevels2[, 2], CoefLMFract))
                    break
                    },
                        error = function(e){
                            CoefLMFract <<- CoefLMFract - 1
                        }
                    )
                }

            Model1YFit <- predict(Model1)
            lines(WindowCpGLevels1[, 2], Model1YFit, col = DMREntropyPolyC1, lwd = DMRLinesW, lty = DMRLinesT)

            Model2YFit <- predict(Model2)
            lines(WindowCpGLevels2[, 2], Model2YFit, col = DMREntropyPolyC2, lwd = DMRLinesW, lty = DMRLinesT)
        }

        if(DMREntropy){
            EntropyCols1 <- Pal[ColInd]
            EntropyCols2 <- Pal[ColInd - 1]
            rect(xright = WindowBLevels[, 3], xleft = WindowBLevels[, 2], ybottom = 0, ytop = WindowBLevels[, 5], col = adjustcolor(ifelse(WindowBLevels[, 5] > 0 , EntropyCols1, EntropyCols2), alpha.f = 0.7), border = ifelse(WindowBLevels[, 5] > 0 , EntropyCols1, EntropyCols2), lwd = 1)
            abline(h = 0)
            ColInd <- ColInd - 2
        }

        if(CoordLines) abline(v = c(WindowBLevels[, 3], WindowBLevels[, 2]), col = 'grey85', lty = 2, lwd = 0.2)

    }

    if(any(FeaturePoints, FeatureLines)){
        par(mar = c(0.7, 4.1, 0.5, 1.1))
        plot(x = 1, y = 1, col = "white", xlim = c(PlotRange[1], PlotRange[2]), ylim = c(0, 1), ylab = expression(paste(Beta, "Levels")), xlab = "", xaxt = "n", bty = "n", main = Input, col.main = ifelse(any(DMRBeta, DMRBetaPoly, DMRPoints, DMRLines), "white", "black"))
        if (Grid) {
            abline(h = seq(-1, 1, 0.1), col = "grey95", lwd = 0.1)
        }
        if (FeaturePoints) {
            FeaturePointsC1 <- Pal[ColInd]
            FeaturePointsC2 <- Pal[ColInd - 1]
            ys <- WindowCpGLevels1[, 5] - WindowCpGLevels2[, 5]
            points(x = WindowCpGLevels2[, 2], y = ys, xlim = c(PlotRange[1], PlotRange[2]), ylim = c(-.1, 1), col = adjustcolor(ifelse(ys < 0, FeaturePointsC1, FeaturePointsC2), alpha.f = .15), pch = 21, cex = PointsCex)
        }
        if (FeatureLines) {
            FeatureLinesC1 <- Pal[ColInd]
            FeatureLinesC2 <- Pal[ColInd - 1]
            ys <- WindowCpGLevels1[, 5] - WindowCpGLevels2[, 5]
            lines(x = WindowCpGLevels2[, 2], y = ys, xlim = c(PlotRange[1], PlotRange[2]), ylim = c(-.1, 1), col = adjustcolor(ifelse(ys < 0, FeatureLinesC1, FeatureLinesC2), alpha.f = .15),  lty = FeatureLinesT, lwd = FeatureLinesW)
        }
    }

    # Plot annotations
    ColInd <- sum(
        DMRBeta,
        (any(DMRBetaPoly, DMRPoints, DMRLines) * 2),
        DMREntropy,
        (any(DMREntropyPoly, FeaturePoints, FeatureLines) * 2)
    )

    N <- 0
    NExtra <- seq(1, sum(CpGPositions, CGI, DNAse, Enhancers, TFBS))
    par(mar = c(0.7, 4.1, 0.9, 1.1))
    plot(x = 1, y = 1, col = "white", xlim = c(PlotRange[1], PlotRange[2]), ylim = c(-.2 * PlotRows, 0), axes = FALSE, ylab = "Symbols", xlab = "")
    if(length(NExtra) == 0){
        axis(labels = TRUE, side = 1, seq(PlotRange[1], PlotRange[2], length.out = 6), cex = 0.5, srt = 45)
    }

    for(N1 in 1:length(TargetGeneUnique)){
        # Feature Plotting Parameters
        FeatCoordLeft <- TargetGeneAnnotation[[N1]][, 2]
        FeatCoordRight <- TargetGeneAnnotation[[N1]][, 3]
        FeatColors <- ifelse(grepl("F?L?.+Exon", TargetGeneAnnotation[[N1]][, 6]), "dodgerblue3", ifelse(grepl("^Exon", TargetGeneAnnotation[[N1]][, 6]), "dodgerblue1", ifelse(grepl("Promoter", TargetGeneAnnotation[[N1]][, 6]), "Yellow", "grey50")))
        FeatShapes <- ifelse(grepl("Prom", TargetGeneAnnotation[[N1]][, 6]), 0.5, ifelse(grepl("Exon", TargetGeneAnnotation[[N1]][, 6]), 1, 0.1))

        # Genes
        CurrentGene <- TargetGeneUnique[N1]
        N = N + 1
        IndCurGene <- which(TargetGene == CurrentGene)
        rect(xleft = FeatCoordLeft, xright = FeatCoordRight, ybottom = (-0.08 * N) + (FeatShapes / 30), ytop = -(0.08 * N) - (FeatShapes / 30), col = FeatColors)
    }
    if(CoordLines) abline(v = c(FeatCoordLeft, FeatCoordRight), col = 'grey85', lty = 2, lwd = 0.2)

    # Render the text over the gene features 
    N = 0
    for(N1 in 1:length(TargetGeneUnique)){
        FeatCoordLeft <- TargetGeneAnnotation[[N1]][, 2]
        CurrentGene <- TargetGeneUnique[N1]
        N <- N + 1
        IndCurGene <- which(TargetGene == CurrentGene)
        CurGeneText <- paste0(CurrentGene, ifelse(unique(TargetStrand[IndCurGene]) == "+", " -->", " <--"))
        text(CurGeneText, x = max(min(FeatCoordLeft), PlotRange[1]), y = (-0.09 * N) - .05, adj = GeneNameOffset, cex = 0.9, col = TxtGeneColor)
    }

    # Plotting extra info
    IndNExtra <- length(NExtra)

    if(length(NExtra)>0){
        par(mar = c(3.1,4.1,0.1,2.1))
        plot(x = 1, y = 1, col = "white", xlim = c(PlotRange[1], PlotRange[2]), ylim = c(-.2, NExtra[IndNExtra]), axes = FALSE, ylab = "Genomic Info", xlab = "")
        axis(labels = TRUE, side = 1, seq(PlotRange[1], PlotRange[2], length.out = 6), cex = 0.5, srt = 45)
    }
    if (CpGPositions) {
        points(x = WindowCpGLevels1[, 2], y = rep(NExtra[IndNExtra], length(WindowCpGLevels1[, 2])), pch = 15, col = Pal[IndNExtra], cex = 0.3)
        text("CpG Positions", x = PlotRange[1], adj = 0, y = (NExtra[IndNExtra] - 0.4))
        IndNExtra <- IndNExtra - 1
        ColInd <- ColInd - 1
    }
    if (CGI & nrow(WindowCGI) > 0) {
        segments(x0 = PlotRange[1], x1 = PlotRange[2], y0 = NExtra[IndNExtra], y1 = NExtra[IndNExtra], col = Pal[IndNExtra], lwd = 0.2)
        rect(xleft = WindowCGI[, 2], xright = WindowCGI[, 3], ybottom = rep(NExtra[IndNExtra] - 0.1, nrow(WindowCGI)), ytop = rep(NExtra[IndNExtra] + 0.1, nrow(WindowCGI)), col = Pal[IndNExtra], xpd = T)
        if(ShowCpGNames) text(WindowCGI[, 4], x = sapply(WindowCGI[, 2], function(x){max(x,PlotRange[1])}), y = NExtra[IndNExtra] - 0.4, xpd = T)
        text("CGIs", x = PlotRange[1], adj = 0, y = NExtra[IndNExtra] - 0.4, xpd = T)
        IndNExtra <- IndNExtra - 1
        ColInd <- ColInd - 1
    }
    if(DNAse) {
        segments(x0 = PlotRange[1], x1 = PlotRange[2], y0 = NExtra[IndNExtra], y1 = NExtra[IndNExtra], col = Pal[IndNExtra], lwd = 0.2)
        rect(xleft = WindowDNAse[, 2], xright = WindowDNAse[, 3], ybottom = rep(NExtra[IndNExtra] - 0.1, nrow(WindowDNAse)), ytop = rep(NExtra[IndNExtra] + 0.1, nrow(WindowDNAse)), col = Pal[IndNExtra])
        text("DNAse Sites", x = PlotRange[1], adj = 0, y = NExtra[IndNExtra] - 0.4, xpd = T)
        IndNExtra <- IndNExtra - 1
        ColInd <- ColInd - 1
    }
    if (Enhancers) {
        segments(x0 = PlotRange[1], x1 = PlotRange[2], y0 = NExtra[IndNExtra], y1 = NExtra[IndNExtra], col = Pal[IndNExtra], lwd = 0.2)
        rect(xleft = WindowEnhancers[, 2], xright = WindowEnhancers[, 3], ybottom = rep(NExtra[IndNExtra] - 0.1, nrow(WindowEnhancers)), ytop = rep(NExtra[IndNExtra] + 0.1, nrow(WindowEnhancers)), col = Pal[IndNExtra])
        text("Enhancers", x = PlotRange[1], adj = 0, y = NExtra[IndNExtra] - 0.4, xpd = T)
        IndNExtra <- IndNExtra - 1
        ColInd <- ColInd - 1
    }
    if (TFBS) {
        segments(x0 = PlotRange[1], x1 = PlotRange[2], y0 = NExtra[IndNExtra], y1 = NExtra[IndNExtra], col = Pal[IndNExtra], lwd = 0.2)
        rect(xleft = WindowTFBS[, 2], xright = WindowTFBS[, 3], ybottom = rep(NExtra[IndNExtra] - 0.1, nrow(WindowTFBS)), ytop = rep(NExtra[IndNExtra] + 0.1, nrow(WindowTFBS)), col = Pal[IndNExtra])
        text("TFBS", x = PlotRange[1], adj = 0, y = NExtra[IndNExtra] - 0.4, xpd = T)
        IndNExtra <- IndNExtra - 1
        ColInd <- ColInd - 1
    }
}
