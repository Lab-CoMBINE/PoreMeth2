#' Function to filter out small segments
#' @NoRd

FilterSeg <- function(TotalPredBreak, FW) {
    controllength <- diff(TotalPredBreak)

    indF <- which(controllength <= FW)
    if (length(indF) != 0) {
        if (indF[1] == 1) {
            indF[1] <- 2
            indF <- unique(indF)
            TotalPredBreak1 <- TotalPredBreak[-(indF)]
        }
        if (indF[1] != 1) {
            TotalPredBreak1 <- TotalPredBreak[-(indF)]
        }
    }
    if (length(indF) == 0) {
        TotalPredBreak1 <- TotalPredBreak
    }
    TotalPredBreak1
}
