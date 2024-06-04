#' Function to calculate the mean value of each segment
#' @NoRd


SegResults <- function(DataSeq, TotalPredBreak) {
    TotalPred <- c()
    NExp <- nrow(DataSeq)
    for (j in 1:NExp)
    {
        s <- rep(0, ncol(DataSeq))
        for (i in 1:(length(TotalPredBreak) - 1))
        {
            s[(TotalPredBreak[i] + 1):TotalPredBreak[i + 1]] <- median(DataSeq[j, (TotalPredBreak[i] + 1):TotalPredBreak[i + 1]])
        }
        TotalPred <- rbind(TotalPred, s)
    }

    Result <- TotalPred
    Result
}
