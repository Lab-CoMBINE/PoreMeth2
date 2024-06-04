#' Segmentation State Sorting
#' @NoRd

SortState <- function(s) {
    l <- 1
    seg <- c()
    brek <- c()
    t <- 1
    for (k in 1:(length(s) - 1))
    {
        if (s[k] != s[k + 1]) {
            brek[t] <- k
            t <- t + 1
            if (length(which(seg == s[k])) == 0) {
                seg[l] <- s[k]
                l <- l + 1
            }
        }
    }
    brek <- c(0, brek, length(s))
    if (length(which(seg == s[length(s)])) == 0) {
        seg <- c(seg, s[length(s)])
    }

    s0 <- c()
    for (k in 1:length(seg))
    {
        s0[which(s == seg[k])] <- k
    }

    SortResult <- list()
    SortResult[[1]] <- s0
    SortResult[[2]] <- seg
    SortResult[[3]] <- brek
    SortResult
}
