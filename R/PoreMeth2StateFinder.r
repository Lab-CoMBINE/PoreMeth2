#' Function for HMM segmentation
#'
#' description description description description description description
#'
#' @param SignalIn
#' @param PosIn numeric,
#' @param muk numeric,
#' @param sepsilon numeric,
#' @param NormDist numeric,
#' @param PTStart numeric,
#'
#' @return what what what what what what what what what
#'
#' @export

PoreMeth2StateFinder <- function(SignalIn, PosIn, muk, sepsilon, NormDist, PTStart) {

    PathSrc <- paste0(.libPaths()[1], "/PoreMeth2/libs/")
    dyn.load(paste0(PathSrc, "PoreMeth2.so"))


    KS <- length(muk)
    CovPos <- diff(PosIn)
    CovDist <- CovPos / NormDist
    CovDist1 <- log(1 - exp(-CovDist))
    T <- length(SignalIn)


    NCov <- length(CovDist)
    TruncCoef <- rep(0, length(muk))
    for (kk in 1:length(muk))
    {
        TruncCoef[kk] <- sepsilon[kk] * (pnorm(1, mean = muk[kk], sd = sepsilon[kk]) - pnorm(0, mean = muk[kk], sd = sepsilon[kk]))
    }


    PT <- log(rep(PTStart, KS))
    P <- matrix(data = 0, nrow = KS, ncol = (KS * NCov))
    emission <- matrix(data = 0, nrow = KS, ncol = T)


    #### Calculates Transition and Emission Probabilities ##
    out <- .Fortran("transemisihmm", as.vector(muk), as.integer(NCov), as.vector(SignalIn), as.integer(KS), as.vector(CovDist1), as.vector(sepsilon), as.integer(T), as.matrix(PT), as.matrix(P), as.matrix(emission), as.vector(TruncCoef), PACKAGE = "PoreMeth2")

    P <- out[[9]]
    emission <- out[[10]]

    ##### Viterbi Algorithm ####
    etav <- log(rep(1, KS) * (1 / KS))
    psi <- matrix(data = 0, nrow = KS, ncol = T)
    path <- c(as.integer(rep(0, T)))
    out2 <- .Fortran("bioviterbii", as.vector(etav), as.matrix(P), as.matrix(emission), as.integer(T), as.integer(KS), as.vector(path), as.matrix(psi), PACKAGE = "PoreMeth2")
    s <- out2[[6]]
    s
}
