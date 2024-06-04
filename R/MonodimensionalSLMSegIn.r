#' Truncated Gaussian Monodimensional SLM Segmentation
#' @NoRd

MonodimensionalSLMSegIn <- function(DataMatrix, muk, mi, smu, sepsilon, Pos, omega, eta, stepeta) {

    PathSrc <- paste0(.libPaths()[1], "/PoreMeth2/libs/")
    dyn.load(paste0(PathSrc, "PoreMeth2.so"))

    ####  Calcolo il vettore delle covariate dipendenti dalla distanza ###
    CovPos <- diff(Pos)
    CovPosNorm <- CovPos / stepeta
    etavec <- eta + ((1 - eta) * exp(log(eta) / CovPosNorm))


    ### Calcolo i parametri del modello SLM
    NCov <- length(etavec)
    K0 <- ncol(muk)
    etav <- log(rep(1, K0) * (1 / K0))
    T <- ncol(DataMatrix)
    NExp <- nrow(DataMatrix)

    TruncCoef <- matrix(0, nrow = NExp, ncol = K0)
    for (jj in 1:NExp)
    {
        for (kk in 1:length(muk))
        {
            TruncCoef[jj, kk] <- sepsilon[jj] * (pnorm(1, mean = muk[jj, kk], sd = sepsilon[jj]) - pnorm(-1, mean = muk[jj, kk], sd = sepsilon[jj]))
        }
    }

    ####  Calcolo le Matrici di Emissione e Transizione ####
    P <- matrix(data = 0, nrow = K0, ncol = (K0 * NCov))
    G <- matrix(data = 0, nrow = K0, ncol = K0)
    emission <- matrix(data = 0, nrow = K0, ncol = T)
    out <- .Fortran("transemisislm", as.vector(muk), as.vector(mi), as.double(etavec), as.integer(NCov), as.matrix(DataMatrix), as.integer(K0), as.integer(NExp), as.vector(smu), as.vector(sepsilon), as.integer(T), as.matrix(G), as.matrix(P), as.matrix(emission), as.matrix(TruncCoef), PACKAGE = "PoreMeth2")

    P <- out[[12]]
    emission <- out[[13]]




    ####### Algoritmo di Viterbi #########
    psi <- matrix(data = 0, nrow = K0, ncol = T)
    path <- c(as.integer(rep(0, T)))
    out2 <- .Fortran("bioviterbii", as.vector(etav), as.matrix(P), as.matrix(emission), as.integer(T), as.integer(K0), as.vector(path), as.matrix(psi), PACKAGE = "PoreMeth2")
    s <- out2[[6]]


    sortResult <- SortState(s)
    TotalPredBreak <- sortResult[[3]]
    TotalPredBreak
}
