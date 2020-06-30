rootFunc <- function(whichtau,
                     est_Z,
                     alldata,
                     alldata_C,
                     M0,
                     rowN = 1) {
    outres <- mono(est_Z, alldata, initau = 0.5, taus = whichtau)
    resid0 <- alldata_C[rowN,] %*% outres$bt - M0[rowN]
    return(resid0)
}
