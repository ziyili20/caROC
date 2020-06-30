GetAllCI <- function(pEst, pSE,
                     alpha = 0.95,
                     logit_CI = TRUE) {

    sai <- qnorm(1 - (1-alpha)/2)
    allCI <- matrix(0, 2, length(pEst))
    if (logit_CI) {
        allCI[1, ] <- revlogit(logit(pEst) - sai * pSE/pEst/(1-pEst))
        allCI[2, ] <- revlogit(logit(pEst) + sai * pSE/pEst/(1-pEst))
        rownames(allCI) <- c("LogitCI_Lower", "LogitCI_Upper")
    } else {
        allCI[1, ] <- pEst - sai * pSE
        allCI[2, ] <- pEst + sai * pSE
        rownames(allCI) <- c("CI_Lower", "CI_Upper")
    }

    return(allCI)
}
