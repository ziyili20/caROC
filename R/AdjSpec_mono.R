AdjSpec_mono <- function(diseaseData,
                         controlData,
                         userFormula,
                         fixSS,
                         whichSE = "numerical",
                         nbootstrap = 100,
                         CI_alpha = 0.95,
                         logit_CI = TRUE,
                         verbose = TRUE) {

    rec_pEst = rec_pEst_mono = rep(0, length(fixSS))
    rec_pSE <- rep(0, length(fixSS))
    for(i in 1:length(fixSS)) {
        EstSE <- GetEst_SE(diseaseData = diseaseData,
                           controlData = controlData,
                           userFormula = userFormula,
                           fixSS = fixSS[i],
                           whichSE = whichSE,
                           nbootstrap = 100,
                           verbose = verbose)
        rec_pEst[i] <- EstSE$Estimates[length(EstSE$Estimates)]
        rec_pSE[i] <- EstSE$SE[length(EstSE$SE)]
    }
    rqfit <- rq(userFormula, tau = -1, data = diseaseData)
    nvar <- dim(EstSE$all_Z_C)[2]
    all_Z <- rbind(EstSE$all_Z_D, EstSE$all_Z_C)
    estparam <- rbind(rqfit$sol[c(4:(4+nvar),1), ])
    monoout <- mono(estparam, as.matrix(cbind(1, all_Z)), taus = 1-fixSS)
    ctrl_threshold_q_mono <- as.matrix(cbind(1, EstSE$all_Z_C)) %*% monoout$bt
    if(length(fixSS) == 1) {
        rec_pEst_mono <- sum(EstSE$M0 <= ctrl_threshold_q_mono)/length(EstSE$M0)
    } else {
        rec_pEst_mono <- colSums(EstSE$M0 <= ctrl_threshold_q_mono)/length(EstSE$M0)
    }

    AllCI <- GetAllCI(pEst = rec_pEst_mono,
                      pSE = rec_pSE,
                      alpha = CI_alpha,
                      logit_CI = logit_CI)
    colnames(AllCI) <- fixSS

    return(list(Estimate = rec_pEst_mono,
                SE = rec_pSE,
                ConfidenceInterval = AllCI))
}
