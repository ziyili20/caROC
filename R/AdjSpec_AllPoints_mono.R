AdjSpec_AllPoints_mono <- function(diseaseData,
                                   controlData,
                                   userFormula) {

    outform <- ParseuserFormula(userFormula = userFormula)
    M0 <- eval(parse(text = paste("controlData$",outform[[1]])))
    M1 <- eval(parse(text = paste("diseaseData$",outform[[1]])))
    if (length(outform[[2]]) == 1) {
        expr2 <- paste0("controlData$", outform[[2]])
        expr3 <- paste0("diseaseData$", outform[[2]])
    } else if (length(outform[[2]]) > 1) {
        expr2 <- paste0("controlData[, ", outform[2], "]")
        expr3 <- paste0("diseaseData[, ", outform[2], "]")
    }
    all_Z_C <- eval(parse(text = expr2))
    all_Z_D <- eval(parse(text = expr3))

    rqfit <- NULL
    expr1 <- paste0("rqfit <- rq(", userFormula, ", tau = -1, data = diseaseData)")
    eval(parse(text = expr1))

    all_tau <- rqfit$sol[1, ]
    sens_vec <- 1 - all_tau
    est_b <- rqfit$sol[4:(nrow(rqfit$sol)),]

    ### run mono
    est_Z <- rbind(est_b, all_tau)
    alldata <- as.matrix(cbind(1, rbind(all_Z_C, all_Z_D)))
    alldata_C <- as.matrix(cbind(1, all_Z_C))

    rec_breaks <- rep(0, nrow(controlData))
    for(i in 1:nrow(controlData)) {
        miii <- rootFunc(0,
                         est_Z,
                         alldata,
                         alldata_C,
                         M0,
                         rowN = i)
        maaa <- rootFunc(1,
                         est_Z,
                         alldata,
                         alldata_C,
                         M0,
                         rowN = i)
        if(miii > 0 & maaa > 0) {
            rec_breaks[i] <- 0
        } else if (miii < 0 & maaa < 0) {
            rec_breaks[i] <- 1
        } else {
            uuout <- uniroot(function(x){rootFunc(x, est_Z,
                                                  alldata,
                                                  alldata_C,
                                                  M0,rowN = i)}, c(0,1))
            rec_breaks[i] <- uuout$root
        }
    }

    all_tau2 <- sort(unique(rec_breaks))
    all_spec_mono <- rep(0, length(all_tau2))
    for(j in 1:length(all_tau2)) {
        all_spec_mono[j] <- sum(rec_breaks < all_tau2[j])/length(rec_breaks)
    }
    sens_vec <- 1 - all_tau2

    return(list(sensitivity = sens_vec,
                specificity = all_spec_mono,
                mono_adj = "mono"))

}
