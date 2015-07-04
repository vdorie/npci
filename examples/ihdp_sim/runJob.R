start   <- as.integer(Sys.getenv("START"))
end     <- as.integer(Sys.getenv("END"))
method  <- Sys.getenv("METHOD")

require(npci)

resultsFileName <- paste0(method, "_", start, "_", end, ".RData")
resultsFile <- file.path("data", resultsFileName)

if (file.exists(resultsFileName)) q("no")

numReps <- end - start + 1

results <- matrix(NA, numReps, 5)

colnames(results) <- 
  c("bias", "cov", "cil", "wrong", "tau.est")

precision <- rep(NA_real_, numReps)

source("data.R")

x.dgp <- cbind(1, as.matrix(x))
w <- matrix(c(0, rep(0.5, ncol(x.dgp) - 1)), nrow(x.dgp), ncol(x.dgp), byrow = TRUE)

for (i in seq_len(numReps)){ 
  iter <- i + start - 1
  
  ## places mu.0, mu.1, y.0, y.1, and y into calling env
  generateDataForIterInCurrentEnvironment(iter, x.dgp, z, w)
    
  meanEffects  <- mu.1[z == 0] - mu.0[z == 0]
  results[i, "tau.est"] <- mean(y.1[z == 0] -  y.0[z == 0])
  
  if (method == "bart") {
    ## bart section
    
    ## te - treatment effect, or calculated difference
    ## each sample is a treatment effect, by person and by replication
    treatmentEffectSamples <- ci.estimate(y, x, z, method = "bart", estimand = "atc")
    atcSamples <- apply(treatmentEffectSamples, 2, mean)          ## average of people first
    results[i, "bias"] <- 4 - mean(atcSamples)
    #ci <- quantile(atcSamples, c(0.025, 0.975))
    ci <- mean(atcSamples) + sd(atcSamples) * qnorm(c(0.025, 0.975))
    results[i, "cov"] <- if (ci[1] < 4 && ci[2] > 4) 1 else 0
    results[i, "cil"] <- ci[2] - ci[1]
    results[i, "wrong"] <- if (ci[1] < 0 && ci[2] > 0) 1 else 0
  
    treatmentEffects <- apply(treatmentEffectSamples, 1, mean)    ## average over samples
    precision[i] <- sqrt(mean((treatmentEffects - meanEffects)^2))
    
  } else if (method == "naive1" || method == "naive2") {
    atc <- ci.estimate(y, x, z, method = method, estimand = "atc")
    est <- mean(atc$te)
    se  <- atc$se
    results[i, "bias"] <- 4 - est
    ci <- est + se * qnorm(c(0.025, 0.975))
    results[i, "cov"] <- if (ci[1] < 4 && ci[2] > 4) 1 else 0
    results[i, "cil"] <- ci[2] - ci[1]
    results[i, "wrong"] <- if (ci[1] < 0 && ci[2] > 0) 1 else 0
    precision[i] <- sqrt(mean((atc$te - meanEffects)^2))
  } else {
    stop("method '", method, "' unrecognized")
  }
}

save(results, precision, file = resultsFile)
