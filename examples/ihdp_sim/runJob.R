start      <- as.integer(Sys.getenv("START"))
end        <- as.integer(Sys.getenv("END"))
method     <- Sys.getenv("METHOD")
overlap    <- if ((overlap <- Sys.getenv("OVERLAP")) == "") TRUE else as.logical(overlap)
covariates <- if ((covariates <- Sys.getenv("COVARIATES")) == "") "full" else covariates
setting    <- if ((setting <- Sys.getenv("SETTING")) == "") "A" else setting
p.score    <- if ((p.score <- Sys.getenv("PROPENSITY")) == "") FALSE else as.logical(p.score)
verbose    <- if ((verbose <- Sys.getenv("VERBOSE")) == "") TRUE else as.logical(verbose)

require(npci)
source("results.R")

prefix <- getPrefix(method, overlap, covariates, setting, p.score)

resultsFileName <- paste0(prefix, "_", start, "_", end, ".RData")
resultsFile <- file.path("data", resultsFileName)

if (file.exists(resultsFileName)) q("no")

numReps <- end - start + 1L

results <- matrix(NA, numReps, 6L)

colnames(results) <- 
  c("bias", "cov", "cil", "wrong", "tau.est", "precision")
opt <- list()
L <- list()
x.0 <- list()
opt.0 <- list()
L.0 <- list()
x.0.0 <- list()
opt.1 <- list()
L.1 <- list()
x.0.1 <- list()

source("data.R")
loadDataInCurrentEnvironment(covariates, p.score)

w <- 0.5

estimand <- if (identical(overlap, TRUE)) "att" else "atc"

for (i in seq_len(numReps)) { 
  iter <- i + start - 1L
  
  if (verbose) cat("  running iter ", iter, "\n", sep = "")
  
  ## places mu.0, mu.1, y.0, y.1, and y into calling env, as well as x.r and z.r
  generateDataForIterInCurrentEnvironment(iter, x, z, w, overlap, covariates, setting)
  ## x.r and z.r are adjusted by simulation if necessary
  
  if (p.score == TRUE)
    x.r <- cbind(x.r, ps.z) 
  
  meanEffects <- if (identical(overlap, TRUE)) mu.1[z.r == 1] - mu.0[z.r == 1] else mu.1[z.r == 0] - mu.0[z.r == 0]
  results[i, "tau.est"] <- if (identical(overlap, TRUE)) mean(y.1[z.r == 1] -  y.0[z.r == 1]) else mean(y.1[z.r == 0] -  y.0[z.r == 0])
  
  if (method == "bart") {
    ## bart section
    
    ## each sample is a treatment effect, by person and by replication
    treatmentEffectSamples <- ci.estimate(y, x.r, z.r, method = "bart", estimand = estimand)
    ## average over people first, produces samples of ATC or ATT
    estimandSamples <- apply(treatmentEffectSamples, 2L, mean)
    results[i,"bias"] <- 4 - mean(estimandSamples)
    #ci <- quantile(atcSamples, c(0.025, 0.975))
    ci <- mean(estimandSamples) + sd(estimandSamples) * qnorm(c(0.025, 0.975))
    results[i,"cov"] <- if (ci[1L] < 4 && ci[2L] > 4) 1 else 0
    results[i,"cil"] <- ci[2L] - ci[1L]
    results[i,"wrong"] <- if (ci[1L] < 0 && ci[2L] > 0) 1 else 0
  
    ## average over samples, produces per-person TE estimates
    treatmentEffects <- apply(treatmentEffectSamples, 1L, mean)
    results[i,"precision"] <- sqrt(mean((treatmentEffects - meanEffects)^2))
    
  } else if (method == "bartipw") {
    treatmentEffectSamples <- ci.estimate(y, x.r, z.r, method = "bart", estimand = estimand, prob.z = ps.z)
    atcSamples <- apply(treatmentEffectSamples, 2L, mean)
    results[i,"bias"] <- 4 - mean(atcSamples)
    #ci <- quantile(atcSamples, c(0.025, 0.975))
    ci <- mean(atcSamples) + sd(atcSamples) * qnorm(c(0.025, 0.975))
    results[i,"cov"] <- if (ci[1L] < 4 && ci[2L] > 4) 1 else 0
    results[i,"cil"] <- ci[2L] - ci[1L]
    results[i,"wrong"] <- if (ci[1L] < 0 && ci[2L] > 0) 1 else 0
  
    treatmentEffects <- apply(treatmentEffectSamples, 1, mean)
    results[i,"precision"] <- sqrt(mean((treatmentEffects - meanEffects)^2))
  } else if (method %in% c("naive1", "naive2")) {
    est <- ci.estimate(y, x.r, z.r, method = method, estimand = estimand)
    te <- mean(est$te)
    se  <- est$se
    results[i,"bias"] <- 4 - te
    ci <- te + se * qnorm(c(0.025, 0.975))
    results[i,"cov"] <- if (ci[1L] < 4 && ci[2L] > 4) 1 else 0
    results[i,"cil"] <- ci[2L] - ci[1L]
    results[i,"wrong"] <- if (ci[1L] < 0 && ci[2L] > 0) 1 else 0
    results[i,"precision"] <- sqrt(mean((est$te - meanEffects)^2))
  } else if (identical(method, "grouped")) {
    est <- ci.estimate(y, x.r, z.r, method = method, estimand = estimand, verbose = verbose)
    te <- mean(est$te)
    se  <- est$se
    results[i,"bias"] <- 4 - te
    ci <- te + se * qnorm(c(0.025, 0.975))
    results[i,"cov"] <- if (ci[1L] < 4 && ci[2L] > 4) 1 else 0
    results[i,"cil"] <- ci[2L] - ci[1L]
    results[i,"wrong"] <- if (ci[1L] < 0 && ci[2L] > 0) 1 else 0
    results[i,"precision"] <- sqrt(mean((est$te - meanEffects)^2))
    
    opt[[i]] <- est$fit@env$opt
    devEnv <- environment(est$fit@deviance)
    L[[i]]   <- environment(devEnv$transformPars)$L
    x.0[[i]] <- environment(devEnv$transformPars)$x.0
  } else if (identical(method, "independent")) {
    est <- ci.estimate(y, x.r, z.r, method = method, estimand = estimand, verbose = verbose)
    te <- mean(est$te)
    se  <- est$se
    results[i,"bias"] <- 4 - te
    ci <- te + se * qnorm(c(0.025, 0.975))
    results[i,"cov"] <- if (ci[1L] < 4 && ci[2L] > 4) 1 else 0
    results[i,"cil"] <- ci[2L] - ci[1L]
    results[i,"wrong"] <- if (ci[1L] < 0 && ci[2L] > 0) 1 else 0
    results[i,"precision"] <- sqrt(mean((est$te - meanEffects)^2))
    
    opt.0[[i]] <- est$fit@env$opt.0
    opt.1[[i]] <- est$fit@env$opt.1
    devEnv <- environment(est$fit@deviance.0)
    L.0[[i]]   <- environment(devEnv$transformPars.0)$L
    L.1[[i]]   <- environment(devEnv$transformPars.1)$L
    x.0.0[[i]] <- environment(devEnv$transformPars.0)$x.0
    x.0.1[[i]] <- environment(devEnv$transformPars.1)$x.0
  } else {
    stop("method '", method, "' unrecognized")
  }
}

if (identical(method, "grouped")) {
  save(results, opt, L, x.0, file = resultsFile)
} else if (identical(method, "independent")) {
  save(results, opt.0, opt.1, L.0, L.1, x.0.0, x.0.1, file = resultsFile)
} else {
  save(results, file = resultsFile)
}
