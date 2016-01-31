loadDataInCurrentEnvironment <- function(covariates = "select", p.score = "none") {
  callingEnv <- parent.frame(1L)
  
  if (is.character(covariates) && covariates != "full") {
    dataFile <- file.path("data", "ihdp.RData")
    if (!file.exists(dataFile)) stop("ihdp data file not found at path: ", dataFile)

    load(dataFile)
    
    ihdp <- subset(ihdp, treat != 1 | momwhite != 0)
    
    covariateNames <- c("bw", "b.head", "preterm", "birth.o", "nnhealth", "momage",
                        "sex", "twin", "b.marr", "mom.lths", "mom.hs", "mom.scoll",
                        "cig", "first", "booze", "drugs", "work.dur", "prenatal",
                        "ark", "ein", "har", "mia", "pen", "tex", "was")

    x <- ihdp[,covariateNames]
    
    trans <- npci:::getTransformations(x)

    x <- npci:::transform(x, trans$standardize)
    z <- ihdp$treat
    
    if (p.score != "none") {
      propFile <- file.path("data", paste0("prop_", p.score, ".RData"))
      if (!file.exists(propFile)) stop("propensity score file not found at path: ", propFile)
      
      load(propFile)
    }
  } else {
    dataFile <- file.path("data", "ihdpFull.RData")
    if (!file.exists(dataFile)) stop("ihdp full data file not found at path: ", dataFile)
    
    load(dataFile)
    
    ihdpFull <- subset(ihdpFull, treat != 1 | momwhite != 0)
    dropCols <- "treat"
    x <- ihdpFull[,colnames(ihdpFull)[match(colnames(ihdpFull), dropCols, nomatch = 0) == 0]]
    
    trans <- npci:::getTransformations(x)

    x <- npci:::transform(x, trans$standardize)
    z <- ihdpFull$treat
    
    if (p.score != "none") {
      propFile <- file.path("data", paste0("propFull_", p.score, ".RData"))
      if (!file.exists(propFile)) stop("full propensity score file not found at path: ", propFile)
      
      load(propFile)
    }
  }
  
  callingEnv$x <- x
  callingEnv$z <- z
  if (p.score != "none")
    callingEnv$ps.z <- ps.z
}

getPropensityScoreInCurrentEnvironment <- function(x, z, p.score = "none")
{
  getQuadraticTerms <- function(x)
  {
    terms <- character()
    isBinary <- sapply(seq_len(ncol(x)), function(j) length(unique(x[,j])) == 2)
    for (i in seq_len(ncol(x))) {
      if (!isBinary[i]) terms[length(terms) + 1L] <- paste0("I(", colnames(x)[i], "^2)")
      
      if (i < ncol(x)) for (j in seq(i + 1L, ncol(x)))
        terms[length(terms) + 1L] <- paste0(colnames(x)[i], ":", colnames(x)[j])
    }
    terms
  }

  callingEnv <- parent.frame(1L)
  
  if (p.score == "none" || p.score == "true") return(invisible(NULL))
  
  if (p.score == "logistic") {
    df <- as.data.frame(x)
    df$z <- z
    m <- glm(z ~ ., data = df, family = binomial())
    
    callingEnv$ps.z <- fitted(m)
  } else if (p.score == "bayesglm") {
    df <- as.data.frame(x)
    df$z <- z
    
    mainEffects <- colnames(x)
    quadEffects <- getQuadraticTerms(x)
    
    formulaString <- paste0("z ~ ",
                        paste0(mainEffects, collapse = " + "),
                        " + ",
                        paste0(quadEffects, collapse = " + "))

    m <- bayesglm(formulaString, data = df, family = binomial())
    
    callingEnv$ps.z <- fitted(m)
  }
}

generateDataForIterInCurrentEnvironment <- function(iter, x, z, w, overlap = TRUE, covariates = "select", setting = "A", p.score = "none") {
  getQuadraticTerms <- function(x)
  {
    terms <- character()
    isBinary <- sapply(seq_len(ncol(x)), function(j) length(unique(x[,j])) == 2)
    for (i in seq_len(ncol(x))) {
      if (!isBinary[i]) terms[length(terms) + 1L] <- paste0("I(", colnames(x)[i], "^2)")
      
      if (i < ncol(x)) for (j in seq(i + 1L, ncol(x)))
        terms[length(terms) + 1L] <- paste0(colnames(x)[i], ":", colnames(x)[j])
    }
    terms
  }

  callingEnv <- parent.frame(1L)
  
  set.seed(iter * 5L + if (iter <= 500L) 565L else 7565L)
  
  if (is.numeric(covariates)) {
    x <- x[,sample(ncol(x), covariates)]
  } else if (covariates == "reduced") {
    x <- x[,sample(ncol(x), 5L)]
  }
  
  x.m <- if (is.data.frame(x)) dbarts::makeModelMatrixFromDataFrame(x) else x
  if (is.integer(x.m))
    x.m <- matrix(as.double(x.m), nrow(x.m), dimnames = dimnames(x.m))
  
  if (covariates == "junk") {
    temp <- cbind(x.m, matrix(rnorm(length(x.m)), nrow(x.m)))
    colnames(temp) <- c(colnames(x.m), paste0("x", seq_along(colnames(x.m))))
    callingEnv$x.r <- temp
  } else {
    callingEnv$x.r <- x.m
  }
  
  n <- nrow(x)
  p <- ncol(x) ## main effects
  
  if (setting == "B" || setting == "C") {
    mainEffects <- colnames(x)
    quadEffects <- getQuadraticTerms(x)
    if (covariates == "full" || (is.numeric(covariates) && covariates > 50)) {
      probs.q.0 <- max(1.0 - 5^(1 - 1 / 8) / (p^(15/16)), 0.4)
      quadEffects <- quadEffects[rbinom(length(quadEffects), 1, 1 - probs.q.0) == 1]
    }
    formulaString <- paste0("y ~ ",
                            paste0(mainEffects, collapse = " + "),
                            " + ",
                            paste0(quadEffects, collapse = " + "))
    
    temp <- x
    temp$y <- numeric(nrow(x))
    mod <- glm(formulaString, data = temp, x = TRUE)
    coefs <- mod$coef[-1L]
    x.m <- mod$x[,-1L]
    x.m <- x.m[,!is.na(coefs)]
  }
  
  x.m <- cbind(1.0, x.m)
  
  sigma.y <- 1.0
  tau <- 4.0
  
  w.full <- matrix(c(0.0, rep_len(w, ncol(x.m) - 1)), n, ncol(x.m), byrow = TRUE)
  
  if (setting == "B" || setting == "C") {
    vals.m  <- c(0.0, 1.0, 2.0)
    probs.m <- max(1.0 - 2.0 / sqrt(p), 0.2) ## hits 0.6 w/25 covariates, 0.8 w/108
    probs.m <- c(probs.m, 0.75 * (1.0 - probs.m), 0.25 * (1.0 - probs.m))
    
    if (covariates != "full" && (!is.numeric(covariates) || covariates <= 50)) {
      vals.q  <- c(0.0, 0.5, 1.0)
      probs.q <- max(1.0 - 5^(1 - 1 / 8) / (p^(15/16)), 0.4) ## hits 0.8 w/25 covariates, 0.95 w/108 
      probs.q <- c(probs.q, 0.75 * (1.0 - probs.q), 0.25 * (1.0 - probs.q))
    } else {
      vals.q  <- c(0.5, 1.0)
      probs.q <- c(0.75, 0.25)
    }
    
    beta.m0 <- sample(vals.m, p + 1, replace = TRUE, prob = probs.m)
    beta.m1 <- sample(vals.m, p + 1, replace = TRUE, prob = probs.m)
    beta.q0 <- sample(vals.q, ncol(x.m) - p - 1, replace = TRUE, prob = probs.q)
    beta.q1 <- sample(vals.q, ncol(x.m) - p - 1, replace = TRUE, prob = probs.q)
  } else {
    vals <- seq(0.0, 0.4, 0.1)
    probs <- max(1.0 - 2 / sqrt(p), 0.2)
    probs <- c(probs, rep(0.25 * (1.0 - probs), 4L))
    
    beta <- c(sample(seq(-1, 1, 0.25), 1), sample(vals, p, replace = TRUE, prob = probs))
  }
  
  if (setting == "B" || setting == "C") {
    mu.0 <- x.m %*% c(beta.m0, beta.q0)
    mu.1 <- x.m %*% c(beta.m1, beta.q1)
  } else {
    mu.0 <- exp((x.m + w.full) %*% beta)
    mu.1 <- x.m %*% beta
  }
  
  if (setting == "C") {
    ## old
    ## gamma <- sample(c(0,.1,.2), ncol(XXXmat), replace=TRUE, c(1/3,1/3,1/3))
    ## ps.z <- invlogit(ps.z %*% gamma - sqrt(ncol(ps.z)) * 5.7 / sqrt(303))
    
    mainEffects <- colnames(x)
    mainEffectColumns <- sample(ncol(x), max(4 * sqrt(ncol(x)) / sqrt(25), 2))
    mainEffects <- mainEffects[mainEffectColumns]
    x.p <- x[,mainEffectColumns]
    quadEffects <- getQuadraticTerms(x.p)
    quadEffects <- quadEffects[sample(length(quadEffects), min(max(5 * sqrt(ncol(x.p)) / sqrt(16), 2), length(quadEffects)))]
    
    formulaString <- paste0("y ~ ",
                            paste0(mainEffects, collapse = " + "),
                            " + ",
                            paste0(quadEffects, collapse = " + "))
    
    temp <- x
    temp$y <- numeric(nrow(x))
    mod <- glm(formulaString, data = temp, x = TRUE)
    coefs <- mod$coef[-1L]
    x.p <- mod$x[,-1L]
    x.p <- x.p[,!is.na(coefs)]
    
    invlogit <- function(x) { e.x <- exp(x); e.x / (1.0 + e.x) }
    
    gamma <- runif(ncol(x.p), -0.5, 0.5)
    lin.pred <- x.p %*% gamma
    lin.pred <- lin.pred - median(lin.pred) - 1.35
    ps.z <- invlogit(lin.pred)
    
    z <- rbinom(n, 1, ps.z)
    if (all(z == 1) || all(z == 0)) browser()
    callingEnv$ps.z <- ps.z
  }
    
  callingEnv$z.r <- z
  
  omega <- if (overlap == TRUE)
    mean(mu.1[z == 1] - mu.0[z == 1]) - tau
  else
    mean(mu.1[z == 0] - mu.0[z == 0]) - tau
    
  mu.1 <- mu.1 - omega

  
  y.0 <- rnorm(n, mu.0, sigma.y)
  y.1 <- rnorm(n, mu.1, sigma.y)
  y <- ifelse(z == 1, y.1, y.0)
  
  if (setting == "A") {
    f.0 <- function(x) exp((x + w) %*% beta)
    f.1 <- function(x) x %*% beta - omega
    environment(f.0) <- npci:::args2env(baseenv(), w = c(0, w), beta, omega)
    environment(f.1) <- environment(f.0)
    
    callingEnv$f.0  <- f.0
    callingEnv$f.1  <- f.1
  }
   
  callingEnv$mu.0 <- mu.0
  callingEnv$mu.1 <- mu.1
  callingEnv$y.0  <- y.0
  callingEnv$y.1  <- y.1
  callingEnv$y    <- y
  invisible(NULL)
}
