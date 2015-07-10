dataFile <- file.path("data", "ihdp.RData")
if (!file.exists(dataFile)) stop("ihdp data file not available")

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

rm(dataFile, ihdp, covariateNames, trans)

generateDataForIterInCurrentEnvironment <- function(iter, x, z, w) {
  callingEnv <- parent.frame(1L)
  
  set.seed(iter * 5 + if (iter <= 500) 564 else 7565)
  
  x <- cbind(1, x)
  
  n <- nrow(x)
  p <- ncol(x)
  sigma.y <- 1
  w.full <- matrix(c(0, w), n, p, byrow = TRUE)
  
  beta <- sample(seq(0.0, 0.4, 0.1), p, replace = TRUE,
                 prob = c(0.6, rep(0.1, 4)))
  
  mu.0 <- exp((x + w.full) %*% beta)
  mu.1 <- x %*% beta
  
  omega <- mean(mu.1[z == 0] - mu.0[z == 0]) - 4
  mu.1 <- mu.1 - omega
  y.0 <- rnorm(n, mu.0, sigma.y)
  y.1 <- rnorm(n, mu.1, sigma.y)
  y <- ifelse(z == 1, y.1, y.0)
  
  f.0 <- function(x) exp((x + w) %*% beta)
  f.1 <- function(x) x %*% beta - omega
  environment(f.0) <- npci:::args2env(baseenv(), w = c(0, w), beta, omega)
  environment(f.1) <- environment(f.0)
  
  callingEnv$mu.0 <- mu.0
  callingEnv$mu.1 <- mu.1
  callingEnv$y.0  <- y.0
  callingEnv$y.1  <- y.1
  callingEnv$y    <- y
  callingEnv$f.0  <- f.0
  callingEnv$f.1  <- f.1
  invisible(NULL)
}
