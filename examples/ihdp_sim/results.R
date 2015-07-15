## point to directory above "data" dir/root of sim folder
collateResults <- function(method, setting = NULL, dir = ".", consolidate = FALSE) {
  prefix <- if (is.null(setting)) method else paste0(method, "_", setting)
  
  files <- list.files("data", paste0(prefix, "_[0-9]+"))
  
  temp <- sapply(strsplit(files, "\\."), function(x) x[1])
  temp <- sapply(strsplit(temp, "_"), function(x) as.integer(x[2:3]))
  start <- temp[1,]
  end   <- temp[2,]
  
  results.t <- matrix(NA_real_, max(end), 5)
  colnames(results.t) <- c("bias", "cov", "cil", "wrong", "tau.est")
  precision.t <- rep(NA_real_, max(end))
  
  for (i in seq_along(files)) {
    load(file.path(dir, "data", files[i]))
    resultsRange <- seq.int(start[i], end[i])
    results.t[resultsRange,] <- results
    precision.t[resultsRange] <- precision
  }
  
  if (consolidate == TRUE) {
    unlink(file.path(dir, "data", files))
    
    numResults <- length(precision.t)
    
    start <- 1
    while (start <= numResults && is.na(precision.t[start])) start <- start + 1
    end  <- start + 1
    
    while (start <= numResults) {
      while (end <= numResults && !is.na(precision.t[end])) end <- end + 1
      
      resultsRange <- seq.int(start, end - 1)
      results <- results.t[resultsRange,]
      precision <- precision.t[resultsRange]
      
      fileName <- paste0(prefix, "_", start, "_", end - 1, ".RData")
      save(results, precision, file = file.path(dir = ".", "data", fileName))
      
      start <- end
      while (start <= numResults && is.na(precision.t[start])) start <- start + 1
      end   <- start + 1
    }
  }
  
  naRows <- is.na(precision.t)
  results.t <- results.t[!naRows,]
  precision.t <- precision.t[!naRows]
  
  return(list(results = results.t, precision = precision.t))
}

getResultIntervals <- function(method, setting = NULL, dir = ".")
{
  prefix <- if (is.null(setting)) method else paste0(method, "_", setting)
  files <- list.files("data", paste0(prefix, "_[0-9]+"))
  
  resultNames <- list(NULL, c("start", "end"))
  
  if (length(files) == 0) return(matrix(integer(), 0, 2, dimnames = resultNames))
  temp <- sapply(strsplit(files, "\\."), function(x) x[1])
  temp <- sapply(strsplit(temp, "_"), function(x) as.integer(x[2:3]))
  start <- temp[1,]
  end   <- temp[2,]
  
  if (length(files) == 1) return(matrix(c(start[1], end[1]), 1, 2, dimnames = resultNames))
  
  for (i in seq.int(2, length(start))) {
    if (start[i] == end[i - 1] - 1)  {
      start[i] <- start[i - 1]
      start[i - 1] <- NA
      end[i - 1] <- NA
    }
  }
  validRows <- !is.na(start)
  matrix(c(start[validRows], end[validRows]), sum(validRows), dimnames = resultNames)
}

## a and b are unions of intervals
intervalSubtraction <- function(a, b)
{
  if (NROW(b) == 0) return(a)
  
  ## total cheese, but we just use ranges
  a.r <- seq.int(a[1L, "start"], a[1L, "end"])
  if (NROW(a) > 1L) for (i in 2L:nrow(a)) a.r <- c(a.r, seq.int(a[i, "start"], a[i, "end"]))
  
  b.r <- seq.int(b[1, "start"], b[1, "end"])
  if (NROW(b) > 1L) for (i in 2L:nrow(b)) b.r <- c(b.r, seq.int(b[i, "start"], b[i, "end"]))
  
  c <- sort(setdiff(a.r, b.r))
  if (length(c) == 0L) return(matrix(integer(), 0L, 2L, dimnames = list(NULL, c("start", "end"))))
  
  m <- 0L
  start <- numeric()
  end   <- numeric()

  lh <- 1L
  rh <- 2L
  n <- length(c)
  
  while (lh <= n) {
    while (rh <= n && c[rh] == c[rh - 1L] + 1L) rh <- rh + 1L
    m <- m + 1L
    start[m] <- c[lh]
    end[m]   <- c[rh - 1]
    lh <- rh
    rh <- lh + 1L
  }
  
  matrix(c(start, end), m, dimnames = list(NULL, c("start", "end")))
}
