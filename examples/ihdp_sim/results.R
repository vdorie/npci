getPrefix <- function(method, overlap = TRUE, covariates = "full") {
  paste0(method, "_", if (overlap) "overlap" else "nonoverlap", "_", if (identical(covariates, "reduced")) "reduced" else "full")
}

## point to directory above "data" dir/root of sim folder
collateResults <- function(method, overlap = TRUE, covariates = "full", dir = ".", consolidate = FALSE) {
  prefix <- getPrefix(method, overlap, covariates)
  
  files <- list.files(file.path(dir, "data"), paste0(prefix, "_[0-9]+"))
  
  temp <- sapply(strsplit(files, "\\."), function(x) x[1L])
  temp <- sapply(strsplit(temp, "_"), function(x) as.integer(x[c(length(x) - 1L, length(x))]))
  start <- temp[1,]
  end   <- temp[2,]
  
  results.t <- matrix(NA_real_, max(end), 6L)
  colnames(results.t) <- c("bias", "cov", "cil", "wrong", "tau.est", "precision")
  
  for (i in seq_along(files)) {
    load(file.path(dir, "data", files[i]))
    resultsRange <- seq.int(start[i], end[i])
    results.t[resultsRange,] <- results
  }
  
  if (consolidate == TRUE) {
    unlink(file.path(dir, "data", files))
    
    numResults <- nrow(results.t)
    
    start <- 1L
    while (start <= numResults && any(is.na(results.t[start,]))) start <- start + 1L
    end  <- start + 1L
    
    while (start <= numResults) {
      while (end <= numResults && !any(is.na(results.t[end,]))) end <- end + 1L
      
      resultsRange <- seq.int(start, end - 1L)
      results <- results.t[resultsRange,]
      
      fileName <- paste0(prefix, "_", start, "_", end - 1L, ".RData")
      save(results, file = file.path(dir = ".", "data", fileName))
      
      start <- end
      while (start <= numResults && any(is.na(results.t[start,]))) start <- start + 1L
      end   <- start + 1L
    }
  }
  
  naRows <- apply(results, 1, function(row) any(is.na(row)))
  results.t <- results.t[!naRows,]
  
  results.t
}

getResultIntervals <- function(method, overlap = TRUE, covariates = "full", dir = ".")
{
  prefix <- getPrefix(method, overlap, covariates)
  files <- list.files(file.path(dir, "data"), paste0(prefix, "_[0-9]+"))
  
  resultNames <- list(NULL, c("start", "end"))
  
  if (length(files) == 0L) return(matrix(integer(), 0L, 2L, dimnames = resultNames))
  
  temp <- sapply(strsplit(files, "\\."), function(x) x[1L])
  temp <- sapply(strsplit(temp, "_"), function(x) as.integer(x[c(length(x) - 1L, length(x))]))
  start <- temp[1L,]
  end   <- temp[2L,]
  
  if (length(files) == 1L) return(matrix(c(start[1L], end[1L]), 1L, 2L, dimnames = resultNames))
  
  for (i in seq.int(2L, length(start))) {
    if (start[i] == end[i - 1L] - 1L)  {
      start[i] <- start[i - 1L]
      start[i - 1L] <- NA
      end[i - 1L] <- NA
    }
  }
  validRows <- !is.na(start)
  matrix(c(start[validRows], end[validRows]), sum(validRows), dimnames = resultNames)
}

## a and b are unions of intervals
intervalSubtraction <- function(a, b)
{
  if (NROW(b) == 0L) return(a)
  
  ## total cheese, but we just use ranges
  a.r <- seq.int(a[1L, "start"], a[1L, "end"])
  if (NROW(a) > 1L) for (i in 2L:nrow(a)) a.r <- c(a.r, seq.int(a[i, "start"], a[i, "end"]))
  
  b.r <- seq.int(b[1L, "start"], b[1L, "end"])
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
    end[m]   <- c[rh - 1L]
    lh <- rh
    rh <- lh + 1L
  }
  
  matrix(c(start, end), m, dimnames = list(NULL, c("start", "end")))
}
