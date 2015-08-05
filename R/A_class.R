setClass("gpci",
         slots = 
         list(data     = "list",
              trans    = "list",
              env      = "environment"),
         
         prototype =
         list(data     = list(y = numeric(0), x = matrix(NA_real_, 0, 0), z = numeric(0), x.mean = matrix(NA_real_, 0, 0)),
              trans    = list(),
              env      = baseenv()),
         
         validity = function(object) {
           if (length(object@data$y) != NROW(object@data$x) || length(object@data$y) != length(object@data$z)) return("length of 'y' must match that of 'z', and number of rows in 'x'")
           TRUE
         }
)
