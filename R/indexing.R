#' @export
`[.BrownianMotion` <- function(bm, x, t = NULL) {
  var <- as.character(sys.call()[-1])[3]
  if(!is.numeric(t) && is.character(x)) {
    # label in x, maybe var in t
    x <- x[x %in% names(bm$labels)]
    t <- c(bm$labels[x], recursive = TRUE)
  }
  missing <- NULL
  if(is.numeric(t)) {
    reorder.t <- order(order(t))
    missing <- t[!(t %in% bm$t)]
    x <- which(bm$t %in% t)
  } else {
    reorder.t <- order(x)
  }

  if(length(x) == 0) {
    if(is.na(var)) {
      return(list(t = NULL,
                  W_t = NULL,
                  W_tm = NULL,
                  layers = bm$layers[c(),]))
    } else if(var == "layers") {
      # Get just layers
      return(bm$layers[c(),])
    } else {
      return(NULL)
    }
  }

  # index in x, maybe var in t
  if(length(x)>0 && any(as.integer(x)!=x)) {
    stop("must be an integer index")
  }
  if(length(x)>0 && any(x < 1 | x > length(bm$t))) {
    stop("index out of bounds")
  }

  if(is.na(var)) {
    # Get everything

    lyrs.x <- (rowSums(sapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t })) > 0)

    res <- matrix(c(bm$t[x], missing,
                    bm$W_t[x], rep(NA, length(missing)),
                    bm$W_tm[x], rep(NA, length(missing))), ncol = 3)
    res <- res[order(res[,1]),,drop=FALSE][reorder.t,,drop=FALSE]

    return(list(t = res[,1],
                W_t = res[,2],
                W_tm = res[,3],
                layers = bm$layers[lyrs.x,]))
  } else if(var == "t") {
    # Get just t
    return(bm$t[x])
  } else if(var == "W_t") {
    # Get just W_t
    return(bm$W_t[x])
  } else if(var == "W_tm") {
    # Get just W_tm
    return(bm$W_tm[x])
  } else if(var == "layers") {
    # Get just layers
    lyrs.x <- (rowSums(sapply(bm$t[x], function(t) { bm$layers$t.l <= t & bm$layers$t.u > t })) > 0)
    return(bm$layers[lyrs.x,])
  }
}

