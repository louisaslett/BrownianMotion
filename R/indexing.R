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
    t <- bm$t[x]
  }

  if(length(x) == 0) {
    if(is.na(var)) {
      if(nrow(bm$layers) == 0) {
        lyrs.x <- NULL
      } else {
        lyrs.x <- (rowSums(sapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t })) > 0)
      }
      return(list(t = t,
                  W_t = rep(NA, length(t)),
                  W_tm = rep(NA, length(t)),
                  layers = bm$layers[lyrs.x,]))
    } else if(var == "t") {
      # Get just t
      return(unname(t))
    } else if(var == "W_t") {
      # Get just W_t
      return(rep(NA, length(t)))
    } else if(var == "W_tm") {
      # Get just W_tm
      return(rep(NA, length(t)))
    } else if(var == "layers") {
      # Get just layers
      if(nrow(bm$layers) == 0) {
        lyrs.x <- NULL
      } else {
        lyrs.x <- (rowSums(sapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t })) > 0)
      }
      return(bm$layers[lyrs.x,])
    }
  }

  # index in x, maybe var in t
  if(length(x)>0 && any(as.integer(x)!=x)) {
    stop("must be an integer index")
  }
  if(length(x)>0 && any(x < 1 | x > length(bm$t))) {
    stop("index out of bounds")
  }

  if(is.na(var) || var == "") {
    # Get everything

    if(nrow(bm$layers) == 0) {
      lyrs.x <- NULL
    } else {
      lyrs.x <- (rowSums(sapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t })) > 0)
    }

    res <- matrix(c(unname(bm$t[x]), missing,
                    unname(bm$W_t[x]), rep(NA, length(missing)),
                    unname(bm$W_tm[x]), rep(NA, length(missing))), ncol = 3)
    res <- res[order(res[,1]),,drop=FALSE][reorder.t,,drop=FALSE]

    return(list(t = res[,1],
                W_t = res[,2],
                W_tm = res[,3],
                layers = bm$layers[lyrs.x,]))
  } else if(var == "t") {
    # Get just t
    return(unname(t))
  } else if(var == "W_t") {
    # Get just W_t
    res <- matrix(c(unname(bm$t[x]), missing,
                    unname(bm$W_t[x]), rep(NA, length(missing))), ncol = 2)
    res <- res[order(res[,1]),,drop=FALSE][reorder.t,,drop=FALSE]
    return(unname(res[,2]))
  } else if(var == "W_tm") {
    # Get just W_tm
    res <- matrix(c(unname(bm$t[x]), missing,
                    unname(bm$W_tm[x]), rep(NA, length(missing))), ncol = 2)
    res <- res[order(res[,1]),,drop=FALSE][reorder.t,,drop=FALSE]
    return(unname(res[,2]))
  } else if(var == "layers") {
    # Get just layers
    if(nrow(bm$layers) == 0) {
      lyrs.x <- NULL
    } else {
      lyrs.x <- (rowSums(sapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t })) > 0)
    }
    return(bm$layers[lyrs.x,])
  }
}

