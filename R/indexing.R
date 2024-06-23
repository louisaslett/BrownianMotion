#' @export
`[.BrownianMotion` <- function(bm, x = NULL, t = NULL) {
  # Determine what variable has been requested (eg bm[1,"W_t"])
  # or is NA if not specified (eg bm[1])
  var <- as.character(sys.call()[-1])[3]

  if(is.null(x)) {
    x <- 1:length(bm$t)
  }

  if(!is.null(t)) {
    t2 <- substitute(t)
    env <- parent.frame()
    if(as.character(t2)[1] == ":") {
      t <- c(eval(t2[[2]], envir = env),
             bm$t[bm$t > eval(t2[[2]], envir = env) &
                  bm$t < eval(t2[[3]], envir = env)],
             eval(t2[[3]], envir = env))
    } else {
      t <- eval(t, envir = env)
    }
  }

  # Do we have label in x, and maybe variable to retrieve in t?
  # If so, get time for each label (and dump variable from t, as we now have in var)
  if(!is.numeric(t) && is.character(x)) {
    x <- x[x %in% names(bm$labels)]
    t <- c(bm$labels[x], recursive = TRUE)
  }

  # Create dummy to hold times that are not in the skeleton for which
  # we will only return the relevant layer info (if any)
  missing <- NULL

  # Do we have a time in t?
  # If so, find index and put in x,
  # otherwise, x must be an index already
  if(is.numeric(t)) {
    t <- unique(t)
    reorder.t <- order(order(t))
    missing <- t[!(t %in% bm$t)]
    x <- which(bm$t %in% t)
  } else {
    reorder.t <- order(x)
    t <- bm$t[x]
  }

  # Deal with the case where none of the times are skeleton points
  if(length(x) == 0) {
    if(is.na(var)) {
      # Get everything
      if(nrow(bm$layers) == 0) {
        lyrs.x <- NULL
      } else {
        # 20/11/23: If there is exactly one layer, the next line is a bug, fixed (I think) by the line after
        #lyrs.x <- (rowSums(sapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t })) > 0)
        lyrs.x <- (rowSums(simplify2array(lapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t }), except = 0)) > 0)
      }
      return(list(t = t,
                  W_t = rep(NA, length(t)),
                  W_tm = rep(NA, length(t)),
                  layers = bm$layers[lyrs.x,],
                  labels = list()))
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
        # 20/11/23: If there is exactly one layer, the next line is a bug, fixed (I think) by the line after
        #lyrs.x <- (rowSums(sapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t })) > 0)
        lyrs.x <- (rowSums(simplify2array(lapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t }), except = 0)) > 0)
      }
      return(bm$layers[lyrs.x,])
    } else if(var == "labels") {
      return(list())
    }
  }

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
      # Extract layers at both skeleton times and non-skeleton times
      # 20/11/23: If there is exactly one layer, the next line is a bug, fixed (I think) by the line after
      #lyrs.x <- (rowSums(as.matrix(sapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t }))) > 0)
      lyrs.x <- (rowSums(simplify2array(lapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t }), except = 0)) > 0)
    }

    res <- matrix(c(unname(bm$t[x]), missing,
                    unname(bm$W_t[x]), rep(NA, length(missing)),
                    unname(bm$W_tm[x]), rep(NA, length(missing))), ncol = 3)
    res <- res[order(res[,1]),,drop=FALSE][reorder.t,,drop=FALSE]
    res.labels <- lapply(bm$labels, function(ts) { ts[ts %in% res[,1]] })

    return(list(t = res[,1],
                W_t = res[,2],
                W_tm = res[,3],
                layers = bm$layers[lyrs.x,],
                labels = res.labels[lapply(res.labels, length)>0]))
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
      # 20/11/23: If there is exactly one layer, the next line is a bug, fixed (I think) by the line after
      #lyrs.x <- (rowSums(sapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t })) > 0)
      lyrs.x <- (rowSums(simplify2array(lapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t }), except = 0)) > 0)
    }
    return(bm$layers[lyrs.x,])
  } else if(var == "labels") {
    res <- matrix(c(unname(bm$t[x]), missing,
                    unname(bm$W_t[x]), rep(NA, length(missing))), ncol = 2)

    res.labels <- lapply(bm$labels, function(ts) { ts[ts %in% res[,1]] })

    return(res.labels[lapply(res.labels, length)>0])
  }
}



#' @export
`[.BrownianMotionNd` <- function(bm, x = NULL, t = NULL) {
  # Determine what variable has been requested (eg bm[1,"W_t"])
  # or is NA if not specified (eg bm[1])
  var <- as.character(sys.call()[-1])[3]

  if(is.null(x)) {
    x <- 1:length(bm$Z.bm[[1]]$t)
  }

  if(!is.null(t)) {
    t2 <- substitute(t)
    env <- parent.frame()
    if(as.character(t2)[1] == ":") {
      t <- c(eval(t2[[2]], envir = env),
             bm$Z.bm[[1]]$t[bm$Z.bm[[1]]$t > eval(t2[[2]], envir = env) &
                            bm$Z.bm[[1]]$t < eval(t2[[3]], envir = env)],
             eval(t2[[3]], envir = env))
    } else {
      t <- eval(t, envir = env)
    }
  }

  # Do we have label in x, and maybe variable to retrieve in t?
  # If so, get time for each label (and dump variable from t, as we now have in var)
  if(!is.numeric(t) && is.character(x)) {
    # TODO: improve this, we're assuming first dimension labels are same as all
    #       other dimensions. May not be true for internally generated labels
    #       from layer routines
    x <- x[x %in% names(bm$Z.bm[[1]]$labels)]
    t <- c(bm$Z.bm[[1]]$labels[x], recursive = TRUE)
  }

  # Create dummy to hold times that are not in the skeleton for which
  # we will only return the relevant layer info (if any)
  missing <- NULL

  # Do we have a time in t?
  # If so, find index and put in x,
  # otherwise, x must be an index already
  if(is.numeric(t)) {
    t <- unique(t)
    reorder.t <- order(order(t))
    # TODO: improve this? This gives correct result, but ...
    #       We've got $t duplicated in all Z.bm objects, but this
    #       makes implementation easier (especially for layers where we fill in
    #       anything missing across dims after simulation)
    missing <- t[!(t %in% bm$Z.bm[[1]]$t)]
    x <- which(bm$Z.bm[[1]]$t %in% t)
  } else {
    reorder.t <- order(x)
    t <- bm$Z.bm[[1]]$t[x]
  }

  # Deal with the case where none of the times are skeleton points
  if(length(x) == 0) {
    if(is.na(var)) {
      if(nrow(bm$Z.bm[[1]]$layers) == 0) {
        # TODO: Murray says he'd love to spend some time allowing
        #       layers in only certain dimensions and for the
        #       path to be unconstrained in others.
        #       Louis places face in palms
        # TODO: Louis thinks it would be lovely to have a rainbow of
        #       colours for different time segments and overlapping transparent
        #       polygons and infilled paths.
        #       Murray places face in palms
        lyrs <- tibble(
          t.l = numeric(),
          t.u = numeric(),
          outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
          outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
          inner.cube = list()
        )
      } else {
        lyrs <- dplyr::bind_rows(lapply(c(bm$Z.bm[[1]]$t[x], missing), function(t) { transformlayersNd_(bm, t) }))
      }
      return(list(t = t,
                  W_t = rep(NA, length(t)),
                  W_tm = rep(NA, length(t)),
                  layers = lyrs,
                  labels = list()))
    } else if(var == "t") {
      # Get just t
      return(unname(t))
    } else if(var == "W_t") {
      # Get just W_t
      return(matrix(NA, nrow = length(t), ncol = bm$dim))
    } else if(var == "W_tm") {
      # Get just W_tm
      return(matrix(NA, nrow = length(t), ncol = bm$dim))
    } else if(var == "layers") {
      # Get just layers
      if(nrow(bm$Z.bm[[1]]$layers) == 0) {
        lyrs <- tibble(
          t.l = numeric(),
          t.u = numeric(),
          outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
          outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
          inner.cube = list()
        )
      } else {
        lyrs <- dplyr::bind_rows(lapply(c(bm$Z.bm[[1]]$t[x], missing), function(t) { transformlayersNd_(bm, t) }))
      }
      return(lyrs)
    } else if(var == "labels") {
      return(list())
    }
  }

  # index in x, maybe var in t
  if(length(x)>0 && any(as.integer(x)!=x)) {
    stop("must be an integer index")
  }
  if(length(x)>0 && any(x < 1 | x > length(bm$Z.bm[[1]]$t))) {
    stop("index out of bounds")
  }

  if(is.na(var) || var == "") {
    # Get everything

    if(nrow(bm$Z.bm[[1]]$layers) == 0) {
      lyrs <- tibble(
        t.l = numeric(),
        t.u = numeric(),
        outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
        outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
        inner.cube = list()
      )
    } else {
      lyrs <- dplyr::bind_rows(lapply(c(bm$Z.bm[[1]]$t[x], missing), function(t) { transformlayersNd_(bm, t) }))
    }

    # res <- matrix(c(unname(bm$t[x]), missing,
    #                 unname(bm$W_t[x]), rep(NA, length(missing)),
    #                 unname(bm$W_tm[x]), rep(NA, length(missing))), ncol = 3)
    # res <- res[order(res[,1]),,drop=FALSE][reorder.t,,drop=FALSE]

    res.t <- c(unname(bm$Z.bm[[1]]$t[x]), missing)
    res.W_t <- rbind(getWtNd_(bm, x), matrix(NA, nrow = length(missing), ncol = bm$dim))
    res.W_tm <- rbind(getWtmNd_(bm, x), matrix(NA, nrow = length(missing), ncol = bm$dim))
    res.labels <- lapply(bm$Z.bm[[1]]$labels, function(ts) { ts[ts %in% res.t] })

    o <- order(res.t)
    return(list(t = res.t[o][reorder.t],
                W_t = res.W_t[o,,drop=FALSE][reorder.t,,drop=FALSE],
                W_tm = res.W_tm[o,,drop=FALSE][reorder.t,,drop=FALSE],
                layers = lyrs,
                labels = res.labels[lapply(res.labels, length)>0]))
  } else if(var == "t") {
    # Get just t
    return(unname(t))
  } else if(var == "W_t") {
    # Get just W_t
    res.t <- c(unname(bm$Z.bm[[1]]$t[x]), missing)
    res.W_t <- rbind(getWtNd_(bm, x), matrix(NA, nrow = length(missing), ncol = bm$dim))

    o <- order(res.t)
    return(res.W_t[o,,drop=FALSE][reorder.t,,drop=FALSE])
  } else if(var == "W_tm") {
    # Get just W_tm
    res.t <- c(unname(bm$Z.bm[[1]]$t[x]), missing)
    res.W_tm <- rbind(getWtmNd_(bm, x), matrix(NA, nrow = length(missing), ncol = bm$dim))

    o <- order(res.t)
    return(res.W_tm[o,,drop=FALSE][reorder.t,,drop=FALSE])
  } else if(var == "layers") {
    # Get just layers
    if(nrow(bm$Z.bm[[1]]$layers) == 0) {
      lyrs <- tibble(
        t.l = numeric(),
        t.u = numeric(),
        outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
        outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
        inner.cube = list()
      )
    } else {
      lyrs <- dplyr::bind_rows(lapply(c(bm$Z.bm[[1]]$t[x], missing), function(t) { transformlayersNd_(bm, t) }))
    }
    return(lyrs)
  } else if(var == "labels") {
    res.t <- c(unname(bm$Z.bm[[1]]$t[x]), missing)
    res.labels <- lapply(bm$Z.bm[[1]]$labels, function(ts) { ts[ts %in% res.t] })

    return(res.labels[lapply(res.labels, length)>0])
  }
}



#' @export
`[.BrownianMotionNdZ` <- function(bm, x = NULL, t = NULL) {
  bm <- bm[[1]]

  # Determine what variable has been requested (eg bm[1,"W_t"])
  # or is NA if not specified (eg bm[1])
  var <- as.character(sys.call()[-1])[3]

  if(is.null(x)) {
    x <- 1:length(bm$Z.bm[[1]]$t)
  }

  if(!is.null(t)) {
    t2 <- substitute(t)
    env <- parent.frame()
    if(as.character(t2)[1] == ":") {
      t <- c(eval(t2[[2]], envir = env),
             bm$Z.bm[[1]]$t[bm$Z.bm[[1]]$t > eval(t2[[2]], envir = env) &
                              bm$Z.bm[[1]]$t < eval(t2[[3]], envir = env)],
             eval(t2[[3]], envir = env))
    } else {
      t <- eval(t, envir = env)
    }
  }

  # Do we have label in x, and maybe variable to retrieve in t?
  # If so, get time for each label (and dump variable from t, as we now have in var)
  if(!is.numeric(t) && is.character(x)) {
    # TODO: improve this, we're assuming first dimension labels are same as all
    #       other dimensions. May not be true for internally generated labels
    #       from layer routines
    x <- x[x %in% names(bm$Z.bm[[1]]$labels)]
    t <- c(bm$Z.bm[[1]]$labels[x], recursive = TRUE)
  }

  # Create dummy to hold times that are not in the skeleton for which
  # we will only return the relevant layer info (if any)
  missing <- NULL

  # Do we have a time in t?
  # If so, find index and put in x,
  # otherwise, x must be an index already
  if(is.numeric(t)) {
    t <- unique(t)
    reorder.t <- order(order(t))
    # TODO: improve this? This gives correct result, but ...
    #       We've got $t duplicated in all Z.bm objects, but this
    #       makes implementation easier (especially for layers where we fill in
    #       anything missing across dims after simulation)
    missing <- t[!(t %in% bm$Z.bm[[1]]$t)]
    x <- which(bm$Z.bm[[1]]$t %in% t)
  } else {
    reorder.t <- order(x)
    t <- bm$Z.bm[[1]]$t[x]
  }

  # Deal with the case where none of the times are skeleton points
  if(length(x) == 0) {
    if(is.na(var)) {
      if(nrow(bm$Z.bm[[1]]$layers) == 0) {
        # TODO: Murray says he'd love to spend some time allowing
        #       layers in only certain dimensions and for the
        #       path to be unconstrained in others.
        #       Louis places face in palms
        # TODO: Louis thinks it would be lovely to have a rainbow of
        #       colours for different time segments and overlapping transparent
        #       polygons and infilled paths.
        #       Murray places face in palms
        lyrs <- tibble(
          t.l = numeric(),
          t.u = numeric(),
          outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
          outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
          inner.cube = list()
        )
      } else {
        lyrs <- dplyr::bind_rows(lapply(c(bm$Z.bm[[1]]$t[x], missing), function(t) { transformlayersNd_(bm, t, do.transform = FALSE) }))
      }
      return(list(t = t,
                  Z_t = rep(NA, length(t)),
                  Z_tm = rep(NA, length(t)),
                  layers = lyrs,
                  labels = list()))
    } else if(var == "t") {
      # Get just t
      return(unname(t))
    } else if(var == "Z_t") {
      # Get just Z_t
      return(matrix(NA, nrow = length(t), ncol = bm$dim))
    } else if(var == "Z_tm") {
      # Get just Z_tm
      return(matrix(NA, nrow = length(t), ncol = bm$dim))
    } else if(var == "layers") {
      # Get just layers
      if(nrow(bm$Z.bm[[1]]$layers) == 0) {
        lyrs <- tibble(
          t.l = numeric(),
          t.u = numeric(),
          outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
          outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
          inner.cube = list()
        )
      } else {
        lyrs <- dplyr::bind_rows(lapply(c(bm$Z.bm[[1]]$t[x], missing), function(t) { transformlayersNd_(bm, t, do.transform = FALSE) }))
      }
      return(lyrs)
    } else if(var == "labels") {
      return(list())
    }
  }

  # index in x, maybe var in t
  if(length(x)>0 && any(as.integer(x)!=x)) {
    stop("must be an integer index")
  }
  if(length(x)>0 && any(x < 1 | x > length(bm$Z.bm[[1]]$t))) {
    stop("index out of bounds")
  }

  if(is.na(var) || var == "") {
    # Get everything

    if(nrow(bm$Z.bm[[1]]$layers) == 0) {
      lyrs <- tibble(
        t.l = numeric(),
        t.u = numeric(),
        outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
        outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
        inner.cube = list()
      )
    } else {
      lyrs <- dplyr::bind_rows(lapply(c(bm$Z.bm[[1]]$t[x], missing), function(t) { transformlayersNd_(bm, t, do.transform = FALSE) }))
    }

    res.t <- c(unname(bm$Z.bm[[1]]$t[x]), missing)
    res.Z_t <- rbind(getZtNd_(bm, x), matrix(NA, nrow = length(missing), ncol = bm$dim))
    res.Z_tm <- rbind(getZtmNd_(bm, x), matrix(NA, nrow = length(missing), ncol = bm$dim))
    res.labels <- lapply(bm$Z.bm[[1]]$labels, function(ts) { ts[ts %in% res.t] })

    o <- order(res.t)
    return(list(t = res.t[o][reorder.t],
                Z_t = res.Z_t[o,,drop=FALSE][reorder.t,,drop=FALSE],
                Z_tm = res.Z_tm[o,,drop=FALSE][reorder.t,,drop=FALSE],
                layers = lyrs,
                labels = res.labels[lapply(res.labels, length)>0]))
  } else if(var == "t") {
    # Get just t
    return(unname(t))
  } else if(var == "Z_t") {
    # Get just Z_t
    res.t <- c(unname(bm$Z.bm[[1]]$t[x]), missing)
    res.Z_t <- rbind(getZtNd_(bm, x), matrix(NA, nrow = length(missing), ncol = bm$dim))

    o <- order(res.t)
    return(res.Z_t[o,,drop=FALSE][reorder.t,,drop=FALSE])
  } else if(var == "Z_tm") {
    # Get just Z_tm
    res.t <- c(unname(bm$Z.bm[[1]]$t[x]), missing)
    res.Z_tm <- rbind(getZtmNd_(bm, x), matrix(NA, nrow = length(missing), ncol = bm$dim))

    o <- order(res.t)
    return(res.Z_tm[o,,drop=FALSE][reorder.t,,drop=FALSE])
  } else if(var == "layers") {
    # Get just layers
    if(nrow(bm$Z.bm[[1]]$layers) == 0) {
      lyrs <- tibble(
        t.l = numeric(),
        t.u = numeric(),
        outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
        outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
        inner.cube = list()
      )
    } else {
      lyrs <- dplyr::bind_rows(lapply(c(bm$Z.bm[[1]]$t[x], missing), function(t) { transformlayersNd_(bm, t, do.transform = FALSE) }))
    }
    return(lyrs)
  } else if(var == "labels") {
    res.t <- c(unname(bm$Z.bm[[1]]$t[x]), missing)
    res.labels <- lapply(bm$Z.bm[[1]]$labels, function(ts) { ts[ts %in% res.t] })

    return(res.labels[lapply(res.labels, length)>0])
  }
}
