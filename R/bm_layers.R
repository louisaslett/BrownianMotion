# Checks if x contains a valid layer index
is.layeridx <- function(bm, x) {
  checkmate::test_integerish(x,
                             lower = 1,
                             upper = n.lyr_(bm),
                             any.missing = FALSE,
                             min.len = 1,
                             unique = TRUE)
}

# Checks if x contains valid layer (bound part only) information
is.layerbounds <- function(x) {
  is_tibble(x) &
    identical(names(x), c("Ld", "Uu", "Lu", "Ud", "Lu.hard", "Ud.hard")) &
    identical(sapply(x, typeof) , c(Ld = "double", Uu = "double", Lu = "double", Ud = "double", Lu.hard = "logical", Ud.hard = "logical"))
}

is.layer <- function(x) {
  is_tibble(x) &
    identical(names(x), c("type", "t.l", "t.u", "Ld", "Uu", "Lu", "Ud", "Lu.hard", "Ud.hard")) &
    identical(sapply(x, typeof) , c(type = "integer", t.l = "double", t.u = "double", Ld = "double",
                                    Uu = "double", Lu = "double", Ud = "double", Lu.hard = "logical",
                                    Ud.hard = "logical")) &
    is.timevector(x$t.l) &
    is.timevector(x$t.u) &
    all(x$t.u > x$t.l)
}

# Retrieve layer index and info for any layers covering a vector of times in t
# NB: the interval is half open [) favouring t.l
get.lyr_ <- function(bm, t) {
  if(!is.bm(bm))
    stop("get.lyr_ received invalid bm argument")
  if(!is.timevector(t))
    stop("get.lyr_ received invalid t argument")

  # outer provides a row per observation in t, with a column for each layer,
  # indicating whether the condition is met
  # apply checks if the layer covers any of the provided times
  idx <- which(apply(outer(t, bm$layers$t.l, ">=") & outer(t, bm$layers$t.u, "<"),
                     2,
                     any))

  if(length(idx) == 0) {
    idx <- NULL
    layer <- NULL
  } else {
    layer <- bm$layers[idx,]
  }

  list(idx = idx,
       layer = layer)
}

# Return number of layers in the Brownian motion
n.lyr_ <- function(bm) {
  if(!is.bm(bm))
    stop("get.lyr_ received invalid bm argument")

  nrow(bm$layers)
}

# Eliminate layer(s) by index
rm.lyr_ <- function(bm, idx) {
  if(!is.bm(bm))
    stop("get.lyr_ received invalid bm argument")
  if(!is.layeridx(bm, idx))
    stop("get.lyr_ received invalid idx argument")

  bm$layers <- bm$layers[-idx,]

  NULL
}

# Update layer(s) by idx with bound information in new
update.lyr_ <- function(bm, idx, new) {
  if(!is.bm(bm))
    stop("update.lyr_ received invalid bm argument")
  if(!is.layeridx(bm, idx))
    stop("update.lyr_ received invalid idx argument")
  if(!is.layerbounds(new))
    stop("update.lyr_ received invalid new argument")

  if(length(idx) != nrow(new))
    stop("update.lyr_ differing number of layers and new layer rows provided")

  bm$layers[idx,4:9] <- new

  NULL
}

split.lyr_ <- function(bm, idx, l, r) {
  if(!is.bm(bm))
    stop("split.lyr_  received invalid bm argument")
  if(!is.layeridx(bm, idx))
    stop("split.lyr_  received invalid idx argument")
  if(!is.layer(l))
    stop("split.lyr_  received invalid l argument")
  if(!is.layer(r))
    stop("split.lyr_  received invalid r argument")

  ## NOTE TO LOUIS: MAKE THIS WORK FOR VECTOR INDEX AND MATRIX L AND R
  ## BY USING IDENTICAL() AND THEN LOOPED INSERTION

  if(l$t.l != bm$layers$t.l[idx])
    stop("split.lyr_ lower time of left layer in split does not match layer to be replaced")
  if(r$t.u != bm$layers$t.u[idx])
    stop("split.lyr_ upper time of right layer in split does not match layer to be replaced")
  if(l$t.u != r$t.l)
    stop("split.lyr_ replacement left/right layers do not meet")

  bm$layers <- rbind(bm$layers[1:(idx-1),],
                     l,
                     r,
                     bm$layers[(idx+1):nrow(bm$layers),])

  NULL
}

add.lyr_ <- function() {

}
