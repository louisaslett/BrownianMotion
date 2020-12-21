#' Refine an intersection layer
#'
#' @param
#'
#'
#'
refine <- function() {

}

refine.intersection_ <- function(bm, lyr.idx, mult) {
  s <- bm$layers[lyr.idx,"t.l",drop=TRUE]
  t <- bm$layers[lyr.idx,"t.u",drop=TRUE]
  x <- bm$W_t[match(s, bm$t)]
  y <- bm$W_t[match(t, bm$t)]
  Ll <- bm$layers[lyr.idx,"Ld",drop=TRUE]
  Lu <- bm$layers[lyr.idx,"Lu",drop=TRUE]
  Ul <- bm$layers[lyr.idx,"Ud",drop=TRUE]
  Uu <- bm$layers[lyr.idx,"Uu",drop=TRUE]

  if(max((Uu-Ul), (Lu-Ll)) > mult*sqrt(t-s)) {
    while(max((Uu-Ul), (Lu-Ll)) > mult*sqrt(t-s)) {
      Ls<-Ll+(Lu-Ll)/2
      Us<-Uu-(Uu-Ul)/2
      mat <- matrix(c(Ll,Ls,Us,Uu,Ls,Lu,Us,Uu,Ll,Ls,Ul,Us,Ls,Lu,Ul,Us), 4, 4, byrow=TRUE)
      bbind <- deind <- 0
      bbct <- 1
      m1 <- 3
      u1 <- runif(1, 0, 1)
      while(bbind == 0) {
        while(deind == 0) {
          debd <- eabe3C(m1,s,t,x,y,Ll,Lu,Ul,Uu)[2:3]
          if(debd[2] <= 0) {
            m1 <- m1+2
          } else {
            deind<-1
          }
        }
        be <- matrix(0, bbct, 2)
        for(i in 1:bbct) {
          be[i,] <- eabetaC(m1,s,t,x,y,mat[i,1],mat[i,2],mat[i,3],mat[i,4])
        }
        bd <- c(sum(be[,1])/debd[1],sum(be[,2])/debd[2])
        if(u1 <= bd[1]) {
          bbind <- 1
        } else {
          if(u1 >= bd[2]) {
            bbct <- bbct+1
            deind <- 0
          } else {
            m1 <- m1+2
            deind <- 0
          }
        }
      }
      Ll <- mat[bbct,1]
      Lu <- mat[bbct,2]
      Ul <- mat[bbct,3]
      Uu <- mat[bbct,4]
    }
  }

  bm$layers[lyr.idx,"Ld"] <- Ll
  bm$layers[lyr.idx,"Lu"] <- Lu
  bm$layers[lyr.idx,"Ud"] <- Ul
  bm$layers[lyr.idx,"Uu"] <- Uu

  bm
}
