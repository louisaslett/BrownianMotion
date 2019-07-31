mpfrthr <- 9 # m threshold for high precision

eagammaC_ <- function(m, s, t, x, y, L, U) {
  if(m >= mpfrthr) {
    pbn <- m*mpfrpbn
    s <- mpfr(s,precBits=pbn)
    t <- mpfr(t,precBits=pbn)
    x <- mpfr(x,precBits=pbn)
    y <- mpfr(y,precBits=pbn)
    L <- mpfr(L,precBits=pbn)
    U <- mpfr(U,precBits=pbn)
  }
  z <- eazetaC_(m,s,t,x,y,L,U)
  c(s1 = as.numeric(1-z[1]), s2 = as.numeric(1-z[2]))
}

eazetaC_ <- function(m, s, t, x, y, L, U) {
  if(max(x-U,y-U,L-x,L-y) >= 0) {
    s1<-s2<-1
  } else {
    j <- 1:((m+1)/2)
    P <- -2/(t-s)
    D <- U-L
    D1 <- D*j+L
    D2 <- D*j-U
    z <- y-x
    s2 <- sum(exp(P*(D1-x)*(D1-y)) +
                exp(P*(D2+x)*(D2+y)) -
                exp(P*j^2*D^2-P*j*D*z) -
                exp(P*j^2*D^2+P*j*D*z))
    s1 <- s2 + exp(P*((m+1)/2)^2*D^2 - P*((m+1)/2)*D*z) +
      exp(P*((m+1)/2)^2*D^2 + P*((m+1)/2)*D*z)
  }
  c(s1=s2,s2=s2)
}
