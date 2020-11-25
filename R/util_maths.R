mpfrthr <- 9 # m threshold for high precision
mpfrpbn <- 10

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
  c(s1=s1,s2=s2)
}

##########

eazeta_ <- function(n, s, t, x, y, L, U) {
  if(max(x-U, y-U, L-x, L-y) >= 0) {
    1
  } else {
    j <- 1:(ceiling(n/2))
    P <- -2/(t-s)
    D <- U-L
    D1 <- D*j + L
    D2 <- D*j - U
    z <- y-x
    if(n%%2 == 0) {
      sum(exp(P*(D1-x)*(D1-y)) +
            exp(P*(D2+x)*(D2+y)) -
            exp(P*j^2*D^2-P*j*D*z) -
            exp(P*j^2*D^2+P*j*D*z))
    } else {
      sum(exp(P*(D1-x)*(D1-y)) +
            exp(P*(D2+x)*(D2+y))) -
        sum(exp(P*j[1:length(j)-1]^2*D^2 - P*j[1:length(j)-1]*D*z) +
              exp(P*j[1:length(j)-1]^2*D^2 + P*j[1:length(j)-1]*D*z))
    }
  }
}

eagamma_ <- function(n, s, t, x, y, L, U) {
  1 - eazeta_(n, s, t, x, y, L, U)
}

eapsi_ <- function(j, s, t, m, xoy, u) {
  P <- -2*abs(u-m)*j/(t-s)
  (2*abs(u-m)*j - (xoy-m)) * exp(P*(abs(u-m)*j - (xoy-m)))
}

eachi_ <- function(j, s, t, m, xoy, u) {
  P <- -2*abs(u-m)*j/(t-s)
  (2*abs(u-m)*j + (xoy-m)) * exp(P*(abs(u-m)*j + (xoy-m)))
}

eadel2_ <- function(n, s, t, m, xoy, u) {
  if(max(xoy-u, m-xoy) >= 0) {
    0
  } else {
    if(n%%2 == 0) {
      j <- 1:(n/2)
      1-(sum(eapsi_(j, s, t, m, xoy, u) - eachi_(j, s, t, m, xoy, u)))/(xoy-m)
    } else {
      if(n>1) {
        j <- 1:((n-1)/2)
        1-(sum(eapsi_(j, s, t, m, xoy, u) - eachi_(j, s, t, m, xoy, u)))/(xoy-m) - eapsi_(max(j)+1, s, t, m, xoy, u)/(xoy-m)
      } else {
        1-eapsi_(1, s, t, m, xoy, u)/(xoy-m)
      }
    }
  }
}

eadelR_ <- function(n, s, t, x, y, m, u) {
  if(x == m) {
    xI <- 1
  } else {
    xI <- 0
  }
  if(y == m) {
    yI <- 1
  } else {
    yI <- 0
  }
  if(max(xI, yI) == 1) {
    delT <- 2
  } else {
    delT <- 1
  }
  if(m > min(x,y)) {
    x <- -x
    y <- -y
    m <- -m
    u <- -u
  }
  if(max(x-u, y-u, m-x, m-y) >= 0) {
    out <- 0
  }
  if(delT == 1) {
    out <- eagamma_(n, s, t, x, y, m, u)/(1 - exp(-2*(x-m)*(y-m)/(t-s)))
  }
  if(delT == 2) {
    if(xI*yI == 0) {
      xoy <- max(x,y)
      out <- eadel2_(n, s, t, m, xoy, u)
    } else {
      out <- 0
    }
  }
  if(out < 0) {
    out <- 0
  }
  if(out > 1) {
    out <- 1
  }
  if((t-s) == 0) {
    out <- 1
  }
  out
}

eadelC_ <- function(mt, s, t, x, y, m, u) {
  if(mt >= mpfrthr) {
    pbn <- mt*mpfrpbn
    s <- mpfr(s, precBits = pbn)
    t <- mpfr(t, precBits = pbn)
    x <- mpfr(x, precBits = pbn)
    y <- mpfr(y, precBits = pbn)
    m <- mpfr(m, precBits = pbn)
    u <- mpfr(u, precBits = pbn)
    c(s1 = eadelR_(mt, s, t, x, y, m, u),
      s2 = eadelR_(mt+1, s, t, x, y, m, u))
  } else {
    scale:::eadel_pair_cpp(mt,s,t,x,y,m,u)
  }
}

