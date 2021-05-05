mpfrthr <- 100 # m threshold for high precision
mpfrpbn <- 5

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
  c(s1 = as.numeric(1-z[1]), s2 = as.numeric(1-z[2])) # s1 is the lower bound, s2 is the upper bound
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
  c(s1=s1,s2=s2) # s1 is the upper bound, s2 is the lower bound
}

##########


# eabetaC_old		<- function(m,s,t,x,y,Ll,Lu,Ul,Uu) {
#   if(m>=mpfrthr) {
#     pbn <- m*mpfrpbn
#     s <- mpfr(s,precBits=pbn)
#     t <- mpfr(t,precBits=pbn)
#     x <- mpfr(x,precBits=pbn)
#     y <- mpfr(y,precBits=pbn)
#     Ll <- mpfr(Ll,precBits=pbn)
#     Lu <- mpfr(Lu,precBits=pbn)
#     Ul <- mpfr(Ul,precBits=pbn)
#     Uu <- mpfr(Uu,precBits=pbn)
#   }
#   z1<-eazetaC_(m,s,t,x,y,Ll,Uu)
#   z2<-c(eazetaC_(m,s,t,x,y,Lu,Uu)[2],eazetaC_(m+2,s,t,x,y,Lu,Uu)[1])
#   z3<-c(eazetaC_(m,s,t,x,y,Ll,Ul)[2],eazetaC_(m+2,s,t,x,y,Ll,Ul)[1])
#   z4<-eazetaC_(m,s,t,x,y,Lu,Ul)
#   c(s1=as.numeric(-z1[1]+z2[1]+z3[1]-z4[1]),s2=as.numeric(-z1[2]+z2[2]+z3[2]-z4[2])) # s1 is the lower bound, s2 is the upper bound
# }




##########

earhoC_ <- function(m,s,q,t,x,w,y,Ll,Lu,Ul,Uu) {
  beta.l.p0 <- eabetaC_(m,s,q,x,w,Ll,Lu,Ul,Uu)
  beta.l.p2 <- eabetaC_(m+2,s,q,x,w,Ll,Lu,Ul,Uu)
  beta.r.p0 <- eabetaC_(m,q,t,w,y,Ll,Lu,Ul,Uu)
  beta.r.p2 <- eabetaC_(m+2,q,t,w,y,Ll,Lu,Ul,Uu)

  z1L.n <- c(beta.l.p0[3],beta.l.p0[4],beta.l.p2[3],beta.l.p2[4])
  z1R.n <- c(beta.r.p0[3],beta.r.p0[4],beta.r.p2[3],beta.r.p2[4])
  z3L.n <- c(beta.l.p0[5],beta.l.p0[6],beta.l.p2[5],beta.l.p2[6])
  z3R.n <- c(beta.r.p0[5],beta.r.p0[6],beta.r.p2[5],beta.r.p2[6])
  z2L.n <- c(beta.l.p0[7],beta.l.p0[8],beta.l.p2[7],beta.l.p2[8])
  z2R.n <- c(beta.r.p0[7],beta.r.p0[8],beta.r.p2[7],beta.r.p2[8])
  z4L.n <- c(beta.l.p0[9],beta.l.p0[10],beta.l.p2[9],beta.l.p2[10])
  z4R.n <- c(beta.r.p0[9],beta.r.p0[10],beta.r.p2[9],beta.r.p2[10])

  # s2.new2 <- beta.l.p0[1] + beta.r.p0[1] + z1L.n[2]*z1R.n[2] - z2L.n[1]*z2R.n[1] - z3L.n[1]*z3R.n[1] +z4L.n[2]*z4R.n[2] # Note, perhaps worth modifying the upper bound to remove the double function calls to beta.

  # s2 lower, s1 upper
  c(s1=(beta.l.p0[2] + beta.r.p0[2] + z1L.n[1]*z1R.n[1] - z2L.n[2]*z2R.n[2] - z3L.n[2]*z3R.n[2] +z4L.n[1]*z4R.n[1]),
    s2=(beta.l.p2[1] + beta.r.p2[1] + z1L.n[2]*z1R.n[2] - z2L.n[3]*z2R.n[3] - z3L.n[3]*z3R.n[3] +z4L.n[2]*z4R.n[2]),
    s3=(beta.l.p2[2] + beta.r.p2[2] + z1L.n[3]*z1R.n[3] - z2L.n[4]*z2R.n[4] - z3L.n[4]*z3R.n[4] +z4L.n[3]*z4R.n[3]))
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
  }
  c(s1 = eadelR_(mt, s, t, x, y, m, u),
    s2 = eadelR_(mt+1, s, t, x, y, m, u))
}

###############

#' @importFrom statmod rinvgauss
eaminmax_ <- function(s,t,x,y,Ll,Lu,Ul,Uu){
  if(rbinom(1,1,0.5)==1){minI<-1}else{minI<--1; x<--x;y<--y;Ll<--Uu;Lu<--Ul}
  e.a <- -2*(Ll-x)*(Ll-y)/(t-s); e.b <- -2*(Lu-x)*(Lu-y)/(t-s)
  u1<-runif(1,0,1); u2<-runif(1,0,1); m<-x-0.5*(sqrt((y-x)^2-2*(t-s)*(e.a+log((1-u1)+u1*(exp(e.b-e.a)))))-(y-x))
  c1 <- (y-m)^2/(2*(t-s)); c2 <- (m-x)^2/(2*(t-s)); I1<-rinvgauss(1,sqrt(c1/c2),2*c1);I2<-1/rinvgauss(1,sqrt(c2/c1),2*c2);V<-if(runif(1,0,1)<=(1/(1+sqrt(c1/c2)))){I1}else{I2}; tau <- (s*V+t)/(1+V)
  if(minI==-1){m<--m};list(m=m,tau=tau,minI=minI)}

eabesmid_ <- function(q,s,tau,t,x,m,y,minI){
  if(minI==-1){x<--x;y<--y;m<--m} # Reflection
  bI<-0;if(q==s){bI<-1;w<-x};if(q==tau){bI<-1;w<-m};if(q==t){bI<-1;w<-y} # Boundary
  t<-t-s;tau<-tau-s;q<-q-s # Rescale time
  if(bI==0){if(q<tau){Ra1<-sqrt(tau);Ra2<-(x-m)*(tau-q)/((tau)^(3/2));Ra3<-(tau-q)/(tau)}else{Ra1<-sqrt(t-tau);Ra2<-(y-m)*(q-tau)/((t-tau)^(3/2));Ra3<-(q-tau)/(t-tau)};BB3<-rnorm(3,0,sqrt(Ra3*(1-Ra3)));w<-m+Ra1*sqrt((Ra2+BB3[1])^2+(BB3[2])^2+(BB3[3])^2)}
  list(w=minI*w)}

matsort_ <- function(mat,n) {mat[rank(mat[,n]),]<- mat[c(1:nrow(mat)),];return(mat)}

eabesex_ <- function(sV,tV,xV,yV,m,B1,B2,minI){ # Vectors of equal length
  if(minI==-1){xV<--xV;yV<--yV;m<--m;B1<--B1;B2<--B2}
  u <- runif(1,0,1); mt<-3;  em<-matrix(0,length(sV),8); em[,1]<-sV; em[,2]<-tV; em[,3]<-xV; em[,4]<-yV
  B1evI<-B2evI<-0; while(B1evI==0){for(i in 1:dim(em)[1]){em[i,5:6]<-eadelC_(mt,em[i,1],em[i,2],em[i,3],em[i,4],m,B1)};if(u<=prod(em[,5])){B1evI<-B2evI<-1;con1I<-1;con2I<-1;ex1I<-0;ex2I<-0}else{if(u>prod(em[,6])){B1evI<-1;con1I<-0;ex1I<-1}else{B1evI<-0;con1I<-0;ex1I<-0;mt<-mt+2}}}
  while(B2evI==0){for(i in 1:dim(em)[1]){em[i,7:8]<-eadelC_(mt,em[i,1],em[i,2],em[i,3],em[i,4],m,B2)};if(u<=prod(em[,7])){B2evI<-1;con2I<-1;ex1I<-0}else{if(u>prod(em[,8])){B2evI<-1;con2I<-0;ex2I<-1}else{B2evI<-0;con2I<-0;ex2I<-0;mt<-mt+2}}}
  if(minI==-1){em[,3]<--em[,3];em[,4]<--em[,4]}; accI<-0; if(con1I==1){accI<-1}else{if(con2I==1){if(rbinom(1,1,0.5)==1){accI<-1}}}
  list(accI=accI,u=u,con1I=con1I,con2I=con2I,em=em)}
