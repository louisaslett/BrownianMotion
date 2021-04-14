easiga <- function(P,z,A,L) {
  P*(A*(A+2*L-z)+L*(L-z))
}

easigb <- function(P,z,A,L) {
  P*(z-A-L)
}

eaphia <- function(P,z,A,si) {
  P*(A^2-si*A*z)
}

eaphib <- function(P,z,A,si) {
  P*si*A
}

eabetaC_ <- function(m,s,t,x,y,Ll,Lu,Ul,Uu) {
  if(Ll==Lu | Ul==Uu | min(x,y) < Lu | max(x,y) > Ul){
    return(c(s1=0,s2=0))
  }

  j <- 2:((m+1)/2) # Evaluated from second term in infinite sum to account for cancellation of terms in the first term
  P <- -2/(t-s)
  z <- y-x
  Dlu <- Uu-Ll
  D1lu <- Dlu*j+Ll
  D2lu <- Dlu*j-Uu
  Dll <- Ul-Ll
  D1ll <- Dll*j+Ll
  D2ll <- Dll*j-Ul
  Duu <- Uu-Lu
  D1uu <- Duu*j+Lu
  D2uu <- Duu*j-Uu
  Dul <- Ul-Lu
  D1ul <- Dul*j+Lu
  D2ul <- Dul*j-Ul

  s1res <- (exp(P*(Uu-x)*(Uu-y)) + exp(P*(-Ll+x)*(-Ll+y)))
  s2res <- - (exp(P*(Ul-x)*(Ul-y)) + exp(P*(-Ll+x)*(-Ll+y)))
  s3res <- - (exp(P*(Uu-x)*(Uu-y)) + exp(P*(-Lu+x)*(-Lu+y)))
  s4res <- (exp(P*(Ul-x)*(Ul-y)) + exp(P*(-Lu+x)*(-Lu+y)))

  # miss <-   s1res + s2res + s3res + s4res

  s1p1 <- + (exp(P*Dlu^2-P*Dlu*z) + exp(P*Dlu^2+P*Dlu*z)) - sum(exp(P*(D1lu-x)*(D1lu-y)) + exp(P*(D2lu+x)*(D2lu+y)) - exp(P*j^2*Dlu^2-P*j*Dlu*z) - exp(P*j^2*Dlu^2+P*j*Dlu*z)) - (exp(P*((m+1)/2)^2*Dlu^2-P*((m+1)/2)*Dlu*z) + exp(P*((m+1)/2)^2*Dlu^2+P*((m+1)/2)*Dlu*z))
  s1p2 <- - (exp(P*Dll^2-P*Dll*z) + exp(P*Dll^2+P*Dll*z)) + sum(exp(P*(D1ll-x)*(D1ll-y)) + exp(P*(D2ll+x)*(D2ll+y)) - exp(P*j^2*Dll^2-P*j*Dll*z) - exp(P*j^2*Dll^2+P*j*Dll*z))
  s1p3 <- - (exp(P*Duu^2-P*Duu*z) + exp(P*Duu^2+P*Duu*z)) + sum(exp(P*(D1uu-x)*(D1uu-y)) + exp(P*(D2uu+x)*(D2uu+y)) - exp(P*j^2*Duu^2-P*j*Duu*z) - exp(P*j^2*Duu^2+P*j*Duu*z))
  s1p4 <- + (exp(P*Dul^2-P*Dul*z) + exp(P*Dul^2+P*Dul*z)) - sum(exp(P*(D1ul-x)*(D1ul-y)) + exp(P*(D2ul+x)*(D2ul+y)) - exp(P*j^2*Dul^2-P*j*Dul*z) - exp(P*j^2*Dul^2+P*j*Dul*z)) - (exp(P*((m+1)/2)^2*Dul^2-P*((m+1)/2)*Dul*z) + exp(P*((m+1)/2)^2*Dul^2+P*((m+1)/2)*Dul*z))

  s2p1 <- s1p1 + (exp(P*((m+1)/2)^2*Dlu^2-P*((m+1)/2)*Dlu*z) + exp(P*((m+1)/2)^2*Dlu^2+P*((m+1)/2)*Dlu*z))
  s2p2 <- s1p2 + (exp(P*(Dll*((m+3)/2)+Ll-x)*(Dll*((m+3)/2)+Ll-y)) + exp(P*(Dll*((m+3)/2)-Ul+x)*(Dll*((m+3)/2)-Ul+y)))
  s2p3 <- s1p3 + (exp(P*(Duu*((m+3)/2)+Lu-x)*(Duu*((m+3)/2)+Lu-y)) + exp(P*(Duu*((m+3)/2)-Uu+x)*(Duu*((m+3)/2)-Uu+y)))
  s2p4 <- s1p4 + (exp(P*((m+1)/2)^2*Dul^2-P*((m+1)/2)*Dul*z) + exp(P*((m+1)/2)^2*Dul^2+P*((m+1)/2)*Dul*z))

  s1z1 <- s1res - s1p1
  s1z2 <- - (s2res - s2p2)
  s1z3 <- - (s3res - s2p3)
  s1z4 <- s4res - s1p4

  s2z1 <- s1res - s2p1
  s2z2 <- - (s2res - s1p2)
  s2z3 <- - (s3res - s1p3)
  s2z4 <- s4res - s2p4

  s1 <- s1p1 + s1p2 + s1p3 + s1p4
  s2 <- s2p1 + s2p2 + s2p3 + s2p4

  c(s1=s1, s2=s2, s1z1=s1z1, s2z1=s2z1, s1z2=s1z2, s2z2=s2z2, s1z3=s1z3, s2z3=s2z3, s1z4=s1z4, s2z4=s2z4) # s1 is the lower bound, s2 is the upper bound
}


