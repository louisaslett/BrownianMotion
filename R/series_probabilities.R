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

  # miss <-
  #   + (exp(P*(Uu-x)*(Uu-y)) + exp(P*(-Ll+x)*(-Ll+y))) -
  #   (exp(P*(Uu-x)*(Uu-y)) + exp(P*(-Lu+x)*(-Lu+y))) -
  #   (exp(P*(Ul-x)*(Ul-y)) + exp(P*(-Ll+x)*(-Ll+y))) +
  #   (exp(P*(Ul-x)*(Ul-y)) + exp(P*(-Lu+x)*(-Lu+y)))

  s1 <-
    + (exp(P*Dlu^2-P*Dlu*z) + exp(P*Dlu^2+P*Dlu*z)) - # Remainder of z1 after cancellation
    (exp(P*Dll^2-P*Dll*z) + exp(P*Dll^2+P*Dll*z)) -  # Remainder of z2 after cancellation
    (exp(P*Duu^2-P*Duu*z) + exp(P*Duu^2+P*Duu*z)) +  # Remainder of z3 after cancellation
    (exp(P*Dul^2-P*Dul*z) + exp(P*Dul^2+P*Dul*z)) -  # Remainder of z4 after cancellation
    sum(exp(P*(D1lu-x)*(D1lu-y)) + exp(P*(D2lu+x)*(D2lu+y)) - exp(P*j^2*Dlu^2-P*j*Dlu*z) - exp(P*j^2*Dlu^2+P*j*Dlu*z)) + # Remaining terms of z1 *A
    sum(exp(P*(D1ll-x)*(D1ll-y)) + exp(P*(D2ll+x)*(D2ll+y)) - exp(P*j^2*Dll^2-P*j*Dll*z) - exp(P*j^2*Dll^2+P*j*Dll*z)) + # Remaining terms of z2 *C
    sum(exp(P*(D1uu-x)*(D1uu-y)) + exp(P*(D2uu+x)*(D2uu+y)) - exp(P*j^2*Duu^2-P*j*Duu*z) - exp(P*j^2*Duu^2+P*j*Duu*z)) - # Remaining terms of z3 *D
    sum(exp(P*(D1ul-x)*(D1ul-y)) + exp(P*(D2ul+x)*(D2ul+y)) - exp(P*j^2*Dul^2-P*j*Dul*z) - exp(P*j^2*Dul^2+P*j*Dul*z)) - # Remaining terms of z4 *B
    (exp(P*((m+1)/2)^2*Dlu^2-P*((m+1)/2)*Dlu*z) + exp(P*((m+1)/2)^2*Dlu^2+P*((m+1)/2)*Dlu*z)) - # We want a lower bound on z1 so removing the final terns from *A
    (exp(P*((m+1)/2)^2*Dul^2-P*((m+1)/2)*Dul*z) + exp(P*((m+1)/2)^2*Dul^2+P*((m+1)/2)*Dul*z)) # We want a lower bound on z4 so removing the final terns from *B

  s2 <- s1 +
    (exp(P*((m+1)/2)^2*Dlu^2-P*((m+1)/2)*Dlu*z) + exp(P*((m+1)/2)^2*Dlu^2+P*((m+1)/2)*Dlu*z)) + # Now a upper bound on z1 by re-introducing removed terms from *A
    (exp(P*((m+1)/2)^2*Dul^2-P*((m+1)/2)*Dul*z) + exp(P*((m+1)/2)^2*Dul^2+P*((m+1)/2)*Dul*z)) + # Now a upper bound on z4 by re-introducing removed terms from *B
    (exp(P*(Dll*((m+3)/2)+Ll-x)*(Dll*((m+3)/2)+Ll-y)) + exp(P*(Dll*((m+3)/2)-Ul+x)*(Dll*((m+3)/2)-Ul+y))) + # Now an upper bound for z2 at next integer evaluation following the form of *C
    (exp(P*(Duu*((m+3)/2)+Lu-x)*(Duu*((m+3)/2)+Lu-y)) + exp(P*(Duu*((m+3)/2)-Uu+x)*(Duu*((m+3)/2)-Uu+y))) # Now an upper bound for z3 at next integer evaluation following the form of *D

  c(s1=s1, s2=s2) # s1 is the lower bound, s2 is the upper bound
}



