#' Simulate Brownian motion conditional on intersection layer
#'
#'
#' @export
sim.condintersection <- function(bm, s, t, q = NULL, q.grid = NULL) {
  # Arg types & combos
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(as.integer(is.null(q)) + as.integer(is.null(q.grid)) != 1) {
    stop("Exactly one of q or q.grid must be specified")
  }
  if(!is.null(q.grid) && !is.intscalar(q.grid)) {
    stop("q.grid must be either NULL or an integer number of grid points.")
  }
  if(!is.null(q) && !is.realvector(q)) {
    stop("q must be either NULL or a vector of times.")
  }
  if(!is.realscalar(t)) {
    stop("t must be a scalar.")
  }
  if(!is.realscalar(s)) {
    stop("s must be a scalar.")
  }

  # Checks
  if(!(s %in% bm$t)) {
    stop("s not found in bm path.")
  }
  if(!(t %in% bm$t)) {
    stop("t not found in bm path.")
  }
  intersection <- bm$layers[bm$layers$type == "intersection",]
  if(!(t %in% intersection$t.u)) {
    stop("t is not the time of a realised maximum/minimum in a intersection layer.")
  }
  # Do we then need to check that s is the start of this same intersection layer?
  # Actually this is taken care of by the following which allows no intermediate obs
  s_idx <- match(s, bm$t)
  t_idx <- match(t, bm$t)
  if(t_idx != s_idx+1) {
    stop("Cannot have other brownian motion observations between times s and t")
  }

  # Make grid if it wasn't supplied
  if(is.null(q)) {
    q <- head(tail(seq(s, t, length.out = q.grid+2), -1), -1)
  }
  if(any(q < s) || any(q > t)) {
    stop("All q must lie between s and t")
  }
  if(is.unsorted(q)) {
    q <- sort(q)
  }
  # Eliminate times we know
  q <- setdiff(q, bm$t)

  for(qq in q) {
    bm.res <- sim.condintersection_(bm, s_idx, qq, t_idx)
    s_idx <- s_idx+1
    t_idx <- t_idx+1
  }

  bm.res$layers <- bm.res$layers[order(bm.res$layers$t.l),]
  invisible(bm.res)
}



sim.condintersection_ <- function(bm, s_idx, q, t_idx) {
  s <- bm$t[s_idx]
  t <- bm$t[t_idx]
  cur.layer <- which(bm$layers$t.l == s & bm$layers$t.u == t)

  # Simulate new point conditional on intersection layer
  w <- sim.condintersection.simpt_(s, q, t, bm$W_t[s_idx], bm$W_t[t_idx],
                                   Ll = bm$layers$Ld[cur.layer],
                                   Lu = bm$layers$Lu[cur.layer],
                                   Ul = bm$layers$Ud[cur.layer],
                                   Uu = bm$layers$Uu[cur.layer])
  bm$t <- c(bm$t[1:s_idx],
            q,
            bm$t[t_idx:length(bm$t)])
  bm$W_t <- c(bm$W_t[1:s_idx],
              w,
              bm$W_t[t_idx:length(bm$W_t)])

  # Update layer info
  # We split the layer in two, adding left and right layers either side of the
  # newly simulated time, then remove the old layer
  newlyrs <- sim.condintersection.simlyr_(s, q, t, bm$W_t[s_idx], w, bm$W_t[t_idx+1],
                                          Ll = bm$layers$Ld[cur.layer],
                                          Lu = bm$layers$Lu[cur.layer],
                                          Ul = bm$layers$Ud[cur.layer],
                                          Uu = bm$layers$Uu[cur.layer])

  bm$layers <- add_row(bm$layers,
                       type = "intersection",
                       t.l = s,
                       t.u = q,
                       Ld = newlyrs[1,5],
                       Uu = newlyrs[1,8],
                       Lu = newlyrs[1,6],
                       Ud = newlyrs[1,7],
                       Lu.hard = TRUE,
                       Ud.hard = TRUE)
  bm$layers <- add_row(bm$layers,
                       type = "intersection",
                       t.l = q,
                       t.u = t,
                       Ld = newlyrs[2,5],
                       Uu = newlyrs[2,8],
                       Lu = newlyrs[2,6],
                       Ud = newlyrs[2,7],
                       Lu.hard = TRUE,
                       Ud.hard = TRUE)
  bm$layers <- bm$layers[-cur.layer,]

  bm
}



sim.condintersection.simpt_ <- function(s, q, t, x, y, Ll, Lu, Ul, Uu) {

  even <- function(x) { x%%2 == 0 }
  odd <- function(x) { x%%2 == 1 }
  brutepbn<-200
  brutesig	<- function(j,s,t,x,y,L,U,side){ # Calculate sigma
    pbn<-brutepbn;s<-mpfr(s,precBits=pbn);t<-mpfr(t,precBits=pbn);x<-mpfr(x,precBits=pbn);y<-mpfr(y,precBits=pbn);L<-mpfr(L,precBits=pbn);U<-mpfr(U,precBits=pbn)
    if(side==-1){a<-(-2/(t-s))*((abs(U-L)*j)^2+2*abs(U-L)*j*min(L,U)+(min(L,U))^2-abs(U-L)*j*x-min(L,U)*x); b<-(-2/(t-s))*(-abs(U-L)*j-min(L,U)+x)}
    if(side==1){a<-(-2/(t-s))*((abs(U-L)*j)^2+2*abs(U-L)*j*min(L,U)+(min(L,U))^2-abs(U-L)*j*y-min(L,U)*y);b<-(-2/(t-s))*(-abs(U-L)*j-min(L,U)+y)}
    list(a=as.numeric(a),b=as.numeric(b))}
  brutetau	<- function(j,s,t,x,y,L,U,side){ # Calculate tau
    pbn<-brutepbn;s<-mpfr(s,precBits=pbn);t<-mpfr(t,precBits=pbn);x<-mpfr(x,precBits=pbn);y<-mpfr(y,precBits=pbn);L<-mpfr(L,precBits=pbn);U<-mpfr(U,precBits=pbn)
    if(side==-1){a<-(-2*j/(t-s))*(j*abs(U-L)^2+abs(U-L)*x);b<-(-2*j/(t-s))*(-abs(U-L))}
    if(side==1){a<-(-2*j/(t-s))*(j*abs(U-L)^2-abs(U-L)*y);b<-(-2*j/(t-s))*(abs(U-L))}
    list(a=as.numeric(a),b=as.numeric(b))}
  brutezet	<- function(n,s,t,x,y,L,U,side){ # Calculate zeta
    zetm<-matrix(0,0,3); colnames(zetm)<-c("sign","a","b")
    if(even(n)){for(i in 1:(n/2)){spos<-brutesig(i,s,t,x,y,L,U,side); sneg<-brutesig(i,s,t,-x,-y,-L,-U,side); tpos<-brutetau(i,s,t,x,y,L,U,side); tneg<-brutetau(i,s,t,-x,-y,-L,-U,side); zetm<-rbind(zetm,c(1,spos$a,spos$b),c(1,sneg$a,-sneg$b),c(-1,tpos$a,tpos$b),c(-1,tneg$a,-tneg$b))}}
    if(odd(n)){if(n>1){for(i in 1:((n-1)/2)){spos<-brutesig(i,s,t,x,y,L,U,side); sneg<-brutesig(i,s,t,-x,-y,-L,-U,side); tpos<-brutetau(i,s,t,x,y,L,U,side); tneg<-brutetau(i,s,t,-x,-y,-L,-U,side); zetm<-rbind(zetm,c(1,spos$a,spos$b),c(1,sneg$a,-sneg$b),c(-1,tpos$a,tpos$b),c(-1,tneg$a,-tneg$b))}};spos<-brutesig((n+1)/2,s,t,x,y,L,U,side); sneg<-brutesig((n+1)/2,s,t,-x,-y,-L,-U,side); zetm<-rbind(zetm,c(1,spos$a,spos$b),c(1,sneg$a,-sneg$b))}
    list(zetm=zetm)}
  brutezetmut <- function(n,s,q,t,x,y,L,U){w <- 1 # Dummy w	# Calculate zeta multiples
  zetm<-matrix(0,0,3); colnames(zetm)<-c("sign","a","b"); zetmL<-brutezet(n,s,q,x,w,L,U,-1)$zetm; zetmR<-brutezet(n,q,t,w,y,L,U,1)$zetm; for(i in 1:(dim(zetmL)[1])){for(j in 1:(dim(zetmR)[1])){zetm<-rbind(zetm,c(zetmL[i,"sign"]*zetmR[j,"sign"],zetmL[i,"a"]+zetmR[j,"a"],zetmL[i,"b"]+zetmR[j,"b"]))}}
  list(zetm=zetm)}
  bruterho	<- function(n,s,q,t,x,y,Ll,Lu,Ul,Uu){w <- 1 # Dummy w	# Calculate rho
  gamm1<-matrix(0,0,3); zetmL1<-brutezet(n+1,s,q,x,w,Ll,Uu,-1)$zetm; zetmL1[,"sign"]<--zetmL1[,"sign"]; zetmR1<-brutezet(n+1,q,t,w,y,Ll,Uu,1)$zetm; zetmR1[,"sign"]<--zetmR1[,"sign"]; zetmLR1<-brutezetmut(n,s,q,t,x,y,Ll,Uu)$zetm; gamm1<-rbind(gamm1,zetmL1,zetmR1,zetmLR1); gamm1<-cbind(gamm1,rep(Ll,dim(gamm1)[1]),rep(Uu,dim(gamm1)[1]))
  gamm2<-matrix(0,0,3); zetmL2<-brutezet(n,s,q,x,w,Lu,Uu,-1)$zetm; zetmR2<-brutezet(n,q,t,w,y,Lu,Uu,1)$zetm; zetmLR2<-brutezetmut(n+1,s,q,t,x,y,Lu,Uu)$zetm; zetmLR2[,"sign"]<--zetmLR2[,"sign"]; gamm2<-rbind(gamm2,zetmL2,zetmR2,zetmLR2); gamm2<-cbind(gamm2,rep(Lu,dim(gamm2)[1]),rep(Uu,dim(gamm2)[1]))
  gamm3<-matrix(0,0,3); zetmL3<-brutezet(n,s,q,x,w,Ll,Ul,-1)$zetm; zetmR3<-brutezet(n,q,t,w,y,Ll,Ul,1)$zetm; zetmLR3<-brutezetmut(n+1,s,q,t,x,y,Ll,Ul)$zetm; zetmLR3[,"sign"]<--zetmLR3[,"sign"]; gamm3<-rbind(gamm3,zetmL3,zetmR3,zetmLR3); gamm3<-cbind(gamm3,rep(Ll,dim(gamm3)[1]),rep(Ul,dim(gamm3)[1]))
  gamm4<-matrix(0,0,3); zetmL4<-brutezet(n+1,s,q,x,w,Lu,Ul,-1)$zetm; zetmL4[,"sign"]<--zetmL4[,"sign"]; zetmR4<-brutezet(n+1,q,t,w,y,Lu,Ul,1)$zetm; zetmR4[,"sign"]<--zetmR4[,"sign"]; zetmLR4<-brutezetmut(n,s,q,t,x,y,Lu,Ul)$zetm; gamm4<-rbind(gamm4,zetmL4,zetmR4,zetmLR4); gamm4<-cbind(gamm4,rep(Lu,dim(gamm4)[1]),rep(Ul,dim(gamm4)[1]))
  rhom<-matrix(0,0,5); rhom<-rbind(rhom,gamm1,gamm2,gamm3,gamm4); colnames(rhom)<-c("sign","a","b","L","U")
  list(rhom=rhom,gamm1=gamm1,gamm2=gamm2,gamm3=gamm3,gamm4=gamm4)}
  highevalab	<- function(mat,w){ # Calculate rho at point
    rmat<-mat[mat[,5]>=w,,drop=FALSE]; rmat<-rmat[rmat[,4]<w,,drop=FALSE]
    if(dim(rmat)[1]==0){ev<-0}else{pbn<-ceiling(100+10*(log(max(abs(rmat[,2]+rmat[,3]*w)))/log(2))); s<-mpfr(rmat[,1],pbn); a<-mpfr(rmat[,2],pbn);b<-mpfr(rmat[,3],pbn);w<-mpfr(w,pbn); ev<-sum(s*exp(a+b*w))}
    c(ev=as.numeric(ev))}
  brutelip	<- function(nmu,bbsd,L,U){ # Calculate normal lipschitz constant in interval
    conI<-0; if(min(U-(nmu+bbsd),(nmu+bbsd)-L)>=0){conI<-1}; if(min(U-(nmu-bbsd),(nmu-bbsd)-L)>=0){conI<-1}
    if(conI==1){lip<-abs((1/bbsd)*exp(-1/2))}else{lip<-max(abs(((nmu-L)/(bbsd^2))*exp(-(1/(2*bbsd^2))*(L-nmu)^2)),abs(((nmu-U)/(bbsd^2))*exp(-(1/(2*bbsd^2))*(U-nmu)^2)))}
    c(lip=lip)}
  brutenett1 	<- function(mat,bbm){ # Nett matrix type 1
    nmat<-matrix(0,0,dim(mat)[2]); dmat<-mat; colnames(nmat)<-colnames(dmat)<-colnames(mat)
    while(dim(dmat)[1]>0){
      eqmat<-rbind(1:dim(dmat)[1],mapply(function(ub,b){if(isTRUE(all.equal(ub,b))==TRUE){1}else{0}},rep(dmat[1,"b"],dim(dmat)[1]),dmat[,"b"])); commat<-dmat[eqmat[1,eqmat[2,]==1],,drop=FALSE]; dmat<-dmat[setdiff(1:dim(dmat)[1],eqmat[1,eqmat[2,]==1]),,drop=FALSE]
      pbn<-brutepbn; sm<-mpfr(commat[,"sign"],pbn); am<-mpfr(commat[,"a"],pbn); bm<-mpfr(commat[1,"b"],pbn); Lm<-mpfr(commat[1,"L"],pbn); Um<-mpfr(commat[1,"U"],pbn); bbmm<-mpfr(bbm,pbn); bbsdm<-mpfr(commat[1,"sd"],pbn);
      wnew<-sum(sm*exp(am+bm*bbmm+0.5*(bm*bbsdm)^2)); if(wnew!=0){
        snew<-sign(wnew); anew<-log(wnew*snew)-bm*bbmm-0.5*(bm*bbsdm)^2; trnew<-pnorm(Um,mean=bbmm+bm*bbsdm^2,sd=bbsdm)-pnorm(Lm,mean=bbmm+bm*bbsdm^2,sd=bbsdm); trwnew<-wnew*trnew; lipnew<-brutelip(bbmm+bm*bbsdm^2,bbsdm,Lm,Um)
        nmat<-rbind(nmat,c(snew,as.numeric(anew),commat[1,"b"],commat[1,"L"],commat[1,"U"],as.numeric(wnew),as.numeric(trnew),as.numeric(trwnew),commat[1,"mean"],commat[1,"sd"],as.numeric(lipnew)))}}
    list(nmat=nmat)}
  brutenett2 	<- function(mat){ # Nett matrix type 2
    nmat<-matrix(0,0,dim(mat)[2]); dmat<-mat; colnames(nmat)<-colnames(dmat)<-colnames(mat)
    while(dim(dmat)[1]>0){
      eqmat<-rbind(1:dim(dmat)[1],mapply(function(ub,b){if(isTRUE(all.equal(ub,b))==TRUE){1}else{0}},rep(dmat[1,"b"],dim(dmat)[1]),dmat[,"b"])); commat<-dmat[eqmat[1,eqmat[2,]==1],,drop=FALSE]; dmat<-dmat[setdiff(1:dim(dmat)[1],eqmat[1,eqmat[2,]==1]),,drop=FALSE]
      a1<-sum(commat[,"a1"]);w1<-sum(commat[,"w1"]);t1<-sum(commat[,"t1"]);wt1<-sum(commat[,"wt1"]);l1<-sum(commat[,"l1"]);a2<-sum(commat[,"a2"]);w2<-sum(commat[,"w2"]);t2<-sum(commat[,"t2"]);wt2<-sum(commat[,"wt2"]);l2<-sum(commat[,"l2"]);a3<-sum(commat[,"a3"]);w3<-sum(commat[,"w3"]);t3<-sum(commat[,"t3"]);wt3<-sum(commat[,"wt3"]);l3<-sum(commat[,"l3"])
      wt<-wt1+wt2+wt3; lip<-max(l1,l2,l3)
      nmat<-rbind(nmat,c(commat[1,"b"],commat[1,"mean"],commat[1,"sd"],a1,w1,t1,wt1,l1,a2,w2,t2,wt2,l2,a3,w3,t3,wt3,l3,wt,lip))}
    list(nmat=nmat)}
  brutenorm 	<- function(mat,bbm,bbsd){ # Calculate normalising costant, Lipshitz constant and normalized Lipshitz constant
    pbn<-brutepbn; sm<-mpfr(mat[,"sign"],pbn);am<-mpfr(mat[,"a"],pbn);bm<-mpfr(mat[,"b"],pbn);Lm<-mpfr(mat[,"L"],pbn);Um<-mpfr(mat[,"U"],pbn);trm<-mpfr(mat[,"trunc"],pbn);lipm<-mpfr(mat[,"lip"],pbn);bbmm<-mpfr(bbm,pbn);bbsdm<-mpfr(bbsd,pbn)
    wei<-sm*exp(am+bm*bbmm+0.5*(bm*bbsdm)^2); mat[,"wei"]<-as.numeric(wei); trwei<-wei*trm; mat[,"trwei"]<-as.numeric(trwei); wlip<-abs(wei)*lipm; mat[,"wlip"]<-as.numeric(wlip); Z<-sum(trwei); glip<-sum(wlip);
    nwei<-wei/Z; mat[,"nwei"]<-as.numeric(nwei); ntrwei<-nwei*trm; mat[,"ntrwei"]<-as.numeric(ntrwei); nwlip<-abs(nwei)*lipm; mat[,"nwlip"]<-as.numeric(nwlip); nglip<-sum(nwlip)
    list(mat=mat,Z=as.numeric(Z),glip=as.numeric(glip),nglip=as.numeric(nglip))}

  bruterhopi <- function(n, s, q, t, x, y, Ll, Lu, Ul, Uu) { # Calculate full density matrix, normalising constant and Lipschitz constant for given parameters
    # Required fields
    mat <- bruterho(n,s,q,t,x,y,Ll,Lu,Ul,Uu)
    rhom <- mat$rhom
    gamm1 <- mat$gamm1
    gamm2 <- mat$gamm2
    gamm3 <- mat$gamm3
    gamm4 <- mat$gamm4 # Rho matrix
    bbm	<- x + (q-s)*(y-x)/(t-s)
    bbsd <- ((t-q)*(q-s)/(t-s))^(0.5) # Mean and Variance
    # Matrix
    prmat <- matrix(0, (dim(rhom)[1]), 11)
    colnames(prmat) <- c("sign", "a", "b", "L","U", "weight", "trunc", "trweight", "mean", "sd", "lip")
    prmat[,"sign"] <- rhom[,"sign"]
    prmat[,"a"] <- rhom[,"a"]
    prmat[,"b"] <- rhom[,"b"]
    prmat[,"L"] <- rhom[,"L"]
    prmat[,"U"] <- rhom[,"U"]
    prmat[,"sd"] <- bbsd
    prmat[,"mean"] <- bbm+prmat[,"b"]*bbsd^2
    prmat[,"weight"] <- prmat[,"sign"]*exp(prmat[,"a"]+prmat[,"b"]*bbm+0.5*(prmat[,"b"]*bbsd)^2)
    prmat[,"trunc"] <- pnorm(prmat[,"U"], mean = prmat[,"mean"], sd = bbsd)-pnorm(prmat[,"L"], mean = prmat[,"mean"], sd = bbsd)
    prmat[,"trweight"] <- prmat[,"weight"]*prmat[,"trunc"]
    for(i in 1:dim(prmat)[1]) {
      prmat[i,"lip"] <- brutelip(prmat[i,"mean"], bbsd, prmat[i,"L"], prmat[i,"U"])
    }
    m1 <- dim(gamm1)[1]
    m2 <- m1 + dim(gamm2)[1]
    m3 <- m2 + dim(gamm3)[1]
    m4 <- m3 + dim(gamm4)[1]
    prmat1 <- prmat[1:m1,]
    prmat2 <- prmat[(m1+1):m2,]
    prmat3 <- prmat[(m2+1):m3,]
    prmat4 <- prmat[(m3+1):m4,]
    # Band split matrix
    splmat1 <- splmat2 <- splmat3 <- matrix(0, 0, 11)
    colnames(splmat1) <- colnames(splmat2) <- colnames(splmat3) <- colnames(prmat)
    splmat1 <- rbind(splmat1, prmat1, prmat3)
    splmat1[,"L"] <- Ll
    splmat1[,"U"] <- Lu
    nmat1 <- brutenett1(splmat1, bbm)$nmat
    splmat2 <- rbind(splmat2, prmat1, prmat2, prmat3, prmat4)
    splmat2[,"L"] <- Lu
    splmat2[,"U"] <- Ul
    nmat2 <- brutenett1(splmat2, bbm)$nmat
    splmat3 <- rbind(splmat3, prmat1, prmat2)
    splmat3[,"L"] <- Ul
    splmat3[,"U"] <- Uu
    nmat3 <- brutenett1(splmat3, bbm)$nmat
    # Net matrix
    nmata <- rbind(nmat1, nmat2, nmat3)
    nmatb <- matrix(0, 0, 20)
    colnames(nmatb) <- c("b", "mean", "sd", "a1", "w1", "t1", "wt1", "l1", "a2", "w2", "t2", "wt2", "l2", "a3", "w3", "t3", "wt3", "l3", "wt", "l")
    nmatb <- rbind(nmatb, cbind(nmat1[,c(3,9,10),drop=FALSE],
                                nmat1[,c(2,6,7,8,11),drop=FALSE],
                                rep(0,dim(nmat1)[1]),
                                rep(0,dim(nmat1)[1]),
                                rep(0,dim(nmat1)[1]),
                                rep(0,dim(nmat1)[1]),
                                rep(0,dim(nmat1)[1]),
                                rep(0,dim(nmat1)[1]),
                                rep(0,dim(nmat1)[1]),
                                rep(0,dim(nmat1)[1]),
                                rep(0,dim(nmat1)[1]),
                                rep(0,dim(nmat1)[1]),
                                rep(0,dim(nmat1)[1]),
                                rep(0,dim(nmat1)[1])))
    nmatb <- rbind(nmatb, cbind(nmat2[,c(3,9,10),drop=FALSE],
                                rep(0,dim(nmat2)[1]),
                                rep(0,dim(nmat2)[1]),
                                rep(0,dim(nmat2)[1]),
                                rep(0,dim(nmat2)[1]),
                                rep(0,dim(nmat2)[1]),
                                nmat2[,c(2,6,7,8,11),drop=FALSE],
                                rep(0,dim(nmat2)[1]),
                                rep(0,dim(nmat2)[1]),
                                rep(0,dim(nmat2)[1]),
                                rep(0,dim(nmat2)[1]),
                                rep(0,dim(nmat2)[1]),
                                rep(0,dim(nmat2)[1]),
                                rep(0,dim(nmat2)[1])))
    nmatb <- rbind(nmatb, cbind(nmat3[,c(3,9,10),drop=FALSE],
                                rep(0,dim(nmat3)[1]),
                                rep(0,dim(nmat3)[1]),
                                rep(0,dim(nmat3)[1]),
                                rep(0,dim(nmat3)[1]),
                                rep(0,dim(nmat3)[1]),
                                rep(0,dim(nmat3)[1]),
                                rep(0,dim(nmat3)[1]),
                                rep(0,dim(nmat3)[1]),
                                rep(0,dim(nmat3)[1]),
                                rep(0,dim(nmat3)[1]),
                                nmat3[,c(2,6,7,8,11),drop=FALSE],
                                rep(0,dim(nmat3)[1]),
                                rep(0,dim(nmat3)[1])))
    nmatb <- brutenett2(nmatb)
    # Normalised matrix
    normm <- cbind(nmata[,c(1:6,8,6,8,7,9,10,11,11,11)])
    colnames(normm) <- c("sign", "a", "b", "L", "U", "wei", "trwei", "nwei", "ntrwei", "trunc", "mean", "sd", "lip", "wlip", "nwlip")
    normmfn <- brutenorm(normm,bbm,bbsd)
    normm <- normmfn$mat
    Z <- normmfn$Z
    glip <- normmfn$glip
    nglip <- normmfn$nglip
    # Output
    list(bbm=bbm,bbsd=bbsd,rhom=rhom,prmat=prmat,nmata=nmata,nmatb=nmatb,normm=normm,Z=Z,gausslip=abs((1/sqrt(2*pi*bbsd^2))*(-1/bbsd)*exp(-1/2)),glip=glip,nglip=nglip,Ll=Ll,Lu=Lu,Ul=Ul,Uu=Uu)
  }

  bruteneval <- function(m, s, q, t, x, y, Ll, Lu, Ul, Uu) { # m odd only # Evaluate ratio of alternating Cauchy sequence normalising constants
    upper <- bruterhopi(m,s,q,t,x,y,Ll,Lu,Ul,Uu)
    lower <- bruterhopi(m+1,s,q,t,x,y,Ll,Lu,Ul,Uu)
    ratio <- min(lower$Z/upper$Z,upper$Z/lower$Z)
    list(upper=upper,lower=lower,ratio=ratio)
  }

  earplays 	<- function(mwt,mmu){ # Nett matrix type 1
    nmat <- matrix(0, 0, 2)
    dmat <- cbind(mmu, mwt)
    colnames(dmat) <- c("mean", "weight")
    colnames(nmat) <- c("weight", "mean")
    while(dim(dmat)[1] > 0) {
      eqmat <- rbind(1:dim(dmat)[1],
                     mapply(function(ub, b) { if(isTRUE(all.equal(ub, b)) == TRUE) { 1 } else { 0 } }, rep(dmat[1,"mean"], dim(dmat)[1]), dmat[,"mean"]))
      commat <- dmat[eqmat[1, eqmat[2,] == 1],, drop = FALSE]
      dmat <- dmat[setdiff(1:dim(dmat)[1], eqmat[1,eqmat[2,] == 1]),, drop = FALSE]
      wnew <- sum(commat[,"weight"])
      if(wnew != 0) {
        nmat <- rbind(nmat, c(wnew, commat[1,"mean"]))
      }
    }
    nmat
  }

  earpbrute <- function(s, q, t, x, y, Ll, Lu, Ul, Uu) { # Bound piecewise constant uniforms
    m <- -1
    R1 <- 0
    while(R1 <= 0.5) {
      m <- m+4
      altmat <- bruteneval(m, s, q, t, x, y, Ll, Lu, Ul, Uu)
      upper <- altmat$upper
      R1 <- altmat$ratio
    } # R1 denotes the ratio of the lower to upper normalising constant
    mesh <- 0
    R2 <- 0
    while(R2 <= 0.5) {
      mesh <- mesh+1
      bndmat <- brutemesh(upper, mesh)
      pmat <- bndmat$pmat
      R2 <- bndmat$aprob
      R3 <- R1*R2
    } # R2 denotes the ratio of the upper to bounding normalising constant, R3 is the lower bound on the acceptance probability
    dind1I <- 0
    while(dind1I == 0) {
      band <- pmat[sample(1:dim(pmat)[1], 1, prob = pmat[,"prob"]),, drop = FALSE]
      draw <- runif(1, band[,"s"], band[,"t"])
      u <- runif(1, 0, band[,"bnd"])
      nev <- 30+m
      dind2I <- 0
      while(dind2I == 0) {
        cauc <- earh3C(nev, s, q, t, x, draw, y, Ll, Lu, Ul, Uu)*dnorm(draw, mean = upper$bbm, sd = upper$bbsd)
        if(u <= cauc[2]) {
          dind1I <- dind2I <- 1
        } else {
          if(u >= cauc[3]) {
            dind2I <- 1
          } else {
            nev <- nev+20
          }
        }
      }
    }

    list(draw = draw, altmat = altmat, upper = upper, bndmat = bndmat, pmat = pmat, R1 = R1, R2 = R2, R3 = R3)
  }

  cnt1I <- 0
  cnt1M <- 250
  cnt2M <- 25 # earpreg breakout count

  bbm <- x + (q-s)*(y-x)/(t-s)
  bbsd <- ((t-q)*(q-s)/(t-s))^(0.5)
  P1 <- -2/(q-s)
  P2 <- -2/(t-q)
  Alu <- Uu-Ll
  Auu <- Uu-Lu
  All <- Ul-Ll
  Aul <- Ul-Lu
  mat <- matrix(c(-1,easiga(P1,x,Alu,Ll),easigb(P1,x,Alu,Ll),
                  -1,easiga(P1,-x,Alu,-Uu),-easigb(P1,-x,Alu,-Uu),
                  1,eaphia(P1,x,Alu,-1),eaphib(P1,x,Alu,-1),
                  1,eaphia(P1,-x,Alu,-1),-eaphib(P1,-x,Alu,-1),
                  -1,easiga(P2,y,Alu,Ll),easigb(P2,y,Alu,Ll),
                  -1,easiga(P2,-y,Alu,-Uu),-easigb(P2,-y,Alu,-Uu),
                  1,eaphia(P2,y,Alu,1),eaphib(P2,y,Alu,1),
                  1,eaphia(P2,-y,Alu,1),-eaphib(P2,-y,Alu,1),
                  1,easiga(P1,x,Alu,Ll)+easiga(P2,y,Alu,Ll),easigb(P1,x,Alu,Ll)+easigb(P2,y,Alu,Ll),
                  1,easiga(P1,x,Alu,Ll)+easiga(P2,-y,Alu,-Uu),easigb(P1,x,Alu,Ll)-easigb(P2,-y,Alu,-Uu),
                  1,easiga(P1,-x,Alu,-Uu)+easiga(P2,y,Alu,Ll),-easigb(P1,-x,Alu,-Uu)+easigb(P2,y,Alu,Ll),
                  1,easiga(P1,-x,Alu,-Uu)+easiga(P2,-y,Alu,-Uu),-easigb(P1,-x,Alu,-Uu)-easigb(P2,-y,Alu,-Uu),
                  1,easiga(P1,x,Auu,Lu),easigb(P1,x,Auu,Lu),
                  1,easiga(P1,-x,Auu,-Uu),-easigb(P1,-x,Auu,-Uu),
                  1,easiga(P2,y,Auu,Lu),easigb(P2,y,Auu,Lu),
                  1,easiga(P2,-y,Auu,-Uu),-easigb(P2,-y,Auu,-Uu),
                  -1,easiga(P1,x,Auu,Lu)+easiga(P2,y,Auu,Lu),easigb(P1,x,Auu,Lu)+easigb(P2,y,Auu,Lu),
                  -1,easiga(P1,x,Auu,Lu)+easiga(P2,-y,Auu,-Uu),easigb(P1,x,Auu,Lu)-easigb(P2,-y,Auu,-Uu),
                  1,easiga(P1,x,Auu,Lu)+eaphia(P2,y,Auu,1),easigb(P1,x,Auu,Lu)+eaphib(P2,y,Auu,1),
                  1,easiga(P1,x,Auu,Lu)+eaphia(P2,-y,Auu,1),easigb(P1,x,Auu,Lu)-eaphib(P2,-y,Auu,1),
                  -1,easiga(P1,-x,Auu,-Uu)+easiga(P2,y,Auu,Lu),-easigb(P1,-x,Auu,-Uu)+easigb(P2,y,Auu,Lu),
                  -1,easiga(P1,-x,Auu,-Uu)+easiga(P2,-y,Auu,-Uu),-easigb(P1,-x,Auu,-Uu)-easigb(P2,-y,Auu,-Uu),
                  1,easiga(P1,-x,Auu,-Uu)+eaphia(P2,y,Auu,1),-easigb(P1,-x,Auu,-Uu)+eaphib(P2,y,Auu,1),
                  1,easiga(P1,-x,Auu,-Uu)+eaphia(P2,-y,Auu,1),-easigb(P1,-x,Auu,-Uu)-eaphib(P2,-y,Auu,1),
                  1,eaphia(P1,x,Auu,-1)+easiga(P2,y,Auu,Lu),eaphib(P1,x,Auu,-1)+easigb(P2,y,Auu,Lu),
                  1,eaphia(P1,x,Auu,-1)+easiga(P2,-y,Auu,-Uu),eaphib(P1,x,Auu,-1)-easigb(P2,-y,Auu,-Uu),
                  -1,eaphia(P1,x,Auu,-1)+eaphia(P2,y,Auu,1),eaphib(P1,x,Auu,-1)+eaphib(P2,y,Auu,1),
                  -1,eaphia(P1,x,Auu,-1)+eaphia(P2,-y,Auu,1),eaphib(P1,x,Auu,-1)-eaphib(P2,-y,Auu,1),
                  1,eaphia(P1,-x,Auu,-1)+easiga(P2,y,Auu,Lu),-eaphib(P1,-x,Auu,-1)+easigb(P2,y,Auu,Lu),
                  1,eaphia(P1,-x,Auu,-1)+easiga(P2,-y,Auu,-Uu),-eaphib(P1,-x,Auu,-1)-easigb(P2,-y,Auu,-Uu),
                  -1,eaphia(P1,-x,Auu,-1)+eaphia(P2,y,Auu,1),-eaphib(P1,-x,Auu,-1)+eaphib(P2,y,Auu,1),
                  -1,eaphia(P1,-x,Auu,-1)+eaphia(P2,-y,Auu,1),-eaphib(P1,-x,Auu,-1)-eaphib(P2,-y,Auu,1),
                  1,easiga(P1,x,All,Ll),easigb(P1,x,All,Ll),
                  1,easiga(P1,-x,All,-Ul),-easigb(P1,-x,All,-Ul),
                  1,easiga(P2,y,All,Ll),easigb(P2,y,All,Ll),
                  1,easiga(P2,-y,All,-Ul),-easigb(P2,-y,All,-Ul),
                  -1,easiga(P1,x,All,Ll)+easiga(P2,y,All,Ll),easigb(P1,x,All,Ll)+easigb(P2,y,All,Ll),
                  -1,easiga(P1,x,All,Ll)+easiga(P2,-y,All,-Ul),easigb(P1,x,All,Ll)-easigb(P2,-y,All,-Ul),
                  1,easiga(P1,x,All,Ll)+eaphia(P2,y,All,1),easigb(P1,x,All,Ll)+eaphib(P2,y,All,1),
                  1,easiga(P1,x,All,Ll)+eaphia(P2,-y,All,1),easigb(P1,x,All,Ll)-eaphib(P2,-y,All,1),
                  -1,easiga(P1,-x,All,-Ul)+easiga(P2,y,All,Ll),-easigb(P1,-x,All,-Ul)+easigb(P2,y,All,Ll),
                  -1,easiga(P1,-x,All,-Ul)+easiga(P2,-y,All,-Ul),-easigb(P1,-x,All,-Ul)-easigb(P2,-y,All,-Ul),
                  1,easiga(P1,-x,All,-Ul)+eaphia(P2,y,All,1),-easigb(P1,-x,All,-Ul)+eaphib(P2,y,All,1),
                  1,easiga(P1,-x,All,-Ul)+eaphia(P2,-y,All,1),-easigb(P1,-x,All,-Ul)-eaphib(P2,-y,All,1),
                  1,eaphia(P1,x,All,-1)+easiga(P2,y,All,Ll),eaphib(P1,x,All,-1)+easigb(P2,y,All,Ll),
                  1,eaphia(P1,x,All,-1)+easiga(P2,-y,All,-Ul),eaphib(P1,x,All,-1)-easigb(P2,-y,All,-Ul),
                  -1,eaphia(P1,x,All,-1)+eaphia(P2,y,All,1),eaphib(P1,x,All,-1)+eaphib(P2,y,All,1),
                  -1,eaphia(P1,x,All,-1)+eaphia(P2,-y,All,1),eaphib(P1,x,All,-1)-eaphib(P2,-y,All,1),
                  1,eaphia(P1,-x,All,-1)+easiga(P2,y,All,Ll),-eaphib(P1,-x,All,-1)+easigb(P2,y,All,Ll),
                  1,eaphia(P1,-x,All,-1)+easiga(P2,-y,All,-Ul),-eaphib(P1,-x,All,-1)-easigb(P2,-y,All,-Ul),
                  -1,eaphia(P1,-x,All,-1)+eaphia(P2,y,All,1),-eaphib(P1,-x,All,-1)+eaphib(P2,y,All,1),
                  -1,eaphia(P1,-x,All,-1)+eaphia(P2,-y,All,1),-eaphib(P1,-x,All,-1)-eaphib(P2,-y,All,1),
                  -1,easiga(P1,x,Aul,Lu),easigb(P1,x,Aul,Lu),
                  -1,easiga(P1,-x,Aul,-Ul),-easigb(P1,-x,Aul,-Ul),
                  1,eaphia(P1,x,Aul,-1),eaphib(P1,x,Aul,-1),
                  1,eaphia(P1,-x,Aul,-1),-eaphib(P1,-x,Aul,-1),
                  -1,easiga(P2,y,Aul,Lu),easigb(P2,y,Aul,Lu),
                  -1,easiga(P2,-y,Aul,-Ul),-easigb(P2,-y,Aul,-Ul),
                  1,eaphia(P2,y,Aul,1),eaphib(P2,y,Aul,1),
                  1,eaphia(P2,-y,Aul,1),-eaphib(P2,-y,Aul,1),
                  1,easiga(P1,x,Aul,Lu)+easiga(P2,y,Aul,Lu),easigb(P1,x,Aul,Lu)+easigb(P2,y,Aul,Lu),
                  1,easiga(P1,x,Aul,Lu)+easiga(P2,-y,Aul,-Ul),easigb(P1,x,Aul,Lu)-easigb(P2,-y,Aul,-Ul),
                  1,easiga(P1,-x,Aul,-Ul)+easiga(P2,y,Aul,Lu),-easigb(P1,-x,Aul,-Ul)+easigb(P2,y,Aul,Lu),
                  1,easiga(P1,-x,Aul,-Ul)+easiga(P2,-y,Aul,-Ul),-easigb(P1,-x,Aul,-Ul)-easigb(P2,-y,Aul,-Ul)),
                nrow = 64, ncol = 3, byrow = TRUE)
  mwt <- mat[,1]*exp(0.5*(bbsd*mat[,3])^2+bbm*mat[,3]+mat[,2])
  mwt <- mwt[c(1:12,33:52,1:64,1:32)]
  mmu <- bbm+(bbsd^2)*mat[,3]
  mmu<-mmu[c(1:12,33:52,1:64,1:32)]
  umat1 <- earplays(mwt[1:32], mmu[1:32])
  umat2<-earplays(mwt[33:96], mmu[33:96])
  umat3<-earplays(mwt[97:128], mmu[97:128])
  nmat <- rbind(cbind(L = Ll, U = Lu, B = 1, umat1), cbind(L = Lu, U = Ul, B = 2, umat2), cbind(L = Ul, U = Uu, B = 3, umat3))
  nmat <- nmat[nmat[,4] != 0,]
  eb <- (nmat[,5]-bbm)/(bbsd^2)
  ea <- log(abs(nmat[,4]))-0.5*(bbsd*eb)^2-bbm*eb
  esi <- sign(nmat[,4])
  Lp <- pnorm(nmat[,1], nmat[,5], sd = bbsd)
  Up <- pnorm(nmat[,2], nmat[,5], sd = bbsd)
  itwt <- nmat[,4]*(Up-Lp)
  # nmat <- cbind(esi, ea, eb, wt = nmat[,4], itwt = (itwt)/sum(itwt), L = nmat[,1], Lp, U = nmat[,2], Up, mmu = nmat[,5], B = nmat[,3])
  # Debugging 10/3/2021
  nmat <- cbind(esi, ea, eb, wt = nmat[,4], itwt = itwt, L = nmat[,1], Lp, U = nmat[,2], Up, mmu = nmat[,5], B = nmat[,3])
  pmat <- nmat[nmat[,1] == 1,]
  dind1 <- 0
  while(dind1 == 0) {
    dind2 <- 0
    cnt2I <- 0
    while(dind2 == 0) {
      cnt2I <- cnt2I+1; if(any(pmat[,5]<0)) browser()
      sp <- pmat[sample(1:(dim(pmat)[1]), 1, replace = TRUE, prob = pmat[,5]),]
      dr <- qnorm(runif(1, sp[7], sp[9]), sp[10], bbsd)
      pmt <- pmat[pmat[,11] == sp[11],]
      nmt <- nmat[nmat[,11] == sp[11],]
      u <- runif(1, 0, sum(pmt[,1]*exp(pmt[,2]+pmt[,3]*dr)))
      if(u <= sum(nmt[,1]*exp(nmt[,2]+nmt[,3]*dr))) {
        dind2 <- 1
      }
      if(dind2 == 0) {
        if(cnt2I >= cnt2M) {
          dind2 <- -1
        }
      }
    }
    if(dind2 == 1) {
      dind3 <- 0
      m <- 3
      while(dind3 == 0) {
        cauc <- earh3C(m, s, q, t, x, dr, y, Ll, Lu, Ul, Uu)
        if(u <= cauc[2]) {
          dind1 <- dind3 <- 1
        } else {
          if(u >= cauc[3]) {
            dind3 <- 1
          } else {
            m <- m+2
          }
        }
      }
    }
    if(dind1 == 0) {
      cnt1I <- cnt1I+1
      if(cnt1I >= cnt1M) {
        dind1 <- 1
        dr <- earpbrute(s, q, t, x, y, Ll, Lu, Ul, Uu)$draw
      }
    }
  } # Index breakout counter and breakout with draw
  dr
}



sim.condintersection.simlyr_ <- function(s, q, t, x, w, y, Ll, Lu, Ul, Uu) {
  Lmin <- min(x,w)
  Lmax <- max(x,w)
  Rmin <- min(w,y)
  Rmax <- max(w,y)
  if(w >= Ul) {
    Ul<-w
  }
  if(w <= Lu) {
    Lu <- w
  }
  Lmat <- matrix(c(Ll,Lu,Ul,Uu,
                   Ll,Lu,Ul,Uu,
                   Lu,Lmin,Ul,Uu,
                   Ll,Lu,Ul,Uu,
                   Ll,Lu,Ul,Uu,
                   Lu,Lmin,Ul,Uu,
                   Ll,Lu,Lmax,Ul,
                   Ll,Lu,Lmax,Ul,
                   Lu,Lmin,Lmax,Ul),
                 nrow = 9, ncol = 4, byrow = TRUE)
  Rmat <- matrix(c(Ll,Lu,Ul,Uu,
                   Lu,Rmin,Ul,Uu,
                   Ll,Lu,Ul,Uu,
                   Ll,Lu,Rmax,Ul,
                   Lu,Rmin,Rmax,Ul,
                   Ll,Lu,Rmax,Ul,
                   Ll,Lu,Ul,Uu,
                   Lu,Rmin,Ul,Uu,
                   Ll,Lu,Ul,Uu),
                 nrow = 9, ncol = 4, byrow = TRUE)
  LRmat <- cbind(Lmat, Rmat)
  bbrind <- deind <- 0
  bbrct <- 1
  m1 <- 3
  u1 <- runif(1, 0, 1)
  while(bbrind == 0) {
    while(deind==0) {
      debd <- earhoC(m1,s,q,t,x,w,y,Ll,Lu,Ul,Uu)
      if(debd[2] <= 0) {
        m1 <- m1+2
      } else {
        deind <- 1
      }
    }
    bebe <- matrix(0, bbrct, 4)
    for(i in 1:bbrct) {
      bebe[i,1:2] <- eabetaC(m1,s,q,x,w,LRmat[i,1],LRmat[i,2],LRmat[i,3],LRmat[i,4])
      bebe[i,3:4] <- eabetaC(m1,q,t,w,y,LRmat[i,5],LRmat[i,6],LRmat[i,7],LRmat[i,8])
    }
    bd <- c(sum(bebe[,1]*bebe[,3])/debd[1], sum(bebe[,2]*bebe[,4])/debd[2])
    if(u1 <= bd[1]) {
      bbrind <- 1
    } else {
      if(u1 >= bd[2]) {
        bbrct <- bbrct+1
        deind <- 0
      } else {
        m1 <- m1+2
        deind <- 0
      }
    }
  }
  matrix(c(s,q,x,w,LRmat[bbrct,1:4],
           q,t,w,y,LRmat[bbrct,5:8]),
         nrow = 2, ncol = 8, byrow = TRUE)
}
