# sumFREGAT (2017) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

sumstat.famSKAT <- function(Z) {
	Q <- sum((Z$w * Z$Z) ^ 2) # statistic
	KKK <- t(Z$U * Z$w) * Z$w # kernel matrix
	if (method != 'hybrid') {
	eig <- eigen(KKK, symmetric = TRUE, only.values = TRUE)
	ev <- eig$values[eig$values > 1e-6 * eig$values[1]]
	if (method == 'kuonen') {
		p <- pchisqsum(Q, rep(1, length(ev)), ev, lower.tail = F, method = 'sad') 
	} else if (method == 'davies') {
		p <- davies(Q, ev, acc = acc, lim = lim)$Qq
	} 
	} else if (method == 'hybrid') {
        lam = svd(KKK, nu=0,nv=0)$d; 	#lam <- lam[lam > 1e-6 * lam[1]]
		p <- KAT.pval(Q, lam)
	}
	return(max(min(p, 1), 0))
}

integrateNEW = function(T0,katint,q1,Pmin,m1) {
 p.value = try({ T0 + integrate(katint,0,q1,subdivisions=1000,abs.tol=1e-25)$val }, silent=TRUE)
 return(p.value)
}

sumstat.famSKATO <- function(Z) {

	Q05 <- Z$w * Z$Z  # weighting of Z-score statistic
	KKK <- t(Z$U * Z$w) * Z$w  # kernel matrix
	m1 <- length(Z$Z)

if (method != 'hybrid')	{
	eig <- eigen(KKK, symmetric = TRUE)
	eig$values[eig$values <= 0] <- 1e-7

	C05 <- eig$vec %*% (t(eig$vec) * sqrt(eig$val))
	Q.all <- c()
	for (rh in rhos) {
		CORB <- matrix(NA, m1, m1)
		if (rh < 1) {
			CORB05 <- chol(diag(m1) * (1 - rh) + rh)  # (correlation matrix for betas) ^ 0.5
		} else { CORB05 <- matrix(1, m1, m1) / sqrt(m1) }  # case of famBT (rh = 1)
		Q05W <- Q05 %*% t(CORB05)
		Q <- sum(Q05W * Q05W)
		Q.all <- c(Q.all, Q)
	}
	Q <- rbind(Q.all, NULL)
	out <- SKAT_Optimal_Get_Pvalue(Q, C05, rhos, method, acc, lim)
	return(out$p.value[1])
	
} else {

  K = length(rhos); K1 = K
  Qs = sum(Q05^2) ;   Qb = sum(Q05)^2 ;   Qw = (1-rhos)*Qs + rhos*Qb
  pval = rep(0,K)
  Rs = rowSums(KKK); R1 = sum(Rs); R2 = sum(Rs^2); R3 = sum(Rs*colSums(KKK * Rs))
  RJ2 = outer(Rs,Rs,'+')/m1

  if(rhos[K]>=1){ K1 = K-1;     pval[K] = pchisq(Qb/R1, 1, lower.tail=FALSE)   }
  Lamk = vector('list', K1);  rho1 = rhos[1:K1]; 
  t1 = sqrt(1-rho1)
  tmp = sqrt(1-rho1 + m1*rho1) - t1
  c1 = t1*tmp;  c2 = ((tmp/m1)^2)*R1

  for(k in 1:K1){
    mk = (1-rhos[k])*KKK + c1[k]*RJ2 + c2[k]
	SSS <- eigen(mk,symmetric=TRUE,only.values = TRUE)$val; SSS <- SSS[SSS > 1e-08]	#SSS1 <- svd(mk, nu=0,nv=0)$d
	Lamk[[k]] = pmax(SSS, 0)
    pval[k] = KAT.pval(Qw[k],Lamk[[k]])
  }
  Pmin = min(pval)
  qval = rep(0,K1)
  for(k in 1:K1) qval[k] = Liu.qval.mod(Pmin, Lamk[[k]])
    
  SVD<-eigen(KKK - outer(Rs,Rs)/R1, symmetric=TRUE, only.values = TRUE)$val #SVD1<-svd(KKK - outer(Rs,Rs)/R1, nu=0,nv=0)$d

  lam = pmax(SVD[-m1], 0) 	
  tauk = (1-rho1)*R2/R1 + rho1*R1;  vp2 = 4*(R3/R1-R2^2/R1^2)
  MuQ = sum(lam);  VarQ = sum(lam^2)*2
  sd1 = sqrt(VarQ/(VarQ+vp2))
  MMM <- MuQ - MuQ*sd1
  
  if(K1<K){
    q1 = qchisq(Pmin,1,lower.tail=FALSE)
    T0 = Pmin
  } else{
    tmp = ( qval-(1-rhos)*MMM/sd1 )/tauk
    q1 = min(tmp)
    T0 = pchisq(q1,1,lower.tail=FALSE)
  }
  
#-------------------------
  katint <- function(xpar){
    eta1 = sapply(xpar, function(eta0) min(X2  - X3 *  eta0 ))
    return(KAT.pval(eta1,lam)*dchisq(xpar,1))
  }
#-----------------------
 
X1<-(1-rho1); X2 <-  sd1 * qval/X1 + MMM; X3 <- sd1 * tauk/X1

p.value <- integrateNEW(T0,katint,q1,Pmin,m1)

min(p.value, Pmin*K)
}
}

