# sumFREGAT (2017-2018) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

integrateNEW = function(T0,katint,q1,Pmin,m1) {
 p.value = try({ T0 + integrate(katint,0,q1,subdivisions=1000,abs.tol=1e-25)$val }, silent=TRUE)
 return(p.value) 
}
#--------------------------------------------------------
search.subregions <- function(x,threshold){
	m = length(x); A = matrix(0,m,m)
	x = ifelse(abs(x) < threshold,1,0)
	A[1,1]=1; K=1
	for (i in 2:m)
		if (x[i]==1 & x[i-1]==1) { A[i,K]=1 }
		else {K=K+1; A[i,K]=1}
	A[,1:K]
}

#-------------------------------------------------------

sumstat.SKATO <- function(obj) {

	with(obj, with(df, { # Z, U, w, method, acc, lim, rhos, p.threshold

		Q05 <- w * Z  # weighting of Z-score statistic
		KKK <- t(U * w) * w  # kernel matrix
		m1 <- length(Z)
	#if (m1 >= 50){  # begin for long region
	#    p.threshold <- 0.95
		if (p.threshold > 0) {
			threshold = qnorm (1. - (p.threshold/2))
			A <- as.matrix(search.subregions(Z, threshold))
			if (dim(A)[2] > 1) {
				Q05 <- t(A) %*% Q05
				KKK <- t(A) %*% KKK %*% A  # compressed kernel matrix
				m1 <- length(Q05)
			} # end
		}
		
		if (method != 'hybrid')	{
			eig <- eigen(KKK, symmetric = TRUE)
			eig$values[eig$values <= 0] <- 1e-7

			C05 <- t(eig$vec) * sqrt(eig$val)   #C05 <- eig$vec %*% (t(eig$vec) * sqrt(eig$val))
			Q.all <- c()
			
			aaa <- sum(Q05^2); bbb <-(sum(Q05))^2
			for (rh in rhos) {
				Q <- aaa*(1-rh) + bbb *rh
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

	}))

}


'SKATO' <- function (score.file, gene.file, genes = 'all', cor.path = 'cor/', anno.type = '', beta.par = c(1, 25), weights.function = NULL,
user.weights = FALSE, method = 'kuonen', acc = 1e-8, lim = 1e+6, rho = TRUE, p.threshold = 0.8, write.file = FALSE) {

	SKAT(score.file, gene.file, genes, cor.path, anno.type, beta.par, weights.function,
	user.weights, gen.var.weights = 'se.beta', method, acc, lim, rho,
	p.threshold, write.file)

}