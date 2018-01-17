# sumFREGAT (2017) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

sumstat.famBT <- function(Z) {
	chi2 <- sum(Z$w * Z$Z) ^ 2
	KKK <- sum(t(Z$U * Z$w) * Z$w) 
	chi2 <- chi2/KKK
	p <- pchisq(chi2, 1, lower.tail = FALSE)
	betaBT <- sum(Z$w * Z$Z)/KKK * length(Z$w)
	se.beta <- sqrt((betaBT ^ 2) / chi2)
	return(c(p, betaBT, se.beta))

}

sumstat.MLR <- function(Z) {
	U_1 <- Z$U %^% (-1)
	m <- qr(Z$U)$rank  # m <- length(Zstat)
	RSS <- n - sum(Z$Z * (U_1 %*% Z$Z))
	Fstat <- ((n - m) / m) * (n - RSS) / RSS  # F-statistic
	p <- pf(Fstat, m, n - m, lower.tail = FALSE)
	return(p)
}

sumstat.famFLM <- function(Z) {

	m1 <- length(Z$Z)
	if (m1 <= kb0) {
		return(c(sumstat.MLR(Z), 'MLR'))
	}

	B <- eval.basis(Z$pos, betabasis0)  # as matrix (m x Kb)
	WB <- B*Z$w
	Mat <- t(WB) %*% Z$U %*% WB
	BUB_1 <- Mat %^% (-1)
	m <- qr(Mat)$rank  # Kb <- length(BZstat)

	BZstat <- as.vector(t(WB) %*% Z$Z)
	RSS <- n - sum(BZstat * (BUB_1 %*% BZstat))

	Fstat <- ((n - m) / m) * (n - RSS) / RSS   # F-statistic
	p <- pf(Fstat, m, n - m, lower.tail = FALSE)
	return(c(p, model0))
}

"%^%" <- function(U, k) {
	UUU <- eigen(U, symmetric = TRUE)  # UUU = Uvec %*% diag(Uval) %*% t(Uvec)
	Uval <- UUU$val
	Uvec <- UUU$vec
	Uvec <- Uvec[, Uval > 1e-7]
	Uval <- Uval[Uval > 1e-7]
	Uvec %*% (t(Uvec) * (Uval ^ k))   #	Uvec %*% (diag(Uval ^ k) %*% t(Uvec))

}

sumstat.PCA <- function(Z) {
	U_05 <- Z$U %^% (-.5)
	yc <- as.vector(U_05 %*% Z$Z)
	Gc <- t((Z$w * sqrt(n + 1)) * t(Z$U %^% .5))
	numberPCA(yc, Gc, n, var.fraction)
}

numberPCA <- function(yc, Gc, n, var.fraction) {
	m <- qr(Gc)$rank
	pCA <- PC(Gc, n) # n = N - 1
	COMPs0 <- pCA$scores[, 1:m]     # matrix of independent scores
	CPV <- pCA$importance[3, 1:m]   # Cumulative Proportion of Variance (CPV)
	comp <- which(CPV >= var.fraction)      # components for which Explained variance fraction is about 85%
	M <- min(comp)
	COMPs <- as.matrix(COMPs0[,1:M])
	m <- qr(COMPs)$rank
	GY <- as.vector(t(COMPs) %*% yc)
	CC <- t(COMPs) %*% COMPs
	if (M > 1) {
		RSS <- n - sum(GY * as.vector((CC %^% (-1)) %*% GY))
	} else { RSS <- n - GY * GY / CC }
	Fstat <- ((n - m) / m) * (n - RSS) / RSS    # F-statistic
	p <- pf(Fstat, m, n - m, lower.tail = FALSE)
	c(p, M, CPV[M])
}

PC <- function(X, n) {
	U <- (t(X) %*% X) / n # n = N - 1
	eX1 <- eigen(U, symmetric = TRUE)
	values <- eX1$values
	vectors <- eX1$vectors
	values[values < 0] <- 0
	sdev <- sqrt(values)
	ind.var <- sdev ^ 2
	total.var <- sum(ind.var)
	scores <- X %*% vectors
	scorenames <- paste('PC', 1:ncol(X), sep = '')
	colnames(scores) <- colnames(vectors) <- scorenames
	rownames(vectors) <- colnames(X)
	prop.var <- ind.var / total.var
	cum.var <- rep(NA, ncol(X))
	for (i in 1:ncol(X)) { cum.var[i] <- sum(prop.var[1:i]) }
	importance <- rbind(sdev, prop.var, cum.var)
	colnames(importance) <- scorenames
	rownames(importance) <- c("Standard Deviation", "Proportion of Variance", "Cumulative Proportion")
	list(values=values,vectors=vectors,scores=scores,importance=importance,sdev=sdev)
}
