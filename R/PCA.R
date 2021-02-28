# sumFREGAT (2017-2018) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

"%^%" <- function(U, k) {
	UUU <- eigen(U, symmetric = TRUE)  # UUU = Uvec %*% diag(Uval) %*% t(Uvec)
	Uval <- UUU$val
	Uvec <- UUU$vec
	Uvec <- Uvec[, Uval > 1e-7]
	Uval <- Uval[Uval > 1e-7]
	Uvec %*% (t(Uvec) * (Uval ^ k))   #	Uvec %*% (diag(Uval ^ k) %*% t(Uvec))

}

sumstat.PCA <- function(obj) {

	with(obj, with(df, { # Z, w, U, n, var.fraction

		WZ  <- as.vector(w * Z) * sqrt(n)
		WUW <- as.matrix(t(U * w) * w) * n
		numberPCA(WZ, WUW, n, var.fraction)

	}))

}

numberPCA <- function(WZ, WUW, n, var.fraction) {
	pCA <- PC(WUW) # n = N - 1
	CPV <- pCA$importance   # Cumulative Proportion of Variance (CPV)
	M <- min(which(CPV >= var.fraction))  # components for which Explained variance fraction is about 85%
	BBB <- as.matrix(pCA$scores[,1:M])
	GY <- as.vector(t(BBB) %*% WZ)
	CC <- as.matrix(t(BBB) %*% WUW %*% BBB)
	m <- qr(BBB)$rank
	if (m > 1) {
		RSS <- (n - sum(GY * as.vector((CC %^% (-1)) %*% GY)))
	} else { RSS <- (n - GY * GY / CC) }
	Fstat <- ((n - m) / m) * (n - RSS) / RSS    # F-statistic
	p <- pf(Fstat, m, n - m, lower.tail = FALSE)
	minP <- 100;
	minM <- 0
	if (p < minP) minM <- M
	minP <- min(p, minP)
	c(minP, minM, CPV[minM])
}

PC <- function(X) {
	eX1 <- eigen(X, symmetric = TRUE)
	#values <- eX1$values
	#vectors <- eX1$vectors
	with (eX1, {
		values[values < 0] <- 0
		#scorenames <- paste('PC', 1:ncol(X), sep = '')
		#colnames(vectors) <- scorenames;
		#rownames(vectors) <- colnames(X)
		prop.var <- values / sum(values)
		cum.var <- cumsum(prop.var)
		list(scores = vectors, importance = cum.var)
	})
}

'PCA' <- function (score.file, gene.file, genes = 'all', cor.path = 'cor/',
anno.type = '', n, beta.par = c(1, 1), weights.function = NULL, user.weights = FALSE,
reference.matrix = TRUE, fun = 'LH', var.fraction = 0.85, write.file = FALSE, quiet = FALSE) {

	do.call(PCA.int, c(as.list(environment()), prob = NA, phred = NA))

}

PCA.int <- function (score.file, gene.file, genes = 'all', cor.path = 'cor/',
anno.type = '', n, beta.par = c(1, 1), weights.function = NULL, user.weights = FALSE,
reference.matrix = TRUE, fun = 'LH', var.fraction = 0.85, write.file = FALSE, quiet = FALSE, prob, phred) {

############ COMMON CHECKS

tmp <- check.input(score.file, cor.path, gene.file, genes)
for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])

############ SPECIFIC CHECKS

tmp <- check.weights(weights.function, beta.par)
for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])

if (!missing(n)) n <- n - 1
if (!is.na(prob)) user.weights <- prob
check.list <- get.check.list('PCA', score.file, anno.type, user.weights, gen.var.weights, fweights, n = n)

############ ANALYSIS

genewise(score.file, gene.file, gf, anno.type, cor.path, cor.file.ext, check.list, write.file, obj0 = list(var.fraction = var.fraction),
	ncl = 5, c('ncomponents', 'explained.variance.fraction'), gen.var.weights, fweights, reference.matrix, fun, n, quiet = quiet, phred = phred, test = 'PCA')

}
