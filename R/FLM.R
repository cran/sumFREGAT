# sumFREGAT (2017-2018) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

sumstat.FLM <- function(obj) {

	with(obj, with(df, { # Z, POS, w, U, n, m, k, basis, model

	if (m <= k) return(c(sumstat.MLR(obj), 'MLR'))

		B <- eval.basis(POS, basis)  # as matrix (m x Kb)
		WB <- B * w
		Mat <- t(WB) %*% U %*% WB
		BUB_1 <- Mat %^% (-1)

		BZstat <- as.vector(t(WB) %*% Z)
		RSS <- n - sum(BZstat * (BUB_1 %*% BZstat))

		Fstat <- ((n - k) / k) * (n - RSS) / RSS   # F-statistic
		p <- pf(Fstat, k, n - k, lower.tail = FALSE)
		

		return(c(p, model))

	}))

}


'FLM' <- function (score.file, gene.file, genes = 'all', cor.path = 'cor/',
anno.type = '', n, beta.par = c(1, 1), weights.function = NULL, user.weights = FALSE,
basis.function = 'fourier', k = 25, order = 4, flip.genotypes = FALSE,
Fan = TRUE, reference.matrix = TRUE, fun = 'LH', write.file = FALSE, quiet = FALSE) {

############ COMMON CHECKS

tmp <- check.input(score.file, cor.path, gene.file, genes)
for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])

############ SPECIFIC CHECKS

tmp <- check.weights(weights.function, beta.par)
for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])

tmp <- check.spec.FLM(basis.function, k, order)
for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])

if (!missing(n)) n <- n - 1
check.list <- get.check.list('FLM', score.file, anno.type, user.weights, gen.var.weights, fweights, n = n)

############ ANALYSIS

obj0 <- sapply(c('k', 'basis', 'model'),
	function(x) get(x), simplify = FALSE, USE.NAMES = TRUE)

genewise(score.file, gene.file, gf, anno.type, cor.path, cor.file.ext, check.list, write.file, obj0, ncl = 4, 'model',
	gen.var.weights, fweights, reference.matrix, fun, n, Fan, flip.genotypes, quiet = quiet, test = 'FLM')

}