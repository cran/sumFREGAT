# sumFREGAT (2017-2018) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

sumstat.SKAT <- function(obj) {

	with(obj, with(df, {# Z, U, w, method, acc, lim

		Q <- sum((w * Z) ^ 2) # statistic
		KKK <- t(U * w) * w # kernel matrix
		if (method != 'hybrid') {
			eig <- eigen(KKK, symmetric = TRUE, only.values = TRUE)
			ev <- eig$values[eig$values > 1e-6 * eig$values[1]]
			if (method == 'kuonen') {
				p <- pchisqsum(Q, rep(1, length(ev)), ev, lower.tail = F, method = 'sad') 
			} else { #if (method == 'davies') {
				p <- davies(Q, ev, acc = acc, lim = lim)$Qq
			} 
		} else { #if (method == 'hybrid') {
			lam <- svd(KKK, nu = 0, nv = 0)$d 	#lam <- lam[lam > 1e-6 * lam[1]]
			p <- KAT.pval(Q, lam)
		}
		return(max(min(p, 1), 0))

	}))

}


'SKAT' <- function (score.file, gene.file, genes = 'all', cor.path = 'cor/', anno.type = '', beta.par = c(1, 25), weights.function = NULL,
user.weights = FALSE, gen.var.weights = 'se.beta', method = 'kuonen', acc = 1e-8, lim = 1e+6, rho = FALSE,
p.threshold = 0.8, write.file = FALSE) {

############ COMMON CHECKS

tmp <- check.input(score.file, cor.path, gene.file, genes)
for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])

############ SPECIFIC CHECKS

tmp <- check.weights(weights.function, beta.par, gen.var.weights)
for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])

tmp <- check.spec.SKAT(method, rho)
for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])

check.list <- get.check.list('SKAT', score.file, anno.type, user.weights, gen.var.weights, fweights, rho)

############ ANALYSIS

obj0 <- sapply(c('method', 'acc', 'lim', 'rhos', 'p.threshold'),
	function(x) get(x), simplify = FALSE, USE.NAMES = TRUE)

genewise(score.file, gene.file, gf, anno.type, cor.path, cor.file.ext, check.list, write.file, obj0, ncl = 3, NULL, gen.var.weights, fweights, rho = rho, test = 'SKAT')

}


