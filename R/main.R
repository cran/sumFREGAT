# sumFREGAT (2017) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

'SKAT' <- function (scoreFile, geneFile, regions, cor.path = '',
annoType = '', beta.par = c(1, 25), weights.function = ifelse(maf > 0, 
dbeta(maf, beta.par[1], beta.par[2]), 0), method = 'kuonen', acc = 1e-8,
lim = 1e+6, rho = FALSE, write.file = FALSE) {

	do.call(do.sumstat, c(make.call(match.call()), test = 'famSKAT'), envir = parent.frame())

}

'BT' <- function (scoreFile, geneFile, regions, cor.path = '',
annoType = '', beta.par = c(1, 25), weights.function = ifelse(maf > 0, 
dbeta(maf, beta.par[1], beta.par[2]), 0), write.file = FALSE) {

	do.call(do.sumstat, c(make.call(match.call()), test = 'famBT'), envir = parent.frame())

}

'MLR' <- function (scoreFile, geneFile, regions, cor.path = '',
annoType = '', n, write.file = FALSE) {

	do.call(do.sumstat, c(make.call(match.call()), test = 'MLR'), envir = parent.frame())

}

'FLM' <- function (scoreFile, geneFile, regions, cor.path = '',
annoType = '', n, beta.par = c(1, 1), weights.function = ifelse(maf > 0,
dbeta(maf, beta.par[1], beta.par[2]), 0), GVF = FALSE, BSF = 'fourier', kg = 30,
kb = 25, order = 4, flip.genotypes = FALSE, Fan = TRUE, write.file = FALSE) {

	do.call(do.sumstat, c(make.call(match.call()), test = 'famFLM'), envir = parent.frame())

}

'PCA' <- function (scoreFile, geneFile, regions, cor.path = '',
annoType = '', n, beta.par = c(1, 1), weights.function = ifelse(maf > 0, 
dbeta(maf, beta.par[1], beta.par[2]), 0), var.fraction = .85, write.file = FALSE) {

	do.call(do.sumstat, c(make.call(match.call()), test = 'PCA'), envir = parent.frame())

}

make.call <- function (cl) {
	cl <- cl0 <- as.list(cl)[-1L]
	names(cl) <- match.arg(names(cl0), names(formals(do.sumstat)), several.ok = TRUE)
	names(cl)[is.na(names(cl))] <- names(cl0)[is.na(names(cl))]
	cl
}

if (getRversion() >= "2.15.1") utils::globalVariables(c('r', 'lgt', 'test', 'rho', 'sumstat.region',
'beta.par', 'kb', 'kg', 'n', 'geneFile', 'annoType', 'nreg', 'flip.genotypes',
'model0', 'g', 'kg0', 'genobasis0', 'J0', 'GVF', 'order0', 'kb0', 'betabasis0', 'BSF', 'stat',
'method', 'acc', 'lim', 'rhos', 'Fan', 'cor.path', 'fweights', 'isDirWritable', 'maf',
'scoreFile', 'var.fraction'))
