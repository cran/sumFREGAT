# sumFREGAT (2017-2018) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

sumstat.MLR <- function(obj) {
	with(obj, with(df, { # Z, U, n, m
		U05 = U %^% (-.5)
		U05Z = U05 %*% Z        # U^-.5 %*% Z
		R2 <- sum(U05Z^2) / n       # assimp m/n
		Fstat = ((n - m) / m) * R2 / (1-R2)  # F-statistic
		p = as.double(pf(Fstat, m, n - m, lower.tail = FALSE)) # F-test
		return(p)
	}))
}


'MLR' <- function (score.file, gene.file, genes = 'all', cor.path = 'cor/',
anno.type = '', n, reference.matrix = TRUE, fun = 'LH', write.file = FALSE, quiet = FALSE) {

############ COMMON CHECKS

tmp <- check.input(score.file, cor.path, gene.file, genes)
for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])

############ SPECIFIC CHECKS

if (!missing(n)) n <- n - 1
check.list <- get.check.list('MLR', score.file, anno.type, n = n)

############ ANALYSIS

genewise(score.file, gene.file, gf, anno.type, cor.path, cor.file.ext, check.list, write.file, reference.matrix = reference.matrix, fun = fun, n = n, quiet = quiet, test = 'MLR')

}