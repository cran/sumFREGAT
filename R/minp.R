# a wrapper to call GBJ minP()

sumstat.minp <- function(obj) {
	p <- GBJ::minP(test_stats = obj$df$Z, cor_mat = obj$U)$minP_pvalue
	return(p)
}

'minp' <- function (score.file, gene.file, genes = 'all', cor.path = 'cor/',
anno.type = '', write.file = FALSE) {

############ COMMON CHECKS

tmp <- check.input(score.file, cor.path, gene.file, genes)
for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])

############ SPECIFIC CHECKS

check.list <- get.check.list('minP', score.file, anno.type)

############ ANALYSIS

genewise(score.file, gene.file, gf, anno.type, cor.path, cor.file.ext, check.list, write.file, NULL, ncl = 3, test = 'minp')

}