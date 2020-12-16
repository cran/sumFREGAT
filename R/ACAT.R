# ACAT() function written by Yaowu Liu

sumstat.ACAT <- function(obj){
	with(obj, with(df, {# Z, w
	#browser()
	p <- pnorm(abs(Z), lower = FALSE) * 2
	p[p == 1] <- 0.999
#	v <- p != 1
#	p <- p[v]
#	w <- w[v]
	
	w <- w / sum(w)
#browser()
    #### check if there are very small non-zero p values
    is.small <- (p < 1e-16)
    if (sum(is.small) == 0) {
        cct.stat <- sum(w * tan((0.5 - p) * pi))
    } else {
        cct.stat <- sum((w[is.small] / p[is.small]) / pi)
        cct.stat <- cct.stat + sum(w[!is.small] * tan((0.5 - p[!is.small]) * pi))
    }
    #### check if the test statistic is very large.
    if (cct.stat > 1e+15){
        pval <- (1 / cct.stat) / pi
    } else {
        pval <- 1 - pcauchy(cct.stat)
    }
    return(pval)
    }))
}


'ACAT' <- function (score.file, gene.file, genes = 'all', anno.type = '',
beta.par = c(1, 1), weights.function = NULL, user.weights = FALSE, gen.var.weights = 'none',
write.file = FALSE, quiet = FALSE) {

############ COMMON CHECKS

tmp <- check.input(score.file, 'do not check cor.path', gene.file, genes)
for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])

############ SPECIFIC CHECKS

tmp <- check.weights(weights.function, beta.par, gen.var.weights)
for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])

check.list <- get.check.list('ACAT', score.file, anno.type, user.weights, gen.var.weights, fweights)

############ ANALYSIS

genewise(score.file, gene.file, gf, anno.type, check.list = check.list, write.file = write.file, gen.var.weights = gen.var.weights, fweights = fweights, quiet = quiet, test = 'ACAT')

}