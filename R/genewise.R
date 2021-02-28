# sumFREGAT (2017-2018) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

genewise <- function(score.file, gene.file, gf, anno.type, cor.path, cor.file.ext, check.list, write.file, obj0 = NULL, ncl = 3, cn = NULL,
	gen.var.weights = FALSE, fweights = NULL, reference.matrix = FALSE, fun, n = NULL, Fan = FALSE, flip.genotypes = FALSE, rho, quiet, phred, test) {

	ngenes <- dim(gf)[1]
	gf <- cbind(gf, matrix(NA, nrow = ngenes, ncol = ncl))
	colnames(gf) <- c('gene', 'chrom', 'start', 'end', 'markers', 'filtered.markers', 'pvalue', cn)
	if (write.file != FALSE) write.table(t(colnames(gf)), file = write.file, col.names = FALSE, row.names = FALSE, quote = FALSE)
	nc <- dim(gf)[2]
	nc2 <- nc - 6

	sumstat.region <- as.function(get(paste('sumstat', test, sep = '.')))
	if (test == 'SKAT') {
		if (rho) sumstat.region <- sumstat.SKATO
	}

	for (i in 1:ngenes) {
		gene <- as.character(gf[i, 1])
		obj <- c(obj0, get.sumstat(score.file, gene.file, gene, anno.type, cor.path, cor.file.ext, check.list,
			reference.matrix, gen.var.weights, fweights, n, Fan, flip.genotypes, fun, quiet, phred, test))
			m1 <- check.sumstat(obj, nc2, test)
		if (length(m1) == 1) {
			obj$m <- m1
			gf[i, 5:nc] <- c(obj$m0, m1, sumstat.region(obj))
		} else {
			gf[i, 5:nc] <- m1
		}
		if (write.file != FALSE) write.table(gf[i, ], file = write.file, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
	}

	gf[, 3:7] <- sapply(3:7, function(x) as.numeric(as.character(gf[, x])))
	if (test %in% c('BT', 'PCA')) gf[, 8:9] <- sapply(8:9, function(x) as.numeric(as.character(gf[, x])))
	gf

}
