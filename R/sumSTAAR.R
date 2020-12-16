# sumFREGAT (2020) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

sumSTAAR <- function(score.file, gene.file, genes = 'all', cor.path = 'cor/', tests = c('BT', 'SKAT', 'ACAT'), annotation_phred = paste0('ANNO', 1:10), beta.par.matrix = rbind(c(1, 1), c(1, 25)), n = NA, write.file = FALSE, quiet = FALSE) {

if (any(c('PCA', 'FLM') %in% tests) & is.na(n)) stop('n must be set for PCA/FLM analyses') 
if ('MLR' %in% tests) {
	warning("Weights do not change MLR statistics, 'MLR' excluded from 'tests'")
	tests <- tests[tests != 'MLR']
}

beta.i <- dim(beta.par.matrix)[2]

pval.all <- c()

for (tt in tests) {

	sumstat.function <- as.function(get(tt))
	my.args0 <- list(score.file = score.file, gene.file = gene.file, genes = genes, quiet = quiet)
	if (tt == 'ACAT') {
		my.args0 <- c(my.args0, gen.var.weights = 'af')
	} else {
		my.args0 <- c(my.args0, cor.path = cor.path)
	}
	if (tt %in% c('PCA', 'FLM')) my.args0 <- c(my.args0, n = n)

	for (i in 1:beta.i) {
		a1 <- beta.par.matrix[i, 1]
		a2 <- beta.par.matrix[i, 2]
		my.args.beta <- my.args0
		my.args.beta$beta.par <- c(a1, a2)
		pval <- c()
		ncyc <- ifelse(is.na(annotation_phred[1]), 1, length(annotation_phred) + 1)
		for (a in 1:ncyc) {
			uw <- ifelse(a == 1, FALSE, annotation_phred[a - 1])
			wf <- ifelse(write.file != FALSE, paste(tt, a1, a2, ifelse(a == 1, 'ANNO0', annotation_phred[a - 1]), write.file, sep = '.'), FALSE)
			my.args <- c(my.args.beta, user.weights = uw, write.file = wf)
			res <- do.call('sumstat.function', my.args)
			pval[a] <- res$pvalue
		}
		#browser()
		pval <- c(pval, ACATO(pval))
		pval <- t(pval)
		if (is.na(annotation_phred[1])) {
			names.tmp <- c('ANNO0', 'STAAR')
		} else {
			names.tmp <- c('ANNO0', annotation_phred, 'STAAR')
		}
		colnames(pval) <- paste(tt, a1, a2, names.tmp, sep = '.')
		pval.all <- cbind(pval.all, pval)
	}
}
#browser()
v <- grepl('ANNO0', colnames(pval.all))
pval.all <- cbind(pval.all, sumSTAAR.ACAT_O = sapply(1:dim(pval.all)[1], function(x) ACATO(pval.all[x, v])))
v <- !grepl('STAAR', colnames(pval.all))
pval.all <- cbind(pval.all, sumSTAAR.STAAR_O = sapply(1:dim(pval.all)[1], function(x) ACATO(pval.all[x, v])))
pval.all <- cbind(gene = res$gene, pval.all)

if (write.file != FALSE) write.table(pval.all, file = write.file, row.names = FALSE, quote = FALSE)

as.data.frame(pval.all)

}
